# -*- coding: utf-8 -*-
"""把 Mermaid 流程图（graph/flowchart + subgraph + 箭头链）渲染为 PNG。

支持子集：
  - graph TB/TD/BT/LR/RL  或  flowchart ...
  - subgraph ID["标题"] ... end  （含内部 direction LR/TB）
  - 节点定义：A["标签"]  A("标签")  A[[标签]]  A{标签}  A(不带引号)
  - 边：A --> B --> C  /  A -->|标签| B  /  A --- B  /  A -.-> B  /  A ==> B
  - 不支持的（sequence/class 等）返回 None，由上层回退占位图。

布局策略：根方向决定各 block（顶层节点组 / subgraph）的排布方向；
subgraph 内部按自身 direction 排成一行/一列并画容器边框与标题；
箭头方向由方向符与边类型共同决定，箭头头指向目标端。
"""
import io
import os
import re
import base64
from PIL import Image, ImageDraw, ImageFont

FONT_R = "C:/Windows/Fonts/msyh.ttc"
FONT_B = "C:/Windows/Fonts/msyhbd.ttc"

PALETTE = [
    ((238, 244, 251), (91, 141, 239)),   # 蓝
    ((228, 245, 235), (46, 160, 120)),   # 绿
    ((253, 241, 234), (224, 122, 80)),   # 橙
    ((245, 238, 253), (150, 120, 200)),  # 紫
]
CONTAINER_FILL = (250, 251, 253)
CONTAINER_BORDER = (170, 185, 205)
TITLE_COLOR = (60, 66, 82)
INK = (60, 66, 82)
ARROW = (120, 132, 152)
DOTTED = (150, 160, 178)

GAP_H = 50      # 同排节点水平间距（含箭头空间）
GAP_V = 42      # 同列节点垂直间距
PAD = 18        # subgraph 容器内边距
TITLE_H = 30    # subgraph 标题高度
MARGIN = 26     # 整图外边距
BLOCK_GAP = 46  # block 之间间距
NODE_MAXW = 168 # 节点最大宽度
NODE_FS = 14
TITLE_FS = 15


def font(size, bold=False):
    try:
        return ImageFont.truetype(FONT_B if bold else FONT_R, size)
    except Exception:
        return ImageFont.load_default()


def _is_word(ch):
    c = ord(ch)
    return (65 <= c <= 90) or (97 <= c <= 122) or (48 <= c <= 57)


def wrap(draw, text, f, max_w):
    """按 max_w 换行；优先在英文/数字词之前断开，避免 SiteMap 这类词被逐字母截断。"""
    lines, cur = [], ""
    for ch in text:
        if draw.textlength(cur + ch, font=f) <= max_w:
            cur += ch
            continue
        # 溢出：优先把末尾的 ASCII 词（含数字）推到下一行
        split = None
        if cur:
            i = len(cur) - 1
            while i >= 0 and _is_word(cur[i]):
                i -= 1
            if i + 1 < len(cur):
                split = i + 1  # ASCII 词起点
        if split is not None and split > 0:
            lines.append(cur[:split])
            cur = cur[split:] + ch
        else:
            if cur:
                lines.append(cur)
            cur = ch
    if cur:
        lines.append(cur)
    return lines


def node_box(draw, label, fs=NODE_FS, max_w=NODE_MAXW):
    f = font(fs)
    lines = wrap(draw, label, f, max_w - 22)
    line_h = fs + 6
    w = max((draw.textlength(ln, font=f) for ln in lines), default=10) + 22
    w = max(w, 60)
    h = len(lines) * line_h + 14
    return int(w), int(h), lines, line_h, f


# ===================== 解析 =====================
def _parse_node_spec(spec):
    spec = spec.strip()
    m = re.match(r'^(\w+)\s*\[\s*"([^"]*)"\s*\]$', spec)
    if m:
        return m.group(1), m.group(2)
    m = re.match(r"^(\w+)\s*\[\s*'([^']*)'\s*\]$", spec)
    if m:
        return m.group(1), m.group(2)
    m = re.match(r'^(\w+)\s*\[\s*([^\]]*?)\s*\]$', spec)
    if m:
        return m.group(1), m.group(2)
    m = re.match(r'^(\w+)\s*\(\s*"([^"]*)"\s*\)$', spec)
    if m:
        return m.group(1), m.group(2)
    m = re.match(r'^(\w+)\s*\(\s*([^)]*?)\s*\)$', spec)
    if m:
        return m.group(1), m.group(2)
    m = re.match(r'^(\w+)\s*\{([^}]*)\}$', spec)
    if m:
        return m.group(1), m.group(2)
    m = re.match(r'^(\w+)$', spec)
    if m:
        return m.group(1), m.group(1)
    return None, None


def _split_line(line):
    """把一行拆成 (节点spec片段列表, 边列表)。边: (kind, label)。"""
    # 链接 token
    link_re = re.compile(r'(-->|---|-\.->|==>)(?:\s*\|\s*([^|]*?)\s*\|)?')
    parts = link_re.split(line)
    # parts: [spec0, kind1, label1, spec1, kind2, label2, spec2, ...]
    specs = [parts[0]] + parts[3::3]
    edges = []
    i = 1
    while i < len(parts) and parts[i] is not None:
        kind = parts[i]
        label = parts[i + 1]
        edges.append((kind, label))
        i += 3
    return [s.strip() for s in specs if s.strip()], edges


def parse_mermaid(content):
    lines = content.splitlines()
    root_dir = "TB"
    nodes = {}          # id -> label
    subgraphs = []      # list of dict
    root_nodes = []     # 顶层节点 id（有序）
    root_edges = []     # (a, b, kind, label)
    sg_ids = set()

    cur_sg = None

    for raw in lines:
        ln = raw.strip()
        if not ln or ln.startswith("%%"):
            continue
        if re.match(r'^(graph|flowchart)\b', ln) and cur_sg is None:
            m = re.match(r'^(?:graph|flowchart)\s+(TB|TD|BT|LR|RL)\b', ln)
            if m:
                root_dir = m.group(1)
            continue
        m = re.match(r'^subgraph\s+(.*)$', ln)
        if m and cur_sg is None:
            hdr = m.group(1).strip()
            sg_id, sg_title = _parse_subgraph_header(hdr)
            cur_sg = {"id": sg_id, "title": sg_title, "dir": None,
                      "node_order": [], "edges": [], "node_labels": {}}
            sg_ids.add(sg_id)
            subgraphs.append(cur_sg)
            continue
        if ln == "end" and cur_sg is not None:
            cur_sg = None
            continue
        if cur_sg is not None:
            if re.match(r'^direction\s+(LR|TB|RL|BT)$', ln):
                cur_sg["dir"] = re.match(r'^direction\s+(LR|TB|RL|BT)$', ln).group(1)
                continue
            specs, edges = _split_line(ln)
            for sp in specs:
                nid, lab = _parse_node_spec(sp)
                if nid:
                    if nid not in cur_sg["node_labels"]:
                        cur_sg["node_labels"][nid] = lab
                    if nid not in cur_sg["node_order"]:
                        cur_sg["node_order"].append(nid)
            for k, (kind, lab) in enumerate(edges):
                a = _parse_node_spec(specs[k])[0]
                b = _parse_node_spec(specs[k + 1])[0]
                if a and b:
                    cur_sg["edges"].append((a, b, kind, lab))
            continue
        # 顶层
        specs, edges = _split_line(ln)
        for sp in specs:
            nid, lab = _parse_node_spec(sp)
            if nid and nid not in sg_ids:
                if nid not in nodes:
                    nodes[nid] = lab
                if nid not in root_nodes:
                    root_nodes.append(nid)
        for k, (kind, lab) in enumerate(edges):
            a = _parse_node_spec(specs[k])[0]
            b = _parse_node_spec(specs[k + 1])[0]
            if a and b:
                root_edges.append((a, b, kind, lab))

    return {
        "root_dir": root_dir,
        "nodes": nodes,
        "subgraphs": subgraphs,
        "root_nodes": root_nodes,
        "root_edges": root_edges,
        "sg_ids": sg_ids,
    }


def _parse_subgraph_header(hdr):
    # subgraph S1["1.诱导采样"]  /  subgraph ["标题"]  /  subgraph 标题
    m = re.match(r'^(\w+)\s*\[\s*"([^"]*)"\s*\]$', hdr)
    if m:
        return m.group(1), m.group(2)
    m = re.match(r'^(\w+)\s*\[\s*([^\]]*?)\s*\]$', hdr)
    if m:
        return m.group(1), m.group(2)
    m = re.match(r'^\[\s*"([^"]*)"\s*\]$', hdr)
    if m:
        return "_sg" + str(id(hdr)), m.group(1)
    # 纯标题
    return "_sg" + str(id(hdr)), hdr


# ===================== 布局 + 绘制 =====================
class Renderer:
    def __init__(self, parsed):
        self.p = parsed
        self.img = None
        self.d = None

    # 在一个 group 内，按 dir 布局节点，返回 boxes(id->[x,y,w,h]) 局部坐标(>=0)
    def _layout_nodes(self, ids, labels, draw, dirn):
        boxes = {}
        if dirn in ("LR", "RL"):
            x, y = 0, 0
            maxh = 0
            order = ids if dirn == "LR" else list(reversed(ids))
            widths = {}
            for nid in order:
                w, h, _, _, _ = node_box(draw, labels.get(nid, nid))
                widths[nid] = (w, h)
                maxh = max(maxh, h)
            cx = 0
            for nid in order:
                w, h = widths[nid]
                boxes[nid] = [cx, (maxh - h) // 2, w, h]
                cx += w + GAP_H
            total_w = cx - GAP_H
            return boxes, total_w, maxh
        else:  # TB / BT
            x, y = 0, 0
            maxw = 0
            order = ids if dirn == "TB" else list(reversed(ids))
            heights = {}
            for nid in order:
                w, h, _, _, _ = node_box(draw, labels.get(nid, nid))
                heights[nid] = (w, h)
                maxw = max(maxw, w)
            cy = 0
            for nid in order:
                w, h = heights[nid]
                boxes[nid] = [(maxw - w) // 2, cy, w, h]
                cy += h + GAP_V
            total_h = cy - GAP_V
            return boxes, maxw, total_h

    def _draw_node(self, draw, x, y, w, h, label, color_idx):
        fill, border = PALETTE[color_idx % len(PALETTE)]
        draw.rounded_rectangle([x, y, x + w, y + h], radius=9,
                               fill=fill, outline=border, width=2)
        f = font(NODE_FS)
        # 用布局时的统一 max_w，避免重算换行导致文字宽度超过框体
        lines = wrap(draw, label, f, NODE_MAXW - 22)
        lh = NODE_FS + 6
        th = len(lines) * lh
        ty = y + (h - th) / 2
        for li, line in enumerate(lines):
            draw.text((x + w / 2, ty + li * lh), line, font=f,
                      fill=INK, anchor="mm")

    def _draw_edge(self, draw, a, b, kind, label, boxes, dirn, color):
        ax, ay, aw, ah = boxes[a]
        bx, by, bw, bh = boxes[b]
        if dirn in ("LR", "RL"):
            if dirn == "LR":
                x1, y1 = ax + aw, ay + ah / 2
                x2, y2 = bx, by + bh / 2
            else:
                x1, y1 = ax, ay + ah / 2
                x2, y2 = bx + bw, by + bh / 2
            self._line(draw, x1, y1, x2, y2, kind, color)
            self._arrow(draw, x2, y2, "right" if x2 > x1 else "left", color)
            if label:
                mx, my = (x1 + x2) / 2, (y1 + y2) / 2 - 10
                draw.text((mx, my), label, font=font(12), fill=ARROW, anchor="mm")
        else:
            if dirn == "TB":
                x1, y1 = ax + aw / 2, ay + ah
                x2, y2 = bx + bw / 2, by
            else:
                x1, y1 = ax + aw / 2, ay
                x2, y2 = bx + bw / 2, by + bh
            self._line(draw, x1, y1, x2, y2, kind, color)
            self._arrow(draw, x2, y2, "down" if y2 > y1 else "up", color)
            if label:
                mx, my = (x1 + x2) / 2 + 8, (y1 + y2) / 2
                draw.text((mx, my), label, font=font(12), fill=ARROW, anchor="mm")

    @staticmethod
    def _line(draw, x1, y1, x2, y2, kind, color):
        if kind == "---":
            draw.line([x1, y1, x2, y2], fill=DOTTED, width=2)
        elif kind == "-.->":
            draw.line([x1, y1, x2, y2], fill=DOTTED, width=2, dash=(6, 5))
        elif kind == "==>":
            draw.line([x1, y1, x2, y2], fill=color, width=4)
        else:
            draw.line([x1, y1, x2, y2], fill=color, width=2)

    @staticmethod
    def _arrow(draw, x, y, direction, color):
        s = 7
        if direction == "right":
            draw.polygon([(x, y), (x - s, y - s), (x - s, y + s)], fill=color)
        elif direction == "left":
            draw.polygon([(x, y), (x + s, y - s), (x + s, y + s)], fill=color)
        elif direction == "down":
            draw.polygon([(x, y), (x - s, y - s), (x + s, y - s)], fill=color)
        else:  # up
            draw.polygon([(x, y), (x - s, y + s), (x + s, y + s)], fill=color)

    def _draw_group(self, draw, origin_x, origin_y, ids, labels, dirn, color_offset=0):
        boxes, _, _ = self._layout_nodes(ids, labels, draw, dirn)
        # 平移到 origin
        for nid in boxes:
            boxes[nid][0] += origin_x
            boxes[nid][1] += origin_y
        for idx, nid in enumerate(ids):
            x, y, w, h = boxes[nid]
            self._draw_node(draw, x, y, w, h, labels.get(nid, nid),
                            color_offset + idx)
        # 组内边（顶层节点间少见，仍处理）
        return boxes

    def render(self):
        p = self.p
        draw = ImageDraw.Draw(Image.new("RGB", (10, 10)))

        # 预计算各 block 尺寸与局部布局
        blocks = []  # dict: kind, id, w, h, draw_fn 元数据
        # 先准备 subgraph 局部布局
        self.sg_layout = {}
        for sg in p["subgraphs"]:
            dirn = sg["dir"] or ("LR" if p["root_dir"] in ("TB", "TD", "BT") else "TB")
            boxes, cw, ch = self._layout_nodes(
                sg["node_order"], sg["node_labels"], draw, dirn)
            # 容器尺寸
            maxx = max((b[0] + b[2] for b in boxes.values()), default=0)
            maxy = max((b[1] + b[3] for b in boxes.values()), default=0)
            cont_w = int(maxx + 2 * PAD)
            cont_h = int(maxy + 2 * PAD + TITLE_H)
            self.sg_layout[sg["id"]] = {
                "dirn": dirn, "boxes": boxes, "edges": sg["edges"],
                "labels": sg["node_labels"], "title": sg["title"],
                "cont_w": cont_w, "cont_h": cont_h,
            }
            blocks.append({"kind": "sg", "id": sg["id"],
                           "w": cont_w, "h": cont_h})

        # 顶层节点组（每个顶层节点各自作为一个小 block；若它们之间在 root_edges 有连接则同组）
        # 简化：把顶层节点各自成 block，按 root_nodes 顺序
        self.top_layout = {}
        for nid in p["root_nodes"]:
            boxes, w, h = self._layout_nodes([nid], p["nodes"], draw, "LR")
            self.top_layout[nid] = {"boxes": boxes, "w": w, "h": h}
            blocks.append({"kind": "top", "id": nid, "w": int(w), "h": int(h)})

        # 决定 block 顺序（按 root_edges 拓扑，失败则声明顺序）
        order = self._block_order(blocks)
        blocks_ordered = sorted(
            [b for b in blocks if b["id"] in order],
            key=lambda b: order.index(b["id"]))

        rd = p["root_dir"]
        vertical = rd in ("TB", "TD", "BT")
        # 计算画布尺寸
        if vertical:
            W = max((b["w"] for b in blocks_ordered), default=0) + 2 * MARGIN
            H = sum(b["h"] for b in blocks_ordered) + \
                BLOCK_GAP * (len(blocks_ordered) - 1) + 2 * MARGIN
        else:
            W = sum(b["w"] for b in blocks_ordered) + \
                BLOCK_GAP * (len(blocks_ordered) - 1) + 2 * MARGIN
            H = max((b["h"] for b in blocks_ordered), default=0) + 2 * MARGIN

        self.img = Image.new("RGB", (int(W), int(H)), "white")
        d = ImageDraw.Draw(self.img)
        self.d = d

        # 放置各 block
        cx = W / 2
        cy = H / 2
        cur = MARGIN
        placed = {}  # id -> (ox, oy, w, h)
        for b in blocks_ordered:
            if vertical:
                ox = cx - b["w"] / 2
                oy = cur
                cur += b["h"] + BLOCK_GAP
            else:
                ox = cur
                oy = cy - b["h"] / 2
                cur += b["w"] + BLOCK_GAP
            placed[b["id"]] = (ox, oy, b["w"], b["h"])
            self._draw_block(d, b, ox, oy)

        # 画 block 间（root）箭头
        for k in range(len(blocks_ordered) - 1):
            a = blocks_ordered[k]["id"]
            b = blocks_ordered[k + 1]["id"]
            ax, ay, aw, ah = placed[a]
            bx, by, bw, bh = placed[b]
            if vertical:
                x1, y1 = ax + aw / 2, ay + ah
                x2, y2 = bx + bw / 2, by
                self._line(d, x1, y1, x2, y2, "-->", ARROW)
                self._arrow(d, x2, y2, "down", ARROW)
            else:
                x1, y1 = ax + aw, ay + ah / 2
                x2, y2 = bx, by + bh / 2
                self._line(d, x1, y1, x2, y2, "-->", ARROW)
                self._arrow(d, x2, y2, "right", ARROW)

        return self.img

    def _block_order(self, blocks):
        # 用 root_edges 拓扑排序，端点可能是 root_nodes 或 subgraph id
        ids = [b["id"] for b in blocks]
        succ = {i: [] for i in ids}
        indeg = {i: 0 for i in ids}
        idset = set(ids)
        for a, b, _, _ in self.p["root_edges"]:
            if a in idset and b in idset and a != b:
                if b not in succ[a]:
                    succ[a].append(b)
                    indeg[b] += 1
        # Kahn
        from collections import deque
        q = deque([i for i in ids if indeg[i] == 0])
        order = []
        while q:
            n = q.popleft()
            order.append(n)
            for m in succ[n]:
                indeg[m] -= 1
                if indeg[m] == 0:
                    q.append(m)
        if len(order) == len(ids):
            return order
        return ids  # 有环则退回声明顺序

    def _draw_block(self, d, b, ox, oy):
        if b["kind"] == "sg":
            L = self.sg_layout[b["id"]]
            # 容器
            d.rounded_rectangle([ox, oy, ox + L["cont_w"], oy + L["cont_h"]],
                               radius=12, fill=CONTAINER_FILL,
                               outline=CONTAINER_BORDER, width=2)
            d.text((ox + 12, oy + 7), L["title"], font=font(TITLE_FS, True),
                   fill=TITLE_COLOR)
            # 节点（平移进容器）
            shift_x = ox + PAD
            shift_y = oy + PAD + TITLE_H
            boxes = {}
            for nid, bb in L["boxes"].items():
                nb = [bb[0] + shift_x, bb[1] + shift_y, bb[2], bb[3]]
                boxes[nid] = nb
            # 画边
            for a, b2, kind, lab in L["edges"]:
                if a in boxes and b2 in boxes:
                    self._draw_edge(d, a, b2, kind, lab, boxes, L["dirn"], ARROW)
            # 画节点（在上层）
            for idx, nid in enumerate(L["boxes"].keys()):
                x, y, w, h = boxes[nid]
                self._draw_node(d, x, y, w, h, L["labels"].get(nid, nid), idx)
        else:
            # 顶层节点
            L = self.top_layout[b["id"]]
            boxes = {}
            for nid, bb in L["boxes"].items():
                nb = [bb[0] + ox, bb[1] + oy, bb[2], bb[3]]
                boxes[nid] = nb
            for nid in L["boxes"]:
                x, y, w, h = boxes[nid]
                self._draw_node(d, x, y, w, h, self.p["nodes"].get(nid, nid), 0)


def render_to_bytes(content, quality=85):
    """返回 (b64, note) 或 (None, None)。"""
    try:
        parsed = parse_mermaid(content)
    except Exception:
        return None, None
    if not parsed["subgraphs"] and not parsed["root_nodes"]:
        return None, None
    try:
        r = Renderer(parsed)
        img = r.render()
    except Exception:
        return None, None
    buf = io.BytesIO()
    img.save(buf, format="JPEG", quality=quality)
    return base64.b64encode(buf.getvalue()).decode(), ""


if __name__ == "__main__":
    import sys
    src = sys.stdin.read()
    b64, _ = render_to_bytes(src)
    if b64:
        print("OK", len(b64))
    else:
        print("FAIL")
