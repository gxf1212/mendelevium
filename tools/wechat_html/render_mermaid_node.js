// render_mermaid_node.js — 用 mermaid.js 官方引擎（mermaid.ai 同款）渲染 Mermaid 图为 PNG。
// 与内置 Pillow 渲染器的区别：官方引擎原生支持 <br/> 换行、subgraph、各种节点形状、
// 官方主题（default/forest/neutral/base），样式更精致。
//
// 依赖：mermaid + puppeteer-core（在 ~/.workbuddy/binaries/node/workspace 安装，
// puppeteer-core 用系统 Chrome/Edge，不下载 Chromium）。
// 运行需设 NODE_PATH 指向 workspace/node_modules。
//
// 输入(JSON, stdin)：{"items":[{"id":"0","code":"graph TD; A-->B","theme":"default",
//   "primaryColor":"#8a7fd1","primaryBorderColor":"#6f63b8","primaryTextColor":"#3f3f3f",
//   "lineColor":"#8a7fd1","fontFamily":"Microsoft YaHei","fontSize":"15px"}, ...]}
// 输出(JSON, stdout)：{"0":"<base64 png>"|null, ...}

const fs = require('fs');
const path = require('path');
const puppeteer = require('puppeteer-core');

// 系统浏览器候选（按优先级）
const CHROME_CANDIDATES = [
  'C:/Program Files/Google/Chrome/Application/chrome.exe',
  'C:/Program Files (x86)/Google/Chrome/Application/chrome.exe',
  'C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe',
  'C:/Program Files/Microsoft/Edge/Application/msedge.exe',
];

function findBrowser() {
  for (const p of CHROME_CANDIDATES) {
    if (fs.existsSync(p)) return p;
  }
  return null;
}

let MERMAID_PATH = null;
try {
  MERMAID_PATH = require.resolve('mermaid/dist/mermaid.min.js');
} catch (e) {
  // 尝试常见相对位置
  const guesses = [
    path.join(__dirname, 'node_modules', 'mermaid', 'dist', 'mermaid.min.js'),
  ];
  for (const g of guesses) {
    if (fs.existsSync(g)) { MERMAID_PATH = g; break; }
  }
}

function defaultThemeVars(item) {
  return {
    theme: item.theme || 'default',
    themeVariables: {
      primaryColor: item.primaryColor || '#ecf3ff',
      primaryBorderColor: item.primaryBorderColor || '#5b8def',
      primaryTextColor: item.primaryTextColor || '#333333',
      lineColor: item.lineColor || '#5b8def',
      fontFamily: item.fontFamily || '"Microsoft YaHei", "PingFang SC", sans-serif',
      fontSize: item.fontSize || '15px',
      // 让背景透明，截图时由页面背景兜底（白底）
      background: 'transparent',
    },
  };
}

async function renderOne(page, code, item) {
  const cfg = defaultThemeVars(item);
  const result = await page.evaluate(async (code, cfg) => {
    try {
      if (!window.mermaid) return { ok: false, err: 'mermaid not loaded' };
      window.mermaid.initialize(Object.assign({ startOnLoad: false }, cfg));
      const { svg } = await window.mermaid.render('mmd-' + Math.random().toString(36).slice(2), code);
      const holder = document.getElementById('mmd');
      holder.innerHTML = svg;
      return { ok: true };
    } catch (e) {
      return { ok: false, err: String(e && e.message ? e.message : e) };
    }
  }, code, cfg);

  if (!result.ok) return null;

  // 等 SVG 完成布局
  await page.evaluate(() => new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r))));

  // 截图 svg 元素；若 svg 尺寸异常则截图整个 holder
  let buf = null;
  try {
    const el = await page.$('#mmd svg');
    const box = await el.boundingBox();
    if (box && box.width > 2 && box.height > 2) {
      buf = await el.screenshot({ type: 'png' });
    }
  } catch (e) { /* ignore */ }

  if (!buf) {
    try {
      const holder = await page.$('#mmd');
      buf = await holder.screenshot({ type: 'png' });
    } catch (e) { /* ignore */ }
  }
  return buf;
}

async function main() {
  let raw = '';
  process.stdin.on('data', d => (raw += d));
  process.stdin.on('end', async () => {
    let req;
    try { req = JSON.parse(raw); } catch (e) {
      process.stderr.write('bad json\n'); process.exit(1);
    }
    const items = req.items || [];
    const out = {};
    const browserPath = findBrowser();
    if (!browserPath || !MERMAID_PATH) {
      // 无浏览器或 mermaid 未装：全部置 null，交由上层回退
      for (const it of items) out[it.id] = null;
      process.stdout.write(JSON.stringify(out));
      return;
    }

    let browser = null;
    try {
      browser = await puppeteer.launch({
        executablePath: browserPath,
        headless: 'new',
        args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage',
               '--hide-scrollbars', '--font-render-hinting=none'],
      });
      const page = await browser.newPage();
      await page.setViewport({ width: 1600, height: 1200, deviceScaleFactor: 2 });
      await page.setContent(
        '<html><body style="margin:0;padding:8px;background:#ffffff;">' +
        '<div id="mmd"></div></body></html>');
      await page.addScriptTag({ path: MERMAID_PATH });

      for (const it of items) {
        try {
          const buf = await renderOne(page, it.code, it);
          out[it.id] = buf ? buf.toString('base64') : null;
        } catch (e) {
          out[it.id] = null;
        }
      }
    } catch (e) {
      for (const it of items) out[it.id] = null;
    } finally {
      if (browser) { try { await browser.close(); } catch (e) { /* ignore */ } }
    }
    process.stdout.write(JSON.stringify(out));
  });
}

main();
