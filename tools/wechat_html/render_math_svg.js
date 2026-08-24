// render_math_svg.js — 用 MathJax(tex-svg) 把 LaTeX 公式转成微信友好的内联 SVG。
// 关键：fontCache:'none' 让每个字形直接输出为 <path>，不产生 <defs>（微信拒绝带 <defs> 的 svg）。
// 输入(JSON, stdin)：{"items":[{"id":"0","tex":"x^2","display":false}, ...]}
// 输出(JSON, stdout)：{"0":"<svg ...>...</svg>","1":null,...}（失败为 null）

const mathjax = require('mathjax-full/js/mathjax.js').mathjax;
const TeX = require('mathjax-full/js/input/tex.js').TeX;
const SVG = require('mathjax-full/js/output/svg.js').SVG;
const liteAdaptor = require('mathjax-full/js/adaptors/liteAdaptor.js').liteAdaptor;
const RegisterHTMLHandler = require('mathjax-full/js/handlers/html.js').RegisterHTMLHandler;
// 显式加载各包配置，否则 packages 里列了也不会自动注册对应宏：
// 之前漏了 ams，导致 \dfrac / \genfrac / \substack 等报 Undefined control sequence。
require('mathjax-full/js/input/tex/ams/AmsConfiguration.js');
require('mathjax-full/js/input/tex/mhchem/MhchemConfiguration.js');
require('mathjax-full/js/input/tex/boldsymbol/BoldsymbolConfiguration.js');
require('mathjax-full/js/input/tex/newcommand/NewcommandConfiguration.js');

const COLOR = '#33312e'; // 与正文同色的深炭灰，写死（微信会丢 currentColor 继承）

const adaptor = liteAdaptor();
RegisterHTMLHandler(adaptor);

const tex = new TeX({ packages: ['base', 'ams', 'mhchem', 'boldsymbol', 'newcommand'] });
const svgOut = new SVG({ fontCache: 'none' });
const doc = mathjax.document('', { InputJax: tex, OutputJax: svgOut });

function sanitize(svgHtml, display) {
  const m = svgHtml.match(/<svg[\s\S]*?<\/svg>/);
  if (!m) return null;
  let s = m[0];
  // 微信会丢 currentColor 继承，写死颜色
  s = s.replace(/currentColor/g, COLOR);
  // 去掉微信不喜欢的属性
  s = s.replace(/\s(jax|display|tabindex|ctxtmenu_counter)=\"[^\"]*\"/g, '');
  // 取出 MathJax 已有的 style（含 vertical-align 基线修正），随后整体移除原 style 属性
  const styleM = s.match(/<svg[\s\S]*?\sstyle=\"([^\"]*)\"/);
  let baseStyle = styleM ? styleM[1] : '';
  if (baseStyle && !baseStyle.endsWith(';')) baseStyle += ';';
  // 把 width/height 属性并入 style（mdnice 同款，粘贴更稳）
  const wm = s.match(/\bwidth=\"([^\"]+)\"/);
  const hm = s.match(/\bheight=\"([^\"]+)\"/);
  let extra = '';
  if (display) {
    extra += 'display:block;margin:16px auto;max-width:96%;';
    if (wm && hm) extra += `width:${wm[1]};height:${hm[1]};`;
    // 块级公式不需要 inline 的 vertical-align
    baseStyle = baseStyle.replace(/vertical-align:[^;]+;?/g, '');
  } else {
    if (wm && hm) extra += `width:${wm[1]};height:${hm[1]};`;
    // 保留 MathJax 的 vertical-align 基线修正，若无则给 middle
    if (!/vertical-align/.test(baseStyle)) extra += 'vertical-align:middle;';
  }
  const newStyle = (baseStyle + extra).replace(/;+$/, '');
  // 先删掉原来的 style 属性（避免重复 style），再在 <svg> 末尾统一写新 style
  s = s.replace(/\sstyle=\"[^\"]*\"/, '');
  s = s.replace(/<svg([\s\S]*?)>/, `<svg$1 style=\"${newStyle}\">`);
  s = s.replace(/\swidth=\"[^\"]+\"/, '').replace(/\sheight=\"[^\"]+\"/, '');
  return s;
}

function main() {
  let raw = '';
  process.stdin.on('data', d => raw += d);
  process.stdin.on('end', () => {
    let req;
    try { req = JSON.parse(raw); } catch (e) { process.stderr.write('bad json\n'); process.exit(1); }
    const items = req.items || [];
    const out = {};
    for (const it of items) {
      try {
        const node = doc.convert(it.tex, { display: !!it.display });
        const html = adaptor.outerHTML(node);
        out[it.id] = sanitize(html, !!it.display);
      } catch (e) {
        out[it.id] = null;
      }
    }
    process.stdout.write(JSON.stringify(out));
  });
}

main();
