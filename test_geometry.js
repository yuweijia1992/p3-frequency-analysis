/* ============================================================================
 * test_geometry.js —— 绘图几何冒烟测试：用记录型 canvas 桩检查坐标合理性
 * 运行：node test_geometry.js
 * ==========================================================================*/
"use strict";

const fs = require("fs");
const path = require("path");
const { JSDOM } = require("jsdom");

const html = fs.readFileSync(path.join(__dirname, "P3频率曲线网页版.html"), "utf8");

/* 记录型 ctx */
const calls = { path: [], arcs: [], texts: [] };
let currentPath = null;
const ctxRec = {
    set fillStyle(v) {}, get fillStyle() { return "#000"; },
    set strokeStyle(v) {}, get strokeStyle() { return "#000"; },
    set lineWidth(v) {}, get lineWidth() { return 1; },
    set font(v) {}, set textAlign(v) {}, set textBaseline(v) {},
    beginPath() { currentPath = []; },
    moveTo(x, y) { if (currentPath) currentPath.push([x, y]); },
    lineTo(x, y) { if (currentPath) currentPath.push([x, y]); },
    arc(x, y, r) { calls.arcs.push({ x, y, r }); },
    stroke() { if (currentPath) { calls.path.push(currentPath); currentPath = null; } },
    fill() {},
    fillText(t, x, y) { calls.texts.push({ t, x, y }); },
    measureText(t) { return { width: String(t).length * 8 }; },
    strokeRect() {}, fillRect() {}, clearRect() {},
    save() {}, restore() {}, rotate() {}, translate() {}, scale() {},
    setTransform() {}, setLineDash() {}, clip() {}, rect() {}, closePath() {},
};

const dom = new JSDOM(html, {
    runScripts: "dangerously",
    pretendToBeVisual: true,
    beforeParse(window) {
        window.HTMLCanvasElement.prototype.getContext = () => ctxRec;
        window.alert = () => {};
    },
});

const doc = dom.window.document;
let failed = 0;
function check(name, cond, detail) {
    if (cond) { console.log(`  [OK] ${name}`); }
    else { failed++; console.log(`  [FAIL] ${name}  ${detail || ""}`); }
}

setTimeout(() => {
    const W = 940, H = 640, PL = 66, PR = 84, PT = 46, PB = 48;
    console.log("== 绘制几何 ==");
    /* 曲线路径：第一个 stroke 的是格纸竖线（很多），找含 >300 点的路径为曲线 */
    const curvePath = calls.path.filter(p => p.length > 300)[0];
    check("存在理论曲线路径（>300 点）", !!curvePath, `路径数 ${calls.path.length}`);
    if (curvePath) {
        const xs = curvePath.map(p => p[0]), ys = curvePath.map(p => p[1]);
        // 曲线设计上从格纸边缘外进入/穿出（画布自然裁剪），应横跨整个绘图区
        check("曲线横跨绘图区", Math.min(...xs) <= PL + 1 && Math.max(...xs) >= W - PR - 1,
              `x∈[${Math.min(...xs).toFixed(1)}, ${Math.max(...xs).toFixed(1)}]`);
        check("曲线单调下降（左高右低）", ys.every((y, i) => i === 0 || y >= ys[i - 1] - 1e-9));
    }
    /* 经验点 r=4.2；图例标记 r=3.6，用 r>3.9 区分 */
    const arcs = calls.arcs.filter(a => a.r > 3.9);
    check("经验点 30 个", arcs.length === 30, `实际 ${arcs.length}`);
    if (arcs.length === 30) {
        const ok = arcs.every(a => a.x >= PL && a.x <= W - PR && a.y >= PT && a.y <= H - PB);
        check("经验点均在绘图区内", ok);
    }
    const pctLabels = calls.texts.filter(t => /^(0\.01|0\.1|1|5|10|20|30|40|50|60|70|80|90|95|99|99\.9)$/.test(String(t.t).trim()));
    check("频率刻度标签", pctLabels.length >= 12, `实际 ${pctLabels.length}`);
    const TLabels = calls.texts.filter(t => t.y < 35 &&
        /^(10000|1000|200|100|50|20|10|5|2)$/.test(String(t.t).trim()));
    check("重现期刻度标签", TLabels.length === 9, `实际 ${TLabels.length}`);
    const title = calls.texts.find(t => String(t.t).includes("P-Ⅲ型频率曲线"));
    check("标题", !!title);
    const legend = calls.texts.find(t => String(t.t).includes("理论曲线"));
    check("图例参数文字", !!legend, legend && legend.t);

    console.log();
    if (failed) { console.log(`${failed} 项失败`); process.exit(1); }
    console.log("绘图几何测试全部通过 ✔");
    process.exit(0);
}, 400);
