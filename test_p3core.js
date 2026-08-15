/* ============================================================================
 * test_p3core.js —— Node 交叉验证：JS 版 p3core 与 Python 版（已验证）比对
 * 运行：node test_p3core.js
 * ==========================================================================*/
"use strict";

const fs = require("fs");
const path = require("path");
const core = require("./p3core.js");

const vec = JSON.parse(fs.readFileSync(path.join(__dirname, "test_vectors.json"), "utf8"));

let failed = 0;
function check(name, cond, detail) {
    if (cond) {
        console.log(`  [OK] ${name}`);
    } else {
        failed++;
        console.log(`  [FAIL] ${name}  ${detail || ""}`);
    }
}

console.log("== 标准正态分位数 norm_ppf ==");
{
    const { p, x } = vec.norm_ppf;
    let maxErr = 0;
    for (let i = 0; i < p.length; i++) {
        const err = Math.abs(core.norm_ppf(p[i]) - x[i]);
        if (err > maxErr) maxErr = err;
    }
    check("vs Python", maxErr < 1e-9, `max err=${maxErr.toExponential(2)}`);
    check("ppf(0.5)=0", Math.abs(core.norm_ppf(0.5)) < 1e-12);
    check("ppf(0.975)≈1.959964", Math.abs(core.norm_ppf(0.975) - 1.959963984540054) < 1e-9);
}

console.log("== 不完全伽马 gammaincc ==");
{
    const { a, x, q } = vec.gammaincc;
    let maxErr = 0;
    for (let i = 0; i < a.length; i++) {
        for (let j = 0; j < x.length; j++) {
            const err = Math.abs(core.gammaincc(a[i], x[j]) - q[i][j]);
            if (err > maxErr) maxErr = err;
        }
    }
    check("vs Python", maxErr < 1e-10, `max err=${maxErr.toExponential(2)}`);
    check("Q(1,x)=e^-x", Math.abs(core.gammaincc(1, 2.0) - Math.exp(-2)) < 1e-12);
}

console.log("== 不完全伽马逆 gammainccinv ==");
{
    const { a, q, t } = vec.gammainccinv;
    let maxRel = 0;
    for (let i = 0; i < a.length; i++) {
        for (let j = 0; j < q.length; j++) {
            const rel = Math.abs(core.gammainccinv(a[i], q[j]) - t[i][j]) / Math.max(t[i][j], 1e-6);
            if (rel > maxRel) maxRel = rel;
        }
    }
    // 容差 1e-8：深尾处两种实现迭代终止条件不同所致（量级 1e-9），非精度问题
    check("vs Python", maxRel < 1e-8, `max rel=${maxRel.toExponential(2)}`);
    // 自洽性
    let maxBack = 0;
    for (const aa of [0.04, 0.44, 1, 4, 16, 400]) {
        for (const qq of [1e-4, 0.01, 0.5, 0.99, 0.999]) {
            const t = core.gammainccinv(aa, qq);
            const back = Math.abs(core.gammaincc(aa, t) - qq);
            if (back > maxBack) maxBack = back;
        }
    }
    check("自洽 Q(a,t)=q", maxBack < 1e-9, `max back=${maxBack.toExponential(2)}`);
}

console.log("== P-Ⅲ 离均系数 phiPearson3 ==");
{
    const { cs, p, phi } = vec.phi;
    let maxErr = 0;
    for (let i = 0; i < cs.length; i++) {
        for (let j = 0; j < p.length; j++) {
            const err = Math.abs(core.phiPearson3(cs[i], p[j]) - phi[i][j]);
            if (err > maxErr) maxErr = err;
        }
    }
    // 容差 1e-8：深尾实现噪声（量级 1e-9）；Φ 表精度仅 2 位小数
    check("vs Python", maxErr < 1e-8, `max err=${maxErr.toExponential(2)}`);
    // 手册值
    check("Cs=1,P=1% ≈ 3.02", Math.abs(core.phiPearson3(1.0, 0.01) - 3.02) < 0.01);
    check("Cs=0.5,P=0.1% ≈ 3.83", Math.abs(core.phiPearson3(0.5, 0.001) - 3.83) < 0.03);
    check("Cs=0 → -ppf", Math.abs(core.phiPearson3(0, 0.01) + core.norm_ppf(0.01)) < 1e-9);
}

console.log("== 矩法统计 ==");
{
    const s = core.sampleStats(vec.stats.x);
    check("mean", Math.abs(s.mean - vec.stats.mean) < 1e-9);
    check("cv", Math.abs(s.cv - vec.stats.cv) < 1e-9);
    check("cs", Math.abs(s.cs - vec.stats.cs) < 1e-9);
}

console.log("== 经验频率 ==");
{
    const e = core.empiricalFreq(vec.stats.x);
    let ok = e.xs.length === vec.emp.xs.length;
    for (let i = 0; i < e.xs.length && ok; i++) {
        if (Math.abs(e.xs[i] - vec.emp.xs[i]) > 1e-9) ok = false;
        if (Math.abs(e.p[i] - vec.emp.p[i]) > 1e-9) ok = false;
    }
    check("排序与 P=m/(n+1)", ok);
}

console.log("== 优化适线 ==");
{
    const r = core.fitCurve(vec.fit_free.x, vec.fit_free.p, null);
    check("均值", Math.abs(r.mean - vec.fit_free.mean) / vec.fit_free.mean < 1e-6,
          `ours=${r.mean.toFixed(4)} py=${vec.fit_free.mean.toFixed(4)}`);
    check("Cv", Math.abs(r.cv - vec.fit_free.cv) / vec.fit_free.cv < 1e-6,
          `ours=${r.cv.toFixed(6)} py=${vec.fit_free.cv.toFixed(6)}`);
    check("Cs", Math.abs(r.cs - vec.fit_free.cs) < 1e-6,
          `ours=${r.cs.toFixed(6)} py=${vec.fit_free.cs.toFixed(6)}`);
    check("objNorm", Math.abs(r.objNorm - vec.fit_free.obj_norm) / vec.fit_free.obj_norm < 1e-6);
    check("rmse", Math.abs(r.rmse - vec.fit_free.rmse) / vec.fit_free.rmse < 1e-6);
    check("r2", Math.abs(r.r2 - vec.fit_free.r2) < 1e-9);
    const rf = core.fitCurve(vec.fit_free.x, vec.fit_free.p, 1200.0);
    check("固定均值 Cs", Math.abs(rf.cs - vec.fit_fixed.cs) < 1e-6);
    check("固定均值 均值", Math.abs(rf.mean - vec.fit_fixed.mean) < 1e-9);
}

console.log("== 性能（曲线重绘所需 phi 调用） ==");
{
    const t0 = Date.now();
    let s = 0;
    for (let i = 0; i < 400; i++) {
        const p = 0.0001 + (0.999 - 0.0001) * i / 399;
        s += core.phiPearson3(1.2, p);
    }
    const dt = Date.now() - t0;
    console.log(`   400 点 phi 耗时 ${dt} ms（${(dt / 400 * 1000).toFixed(1)} µs/点）`);
    const t1 = Date.now();
    core.fitCurve(vec.fit_free.x, vec.fit_free.p, null);
    console.log(`   一次优化适线耗时 ${Date.now() - t1} ms`);
    check("曲线计算 < 100ms", dt < 100);
}

console.log();
if (failed) {
    console.log(`共 ${failed} 项失败`);
    process.exit(1);
}
console.log("JS 核心算法全部通过 ✔");
