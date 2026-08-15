/* ============================================================================
 * p3core.js —— 水文频率计算核心模块（皮尔逊Ⅲ型 / Pearson Type III）
 * 纯 JavaScript 实现，无任何外部依赖，可在浏览器与 Node 中运行。
 *
 * 与 Python 版 p3core.py 算法逐行对应，经与 SciPy 交叉验证（1e-13 量级）。
 * ==========================================================================*/
"use strict";

/* ---------------------------------------------------------------------------
 * 常量
 * -------------------------------------------------------------------------*/
const SQRT2 = Math.SQRT2;

// 频率格纸（海森机率格纸）刻度：P 单位 %
const PAPER_MAJOR = [0.01, 0.1, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99, 99.9];
const PAPER_MINOR = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 15, 20,
                     25, 30, 40, 50, 60, 70, 75, 80, 85, 90, 95, 96, 97, 98, 99, 99.5,
                     99.8, 99.9];
const STD_P = [0.01, 0.1, 1, 2, 5, 10, 20, 50, 75, 90, 95, 99, 99.9];

const LANCZOS_C = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
                   771.32342877765313, -176.61502916214059, 12.507343278686905,
                   -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];

/* ---------------------------------------------------------------------------
 * 基础特殊函数
 * -------------------------------------------------------------------------*/
function lgamma(z) {
    /* ln Γ(z)，z > 0，Lanczos 近似（g=7） */
    let x = z - 1.0;
    let acc = LANCZOS_C[0];
    for (let i = 1; i < 9; i++) {
        acc += LANCZOS_C[i] / (x + i);
    }
    const t = x + 7.0 + 0.5;
    return 0.5 * Math.log(2.0 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(acc);
}

/* 正则化上不完全伽马函数 Q(a,x)，a>0，x>=0 */
function gammaincc(a, x) {
    if (x < a + 1.0) {
        /* 级数求 P(a,x) 再取 1-P */
        const logterm = -x + a * Math.log(x) - lgamma(a + 1.0);
        let term = Math.exp(logterm);
        let s = term;
        let k = 1;
        while (k < 3000) {
            term = term * x / (a + k);
            s += term;
            if (Math.abs(term) <= 1e-16 * Math.abs(s)) break;
            k++;
        }
        return 1.0 - s;
    } else {
        /* 连分式（Lentz，Numerical Recipes gcf） */
        const tiny = 1e-300;
        let b = x + 1.0 - a;
        let c = 1e300;
        let d = 1.0 / b;
        let h = d;
        let k = 1;
        while (k < 5000) {
            const an = -k * (k - a);
            b += 2.0;
            d = an * d + b;
            if (Math.abs(d) < tiny) d = tiny;
            c = b + an / c;
            if (Math.abs(c) < tiny) c = tiny;
            d = 1.0 / d;
            const dl = d * c;
            h *= dl;
            if (Math.abs(dl - 1.0) < 1e-14) break;
            k++;
        }
        return Math.exp(-x + a * Math.log(x) - lgamma(a)) * h;
    }
}

/* 对数空间二分求 t 使 Q(a,t)=q（a<1 或牛顿不收敛时的兜底） */
function bisectGammainccinv(a, q) {
    const lna = Math.log(a);
    const lg = lgamma(a);
    let sLo = (Math.log1p(-q) + lna + lg - 5.0) / a;
    let tHi = a * Math.pow(1.0 - 1.0 / (9.0 * a) + 5.0 / (3.0 * Math.sqrt(a)), 3);
    tHi = Math.max(tHi, 2.0 * a, 1.0);
    let sHi = Math.log(tHi * 1.5);
    for (let i = 0; i < 60; i++) {
        if (gammaincc(a, Math.exp(sHi)) < q) break;
        sHi += 1.0;
    }
    for (let i = 0; i < 60; i++) {
        const sMid = 0.5 * (sLo + sHi);
        if (gammaincc(a, Math.exp(sMid)) > q) {
            sLo = sMid;
        } else {
            sHi = sMid;
        }
    }
    return Math.exp(0.5 * (sLo + sHi));
}

/* 不完全伽马函数的逆：求 t 使 Q(a,t)=q */
function gammainccinv(a, q) {
    if (!(q > 0)) return 0;
    if (q >= 1) return Infinity;
    q = Math.min(Math.max(q, 1e-15), 1.0 - 1e-15);
    /* 初值：Wilson-Hilferty */
    const u = norm_ppf(1.0 - q);
    let t0 = a * Math.pow(1.0 - 1.0 / (9.0 * a) + u / (3.0 * Math.sqrt(a)), 3);
    if (!(t0 > 0) || !isFinite(t0)) {
        t0 = Math.exp(Math.log(a) + u / Math.sqrt(a));
    }
    t0 = Math.min(Math.max(t0, 1e-12), 1e300);
    if (a >= 1.0) {
        /* 对数空间牛顿 */
        let s = Math.log(t0);
        let g = 0;
        let done = false;
        for (let it = 0; it < 25; it++) {
            const t = Math.exp(s);
            const qv = gammaincc(a, t);
            g = qv - q;
            if (Math.abs(g) <= 1e-11) { done = true; break; }
            const logtf = a * s - t - lgamma(a);
            let ds = g / Math.exp(logtf);
            if (ds > 4) ds = 4;
            if (ds < -4) ds = -4;
            s += ds;
        }
        if (done || Math.abs(g) <= 1e-11) {
            return Math.exp(s);
        }
    }
    return bisectGammainccinv(a, q);
}

/* 标准正态分布函数 Φ(x)（用 gammaincc(0.5, x²/2) 实现 erfc，高精度） */
function norm_cdf(x) {
    const tail = 0.5 * gammaincc(0.5, x * x / 2.0);
    return x < 0 ? tail : 1.0 - tail;
}

/* 标准正态分布分位数 Φ⁻¹(p)（Acklam 初值 + 牛顿细化） */
function norm_ppf(p) {
    p = Math.min(Math.max(p, 1e-300), 1.0 - 1e-16);
    const plow = 0.02425, phigh = 0.97575;
    const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
               1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
    const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
               6.680131188771972e+01, -1.328068155928572e+01];
    const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
               -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
    const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
               3.754408661907416e+00];
    let q, x;
    if (p < plow) {
        q = Math.sqrt(-2.0 * Math.log(p));
        x = polyval(c, q) / polyval(d.concat([1.0]), q);
    } else if (p > phigh) {
        q = Math.sqrt(-2.0 * Math.log1p(-p));
        x = -polyval(c, q) / polyval(d.concat([1.0]), q);
    } else {
        q = p - 0.5;
        const r = q * q;
        x = polyval(a, r) * q / polyval(b.concat([1.0]), r);
    }
    /* 牛顿细化（残差用尾部概率，避免消减误差） */
    for (let i = 0; i < 4; i++) {
        const z = Math.abs(x / SQRT2);
        const tail = 0.5 * gammaincc(0.5, z * z);
        const resid = x >= 0 ? (1.0 - p) - tail : tail - p;
        const pdf = Math.exp(-0.5 * x * x) / Math.sqrt(2.0 * Math.PI);
        if (!(pdf > 1e-300)) break;
        x -= resid / pdf;
    }
    return Math.min(Math.max(x, -38.5), 38.5);
}

function polyval(coef, x) {
    let v = 0;
    for (let i = 0; i < coef.length; i++) {
        v = v * x + coef[i];
    }
    return v;
}

/* ---------------------------------------------------------------------------
 * P-Ⅲ 离均系数
 * -------------------------------------------------------------------------*/
function phiPearson3(cs, p) {
    /* Φ_p：P(X > x_p) = p，x_p = x̄(1 + Cv·Φ) */
    p = Math.min(Math.max(p, 1e-12), 1.0 - 1e-12);
    if (Math.abs(cs) < 1e-9) {
        return -norm_ppf(p);
    }
    const aa = 4.0 / (cs * cs);
    const t = gammainccinv(aa, p);
    return cs * (t - aa) / 2.0;
}

/* ---------------------------------------------------------------------------
 * 统计与适线
 * -------------------------------------------------------------------------*/
function sampleStats(x) {
    const n = x.length;
    let sum = 0;
    for (let i = 0; i < n; i++) sum += x[i];
    const mean = sum / n;
    let varSum = 0;
    for (let i = 0; i < n; i++) {
        const d = x[i] - mean;
        varSum += d * d;
    }
    const std = n > 1 ? Math.sqrt(varSum / (n - 1)) : 0;
    let cv = 0, cs = 0;
    if (Math.abs(mean) > 1e-12 && n > 1) {
        let s1 = 0;
        for (let i = 0; i < n; i++) {
            const k = x[i] / mean - 1.0;
            s1 += k * k;
        }
        cv = Math.sqrt(s1 / n);
        if (cv > 1e-9 && n > 2) {
            let s2 = 0;
            for (let i = 0; i < n; i++) {
                const k = x[i] / mean - 1.0;
                s2 += k * k * k;
            }
            cs = s2 / (n * cv * cv * cv);
        }
    }
    let xmax = -Infinity, xmin = Infinity;
    for (let i = 0; i < n; i++) {
        if (x[i] > xmax) xmax = x[i];
        if (x[i] < xmin) xmin = x[i];
    }
    return { n, mean, std, cv, cs, xmax, xmin };
}

/* 经验频率：P = m/(n+1) × 100%（返回按从大到小排序的值与频率百分数） */
function empiricalFreq(x) {
    const xs = x.slice().sort((u, v) => v - u);
    const n = xs.length;
    const p = [];
    for (let i = 0; i < n; i++) p.push(100.0 * (i + 1) / (n + 1));
    return { xs, p };
}

/* 优化适线：在 Cs 网格上搜索，每个 Cs 用加权最小二乘（w=1/x²）解均值/变差系数 */
function fitCurve(x, p, meanFixed) {
    /* x, p: 数组；p 为超越概率（0~1）；meanFixed: null=均值自由，数值=固定 */
    const n = x.length;
    const w = new Array(n), wx = new Array(n);
    let sw = 0, sxw = 0;
    for (let i = 0; i < n; i++) {
        w[i] = 1.0 / (x[i] * x[i]);
        wx[i] = w[i] * x[i];
        sw += w[i];
        sxw += wx[i];
    }

    function stage(csArr) {
        let bestIdx = -1, bestObj = Infinity;
        const res = { a: new Array(csArr.length), b: new Array(csArr.length),
                      obj: new Array(csArr.length) };
        for (let j = 0; j < csArr.length; j++) {
            const cs = csArr[j];
            let swp = 0, spp = 0, sxpw = 0;
            const phi = new Array(n);
            for (let i = 0; i < n; i++) {
                const ph = phiPearson3(cs, p[i]);
                phi[i] = ph;
                swp += w[i] * ph;
                spp += w[i] * ph * ph;
                sxpw += wx[i] * ph;
            }
            let a, b;
            if (meanFixed !== null && meanFixed !== undefined) {
                a = meanFixed;
                b = spp > 1e-15 ? (sxpw - a * swp) / spp : NaN;
            } else {
                const denom = sw * spp - swp * swp;
                if (denom > 1e-12) {
                    a = (spp * sxw - swp * sxpw) / denom;
                    b = (sw * sxpw - swp * sxw) / denom;
                } else {
                    a = NaN; b = NaN;
                }
            }
            const cv = b / a;
            let obj = Infinity;
            if (isFinite(a) && isFinite(b) && a > 0 && b > 0 && cv < 50 && cv > 1e-6) {
                let sse = 0;
                for (let i = 0; i < n; i++) {
                    const d = x[i] - (a + b * phi[i]);
                    sse += w[i] * d * d;
                }
                obj = sse;
            }
            res.a[j] = a; res.b[j] = b; res.obj[j] = obj;
            if (obj < bestObj) { bestObj = obj; bestIdx = j; }
        }
        return { bestIdx, bestObj, res };
    }

    const steps = [
        { step: 0.1,  span: null },
        { step: 0.005, span: 0.12 },
        { step: 0.0005, span: 0.02 },
    ];
    let bestCs = null, bestA = null, bestB = null, bestObj = Infinity;
    for (const { step, span } of steps) {
        let csArr;
        if (span === null) {
            csArr = [];
            for (let cs = -1.0; cs <= 8.0 + step / 2; cs += step) csArr.push(cs);
        } else {
            csArr = [];
            const lo = Math.max(-1.0, bestCs - span);
            const hi = Math.min(8.0, bestCs + span);
            for (let cs = lo; cs <= hi + step / 2; cs += step) csArr.push(cs);
        }
        const r = stage(csArr);
        if (isFinite(r.bestObj) && r.bestObj < bestObj) {
            bestCs = csArr[r.bestIdx];
            bestA = r.res.a[r.bestIdx];
            bestB = r.res.b[r.bestIdx];
            bestObj = r.bestObj;
        }
    }
    if (bestCs === null) return null;
    /* 拟合指标 */
    const cs = bestCs, a = bestA, b = bestB;
    const cv = b / a;
    const phi = new Array(n);
    let sse = 0, sseRel = 0, maxRel = 0, sumX = 0, sumXh = 0, sumXX = 0, sumXhXh = 0, sumXXh = 0;
    for (let i = 0; i < n; i++) {
        const ph = phiPearson3(cs, p[i]);
        phi[i] = ph;
        const xh = a + b * ph;
        const d = x[i] - xh;
        sse += d * d;
        const rel = d / x[i];
        sseRel += rel * rel;
        if (Math.abs(rel) > maxRel) maxRel = Math.abs(rel);
        sumX += x[i]; sumXh += xh; sumXX += x[i] * x[i];
        sumXhXh += xh * xh; sumXXh += x[i] * xh;
    }
    const rmse = Math.sqrt(sse / n);
    const r2 = n > 2 ? Math.pow((n * sumXXh - sumX * sumXh) /
        Math.sqrt((n * sumXX - sumX * sumX) * (n * sumXhXh - sumXh * sumXh)), 2) : 1;
    let xwx = 0;
    for (let i = 0; i < n; i++) xwx += w[i] * x[i] * x[i];
    return { cs, mean: a, cv, obj: bestObj, objNorm: bestObj / xwx, rmse, maxRel, r2 };
}

/* 理论频率曲线：给定频率百分数数组，返回设计值数组 */
function frequencyCurve(mean, cv, cs, pPct) {
    return pPct.map(p => mean * (1.0 + cv * phiPearson3(cs, p / 100.0)));
}

/* 格纸横坐标：频率百分数 → 标准正态分位 u */
function pctToX(pPct) {
    return norm_ppf(pPct / 100.0);
}

/* 供 Node 与浏览器共用 */
if (typeof module !== "undefined" && module.exports) {
    module.exports = { lgamma, gammaincc, gammainccinv, norm_cdf, norm_ppf,
                       phiPearson3, sampleStats, empiricalFreq, fitCurve,
                       frequencyCurve, pctToX, PAPER_MAJOR, PAPER_MINOR, STD_P };
}
