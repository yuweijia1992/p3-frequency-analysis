# -*- coding: utf-8 -*-
"""
P3Core —— 水文频率计算核心模块（皮尔逊Ⅲ型 / Pearson Type III）
================================================================
提供：
  * 矩法统计量             sample_stats
  * 经验频率               empirical_frequency   （P = m/(n+1) × 100%，数学期望公式）
  * 离均系数 Φ_p           phi_pearson3          （皮尔逊Ⅲ型频率曲线离均系数）
  * 理论频率曲线           frequency_curve / curve_sample_u
  * 优化适线（最小二乘）   fit_curve
  * 频率格纸（海森格纸）   pct_to_x / PAPER_MAJOR / PAPER_MINOR
  * 标准设计频率           STD_P

原理：
    皮尔逊Ⅲ型密度函数
        f(x) = β^α / Γ(α) · (x - a₀)^(α-1) · e^(-β(x-a₀))
    以均值 x̄、变差系数 Cv、偏态系数 Cs 表示：
        α = 4/Cs²,  β = 2/(x̄·Cv·Cs),  a₀ = x̄(1 - 2Cv/Cs)
    令 Y = β(X - a₀) ~ Γ(α,1)，标准化变量 Φ = (x - x̄)/(x̄·Cv)，
    则 超越概率  P(X > x_p) = Q(α, t)，其中 t = 2Φ/Cs + 4/Cs²，Q 为上不完全伽马函数。
    给定 P 与 Cs，先求 t = Q⁻¹(α, P)，再得 Φ = Cs(t - α)/2。
    当 Cs = 0 时退化为正态分布：Φ = -Φ⁻¹(P)（标准正态分位数的相反数）。

    本模块自行实现上不完全伽马函数 Q(a,x) 及其逆（连分式 + 级数 + 牛顿/二分法），
    不依赖 scipy，便于打包为轻量级 Windows 可执行程序。
"""

import math

import numpy as np

# ---------------------------------------------------------------------------
# 常数
# ---------------------------------------------------------------------------
_SQRT2 = math.sqrt(2.0)

# 频率格纸（海森机率格纸）刻度：P 单位 %
PAPER_MAJOR = [0.01, 0.1, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99, 99.9]
PAPER_MINOR = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 15, 20,
               25, 30, 40, 50, 60, 70, 75, 80, 85, 90, 95, 96, 97, 98, 99, 99.5,
               99.8, 99.9]

# 成果表中给出的标准设计频率（%）
STD_P = [0.01, 0.1, 1, 2, 5, 10, 20, 50, 75, 90, 95, 99, 99.9]

_LANCZOS_C = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
              771.32342877765313, -176.61502916214059, 12.507343278686905,
              -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]


# ---------------------------------------------------------------------------
# 基础特殊函数（全部向量化，仅依赖 numpy / math）
# ---------------------------------------------------------------------------
def _lgamma(z):
    """ln Γ(z)，z > 0，Lanczos 近似（g=7），精度 ~1e-14。"""
    z = np.asarray(z, float)
    x = z - 1.0
    acc = np.full_like(x, _LANCZOS_C[0])
    for i in range(1, 9):
        acc = acc + _LANCZOS_C[i] / (x + i)
    t = x + 7.0 + 0.5
    return 0.5 * math.log(2.0 * math.pi) + (x + 0.5) * np.log(t) - t + np.log(acc)


_erf_vec = np.frompyfunc(math.erf, 1, 1)
_erfc_vec = np.frompyfunc(math.erfc, 1, 1)


def norm_cdf(x):
    """标准正态分布函数 Φ(x)（向量化，尾部用 erfc 避免消减误差）。"""
    x = np.asarray(x, float)
    scalar = x.ndim == 0
    x = np.atleast_1d(x)
    z = x / _SQRT2
    y = np.where(z < 0.0, 0.5 * _erfc_vec(-z).astype(float),
                 1.0 - 0.5 * _erfc_vec(z).astype(float))
    return y[0] if scalar else y


def norm_ppf(p):
    """标准正态分布分位数 Φ⁻¹(p)（Acklam 有理逼近 + 牛顿细化，双精度，向量化）。"""
    p = np.asarray(p, float)
    scalar = p.ndim == 0
    p = np.atleast_1d(p)
    p = np.clip(p, 1e-300, 1.0 - 1e-16)
    plow, phigh = 0.02425, 0.97575
    a = np.array([-3.969683028665376e+01, 2.209460984245205e+02,
                  -2.759285104469687e+02, 1.383577518672690e+02,
                  -3.066479806614716e+01, 2.506628277459239e+00])
    b = np.array([-5.447609879822406e+01, 1.615858368580409e+02,
                  -1.556989798598866e+02, 6.680131188771972e+01,
                  -1.328068155928572e+01])
    c = np.array([-7.784894002430293e-03, -3.223964580411365e-01,
                  -2.400758277161838e+00, -2.549732539343734e+00,
                  4.374664141464968e+00, 2.938163982698783e+00])
    d = np.array([7.784695709041462e-03, 3.224671290700398e-01,
                  2.445134137142996e+00, 3.754408661907416e+00])

    q = np.where(p < plow, np.sqrt(-2.0 * np.log(p)),
                 np.where(p > phigh, np.sqrt(-2.0 * np.log1p(-p)), p - 0.5))
    r = q * q
    xc = np.polyval(a, r) * q / np.polyval(np.append(b, 1.0), r)
    xt = np.polyval(c, q) / np.polyval(np.append(d, 1.0), q)
    x = np.where(p < plow, xt, np.where(p > phigh, -xt, xc))
    # 牛顿细化至双精度（残差用尾部概率直接计算，避免 1-cdf 消减误差）
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        for _ in range(4):
            z = np.abs(x / _SQRT2)
            tail = 0.5 * _erfc_vec(z).astype(float)   # 尾部概率，高精度
            resid = np.where(x >= 0.0, (1.0 - p) - tail, tail - p)
            pdf = np.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)
            ok = np.isfinite(pdf) & (pdf > 1e-300)
            x = np.where(ok, x - resid / np.where(ok, pdf, 1.0), x)
        x = np.clip(x, -38.5, 38.5)
    return x[0] if scalar else x


def gammaincc(a, x):
    """正则化上不完全伽马函数 Q(a,x) = ∫ₓ^∞ t^(a-1) e^(-t) dt / Γ(a)。

    参数：a, x 为同形状数组（可广播），a > 0，x >= 0。
    对 x < a+1 用级数求下不完全伽马 P 再取 1-P；
    对 x >= a+1 用连分式（Lentz 算法）直接求 Q。
    """
    a = np.asarray(a, float)
    x = np.asarray(x, float)
    a, x = np.broadcast_arrays(a, x)
    shape = a.shape
    a = a.ravel()
    x = x.ravel()
    out = np.empty(a.shape)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        m_series = x < a + 1.0
        idx = np.nonzero(m_series)[0]
        if idx.size:
            aa, xx = a[idx], x[idx]
            # 级数：P(a,x) = Σ term_k, term_0 = e^{-x}·x^a/Γ(a+1),
            #       term_k = term_{k-1}·x/(a+k)
            logterm = -xx + aa * np.log(xx) - _lgamma(aa + 1.0)
            term = np.exp(logterm)
            s = term.copy()
            k = 1
            conv = np.zeros_like(term, dtype=bool)
            while not conv.all() and k < 3000:
                term = term * xx / (aa + k)
                s = s + term
                conv = np.abs(term) <= 1e-16 * np.abs(s)
                k += 1
            out[idx] = 1.0 - s
        idx = np.nonzero(~m_series)[0]
        if idx.size:
            aa, xx = a[idx], x[idx]
            # 连分式（Numerical Recipes: gcf）
            tiny = 1e-300
            b = xx + 1.0 - aa
            c = np.full_like(b, 1e300)
            d = 1.0 / b
            h = d.copy()
            k = 1
            conv = np.zeros_like(b, dtype=bool)
            while not conv.all() and k < 5000:
                an = -k * (k - aa)
                b = b + 2.0
                d = an * d + b
                d = np.where(np.abs(d) < tiny, tiny, d)
                c = b + an / c
                c = np.where(np.abs(c) < tiny, tiny, c)
                d = 1.0 / d
                dl = d * c
                h = h * dl
                conv = np.abs(dl - 1.0) < 1e-14
                k += 1
            logpre = -xx + aa * np.log(xx) - _lgamma(aa)
            out[idx] = np.exp(logpre) * h
    return out.reshape(shape)


def _bisect_gammainccinv(a, q):
    """在 log t 空间对 Q(a,t) = q 二分求解（向量化，稳健兜底）。"""
    lna = np.log(a)
    lg = _lgamma(a)
    # 下界：小 t 渐近 P(a,t) ≈ t^a/(a·Γ(a))，取 P = (1-q)/e^5
    s_lo = (np.log1p(-q) + lna + lg - 5.0) / a
    # 上界：Wilson-Hilferty 近似，u=5 留足裕度
    t_hi = a * (1.0 - 1.0 / (9.0 * a) + 5.0 / (3.0 * np.sqrt(a))) ** 3
    t_hi = np.maximum(t_hi, 2.0 * a)
    t_hi = np.maximum(t_hi, 1.0)
    s_hi = np.log(t_hi * 1.5)
    with np.errstate(over="ignore", invalid="ignore"):
        # 确保 Q(s_hi) < q，不足则外扩
        for _ in range(60):
            expand = gammaincc(a, np.exp(s_hi)) >= q
            if not expand.any():
                break
            s_hi = np.where(expand, s_hi + 1.0, s_hi)
        for _ in range(60):
            s_mid = 0.5 * (s_lo + s_hi)
            qm = gammaincc(a, np.exp(s_mid))
            move_lo = qm > q
            s_lo = np.where(move_lo, s_mid, s_lo)
            s_hi = np.where(move_lo, s_hi, s_mid)
    return np.exp(0.5 * (s_lo + s_hi))


def gammainccinv(a, q, tol=1e-11):
    """正则化上不完全伽马函数的逆：求 t 使 Q(a,t) = q（向量化）。

    a >= 1 时用对数空间阻尼牛顿法（Wilson-Hilferty 初值），未收敛则二分；
    a <  1 时直接用对数空间二分。
    """
    a = np.asarray(a, float)
    q = np.asarray(q, float)
    a, q = np.broadcast_arrays(a, q)
    shape = a.shape
    a = a.ravel()
    q = q.ravel()
    q = np.clip(q, 1e-15, 1.0 - 1e-15)
    t = np.empty_like(a)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        # ---- 初值：Wilson-Hilferty ----
        u = norm_ppf(1.0 - q)
        t0 = a * (1.0 - 1.0 / (9.0 * a) + u / (3.0 * np.sqrt(a))) ** 3
        bad = ~(np.isfinite(t0) & (t0 > 0.0))
        t0 = np.where(bad, np.exp(np.log(a) + u / np.sqrt(a)), t0)
        t0 = np.clip(t0, 1e-12, 1e300)
        # ---- a >= 1：对数空间牛顿 ----
        m_n = a >= 1.0
        if m_n.any():
            am, qm = a[m_n], q[m_n]
            s = np.log(t0[m_n])
            done = np.zeros_like(s, dtype=bool)
            g = np.zeros_like(s)
            for _ in range(25):
                if done.all():
                    break
                tv = np.exp(s)
                qv = gammaincc(am, tv)
                g = qv - qm
                logtf = am * s - tv - _lgamma(am)
                ds = np.clip(g / np.exp(logtf), -4.0, 4.0)
                s = s + ds
                done = np.abs(g) <= tol
            t[m_n] = np.exp(s)
            # 未收敛的用二分兜底
            rem = ~done
            if rem.any():
                idx_all = np.nonzero(m_n)[0]
                idx_rem = idx_all[rem]
                t[idx_rem] = _bisect_gammainccinv(a[idx_rem], q[idx_rem])
        # ---- a < 1：二分 ----
        m_b = ~m_n
        if m_b.any():
            t[m_b] = _bisect_gammainccinv(a[m_b], q[m_b])
    return t.reshape(shape)


# ---------------------------------------------------------------------------
# P-Ⅲ 离均系数与频率曲线
# ---------------------------------------------------------------------------
def phi_pearson3(cs, p):
    """皮尔逊Ⅲ型离均系数 Φ_p。

    参数
    ----
    cs : 偏态系数（标量或数组，|cs| 不能为 0，极小值时按正态近似）
    p  : 超越概率（0 < p < 1，标量或数组，可与 cs 广播）

    返回 Φ，满足 x_p = x̄·(1 + Cv·Φ)，且 P(X > x_p) = p。
    """
    cs = np.asarray(cs, float)
    p = np.asarray(p, float)
    cs, p = np.broadcast_arrays(cs, p)
    shape = cs.shape
    cs = cs.ravel()
    p = p.ravel()
    p = np.clip(p, 1e-12, 1.0 - 1e-12)
    out = np.empty_like(cs)
    m0 = np.abs(cs) < 1e-9
    if m0.any():
        out[m0] = -norm_ppf(p[m0])
    m1 = ~m0
    if m1.any():
        c = cs[m1]
        aa = 4.0 / (c * c)
        t = gammainccinv(aa, p[m1])
        out[m1] = c * (t - aa) / 2.0
    return out.reshape(shape)


def phi_grid(cs_1d, p_1d):
    """二维离均系数表：cs_1d (N,) × p_1d (M,) → (N, M)。"""
    cs = np.asarray(cs_1d, float)
    p = np.asarray(p_1d, float)
    return phi_pearson3(cs[:, None], p[None, :])


def frequency_curve(mean, cv, cs, p_pct):
    """给定均值、Cv、Cs，求指定频率（%）对应的理论值数组。"""
    p = np.asarray(p_pct, float) / 100.0
    return mean * (1.0 + cv * phi_pearson3(cs, p))


def curve_sample_u(mean, cv, cs, u_lo=-4.3, u_hi=3.6, n=400):
    """在频率格纸横坐标 u（标准正态分位）上等距采样，返回 (u, y) 用于绘图。

    u 与频率 P 的关系：P = Φ(u)（超越概率）；格纸横坐标即 u。
    """
    u = np.linspace(u_lo, u_hi, n)
    p = norm_cdf(u)
    return u, mean * (1.0 + cv * phi_pearson3(cs, p))


# ---------------------------------------------------------------------------
# 统计与适线
# ---------------------------------------------------------------------------
def sample_stats(x):
    """矩法统计量。Cv、Cs 采用水文矩法公式（分母用 n）：
        x̄ = Σx/n
        Cv = sqrt(Σ(Ki-1)²/n)，Ki = xi/x̄
        Cs = Σ(Ki-1)³/(n·Cv³)
    同时给出样本标准差 σ（n-1 无偏）供参考。
    """
    x = np.asarray(x, float)
    n = x.size
    mean = float(x.mean())
    std = float(x.std(ddof=1)) if n > 1 else 0.0
    if abs(mean) > 1e-12 and n > 1:
        k = x / mean
        cv = float(np.sqrt(((k - 1.0) ** 2).sum() / n))
    else:
        cv = 0.0
    if cv > 1e-9 and n > 2:
        k = x / mean
        cs = float(((k - 1.0) ** 3).sum() / (n * cv ** 3))
    else:
        cs = 0.0
    return dict(n=n, mean=mean, std=std, cv=cv, cs=cs,
                xmax=float(x.max()), xmin=float(x.min()))


def empirical_frequency(x):
    """经验频率（数学期望公式，Weibull）：P_m = m/(n+1) × 100%。

    返回 (x_sorted_desc, p_pct)，x 按从大到小排列，p 为超越概率百分数。
    """
    xs = np.sort(np.asarray(x, float))[::-1]
    n = xs.size
    p = 100.0 * np.arange(1, n + 1) / (n + 1.0)
    return xs, p


def fit_curve(x, p, mean_fixed=None, cs_lo=-1.0, cs_hi=8.0,
              cs_step=0.1, refine_step=0.005, refine2_step=0.0005):
    """优化适线：在 Cs 网格上搜索，对每个 Cs 用加权最小二乘（权重 w=1/x²）
    解出最优均值/变差系数（x̂ = a + b·Φ，a = x̄，b = x̄·Cv），
    目标函数为相对残差平方和 Σ((x - x̂)/x)²，取全局最优。

    参数
    ----
    x, p : 经验点（p 为超越概率，0~1）
    mean_fixed : 若给定，则均值固定为该值（经典适线法做法，默认 None 表示三参数均可变）
    """
    x = np.asarray(x, float)
    p = np.asarray(p, float)
    w = 1.0 / (x * x)
    sw = float(w.sum())
    sxw = float((w * x).sum())
    xwx = float((w * x * x).sum())

    def stage(cs_arr):
        cs_arr = np.asarray(cs_arr, float)
        phi = phi_pearson3(cs_arr[:, None], p[None, :])          # (N, M)
        swp = phi @ w
        spp = (phi * phi) @ w
        sxpw = phi @ (w * x)
        if mean_fixed is not None:
            a = np.full(cs_arr.shape, float(mean_fixed))
            b = (sxpw - a * swp) / np.where(spp > 1e-15, spp, 1.0)
        else:
            denom = sw * spp - swp * swp
            ok = denom > 1e-12
            den = np.where(ok, denom, 1.0)
            a = (spp * sxw - swp * sxpw) / den
            b = (sw * sxpw - swp * sxw) / den
            a = np.where(ok, a, np.nan)
            b = np.where(ok, b, np.nan)
        cv = b / a
        valid = (np.isfinite(a) & np.isfinite(b) & (a > 0.0) & (b > 0.0)
                 & (cv < 50.0) & (cv > 1e-6))
        xhat = a[:, None] + b[:, None] * phi
        resid = x[None, :] - xhat
        obj = np.where(valid, (w[None, :] * resid * resid).sum(axis=1), np.inf)
        return a, b, obj

    best = None
    for step, span in ((cs_step, None), (refine_step, 0.12), (refine2_step, 0.02)):
        if span is None:
            cs_arr = np.arange(cs_lo, cs_hi + step * 0.5, step)
        else:
            lo = max(cs_lo, best_cs - span)
            hi = min(cs_hi, best_cs + span)
            cs_arr = np.arange(lo, hi + step * 0.5, step)
        a, b, obj = stage(cs_arr)
        i = int(np.argmin(obj))
        if not np.isfinite(obj[i]):
            continue
        best_cs = float(cs_arr[i])
        best_a, best_b, best_obj = float(a[i]), float(b[i]), float(obj[i])
        best = (best_cs, best_a, best_b, best_obj)

    if best is None:
        return None
    cs, a, b, obj = best
    cv = b / a
    phi = phi_pearson3(cs, p)
    xhat = a + b * phi
    resid = x - xhat
    rmse = float(np.sqrt((resid ** 2).mean()))
    max_rel = float(np.max(np.abs(resid / x)))
    r2 = float(np.corrcoef(x, xhat)[0, 1] ** 2) if x.size > 2 else 1.0
    return dict(cs=cs, mean=a, cv=cv, obj=obj, obj_norm=obj / xwx,
                phi=phi, xhat=xhat, resid=resid,
                rmse=rmse, max_rel=max_rel, r2=r2)


def pct_to_x(p_pct):
    """频率格纸横坐标：频率百分数 → 标准正态分位 u。"""
    return norm_ppf(np.asarray(p_pct, float) / 100.0)
