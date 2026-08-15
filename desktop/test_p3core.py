# -*- coding: utf-8 -*-
"""p3core 单元测试：验证离均系数、伽马函数及其逆、矩法、适线的数值精度。"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from p3core import (norm_ppf, norm_cdf, gammaincc, gammainccinv, phi_pearson3,
                    phi_grid, frequency_curve, curve_sample_u, sample_stats,
                    empirical_frequency, fit_curve, pct_to_x, STD_P)

FAIL = []


def check(name, cond, detail=""):
    if cond:
        print(f"  [OK] {name}")
    else:
        FAIL.append(name)
        print(f"  [FAIL] {name}  {detail}")


def test_norm():
    print("== 标准正态分位数 norm_ppf ==")
    check("ppf(0.5)=0", abs(norm_ppf(0.5)) < 1e-12)
    check("ppf(0.975)≈1.959964", abs(norm_ppf(0.975) - 1.959963984540054) < 1e-7)
    check("ppf(0.01)≈-2.3263", abs(norm_ppf(0.01) + 2.3263478740408408) < 1e-7)
    check("ppf(1e-4)≈-3.719", abs(norm_ppf(1e-4) + 3.7190164854557088) < 1e-6)
    check("ppf(0.9999)≈3.719", abs(norm_ppf(0.9999) - 3.7190164854557088) < 1e-6)
    # 逆一致性
    u = np.array([-4.0, -2.5, -1.0, 0.0, 0.7, 2.2, 3.8])
    check("ppf-cdf 互逆", np.max(np.abs(norm_ppf(norm_cdf(u)) - u)) < 1e-9)


def test_gamma():
    print("== 不完全伽马函数 ==")
    # Q(4,10) = e^-10 · Σ_{k=0..3} 10^k/k!
    q410 = np.exp(-10.0) * (1 + 10 + 50 + 1000 / 6.0)
    check("Q(4,10)≈0.01034", abs(gammaincc(4.0, 10.0) - q410) < 1e-12)
    check("Q(1,x)=e^-x", np.max(np.abs(gammaincc(1.0, np.array([0.5, 2.0, 5.0]))
                                - np.exp(-np.array([0.5, 2.0, 5.0])))) < 1e-12)
    # 端点
    check("Q(a,0)=1", abs(gammaincc(2.5, 0.0) - 1.0) < 1e-12)
    # 大参数
    check("Q(400, 401)≈0.47344", abs(gammaincc(400.0, 401.0) - 0.473441079) < 1e-6)
    check("Q(16, 26.8)≈0.01", abs(gammaincc(16.0, 26.8) - 0.01) < 2e-3)

    print("== 不完全伽马逆 gammainccinv ==")
    for a in [0.04, 0.44, 1.0, 4.0, 16.0, 400.0]:
        for q in [1e-4, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999]:
            t = gammainccinv(a, q)
            err = abs(gammaincc(a, t) - q)
            check(f"roundtrip a={a} q={q}", err < 1e-9, f"err={err:.2e}")


def test_phi():
    print("== P-Ⅲ 离均系数 phi_pearson3 ==")
    # Cs=0 → 正态
    p = np.array([0.0001, 0.01, 0.05, 0.5, 0.95, 0.99])
    check("Cs=0 → -norm_ppf", np.max(np.abs(phi_pearson3(0.0, p) + norm_ppf(p))) < 1e-6)
    # 查表值：Cs=1.0, P=1% → Φ≈3.02（水文手册离均系数表）
    v = phi_pearson3(1.0, 0.01)
    check("Cs=1, P=1% ≈ 3.02", abs(v - 3.02) < 0.01, f"got {v:.4f}")
    # Cs=0.5, P=0.1% → Φ≈3.83
    v = phi_pearson3(0.5, 0.001)
    check("Cs=0.5, P=0.1% ≈ 3.83", abs(v - 3.83) < 0.03, f"got {v:.4f}")
    # Cs=1, P=0.1% → Φ≈4.53
    v = phi_pearson3(1.0, 0.001)
    check("Cs=1, P=0.1% ≈ 4.53", abs(v - 4.53) < 0.05, f"got {v:.4f}")
    # Cs=2, P=1% → Φ≈3.61
    v = phi_pearson3(2.0, 0.01)
    check("Cs=2, P=1% ≈ 3.61", abs(v - 3.61) < 0.03, f"got {v:.4f}")
    # Cs=0.5, P=50% → Φ≈-0.08
    v = phi_pearson3(0.5, 0.5)
    check("Cs=0.5, P=50% ≈ -0.08", abs(v + 0.08) < 0.02, f"got {v:.4f}")
    # 负 Cs
    check("Cs=-1, P=1% ≈ -3.02", abs(phi_pearson3(-1.0, 0.01) + 3.02) < 0.01)
    # 下界：P→1 时 Φ→-2/Cs
    check("下限 Φ→-2/Cs", abs(phi_pearson3(1.5, 1 - 1e-6) + 2 / 1.5) < 0.01)
    # 网格形状
    g = phi_grid(np.array([0.5, 1.0]), np.array([0.01, 0.5, 0.99]))
    check("phi_grid 形状", g.shape == (2, 3))


def test_moments():
    print("== 矩法统计 ==")
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    s = sample_stats(x)
    check("mean=3", abs(s["mean"] - 3.0) < 1e-12)
    check("std≈1.5811", abs(s["std"] - 1.581138830) < 1e-6)
    k = np.array([1.0, 2.0, 3.0, 4.0, 5.0]) / 3.0
    cv = np.sqrt(((k - 1) ** 2).sum() / 5)
    check("Cv(矩法)", abs(s["cv"] - cv) < 1e-12)
    cs = ((k - 1) ** 3).sum() / (5 * cv ** 3)
    check("Cs(矩法)", abs(s["cs"] - cs) < 1e-9)
    xs, p = empirical_frequency(x)
    check("经验频率排序", np.allclose(xs, [5, 4, 3, 2, 1]))
    check("经验频率 P=m/(n+1)", np.allclose(p, [100 / 6, 200 / 6, 300 / 6, 400 / 6, 500 / 6]))


def test_fit():
    print("== 优化适线 ==")
    rng = np.random.default_rng(42)
    mean, cv, cs = 1200.0, 0.45, 1.2
    p = np.arange(1, 41) / 41.0
    phi = phi_pearson3(cs, p)
    x = mean * (1 + cv * phi) * (1 + rng.normal(0, 0.01, 40))
    r = fit_curve(x, p)
    check("适线 Cs 恢复", abs(r["cs"] - cs) < 0.15, f"cs={r['cs']:.3f}")
    check("适线 Cv 恢复", abs(r["cv"] - cv) < 0.05, f"cv={r['cv']:.3f}")
    check("适线 均值恢复", abs(r["mean"] - mean) / mean < 0.03, f"mean={r['mean']:.1f}")
    # 均值固定
    r2 = fit_curve(x, p, mean_fixed=mean)
    check("固定均值后均值不变", abs(r2["mean"] - mean) < 1e-9)
    # 曲线与表
    y = frequency_curve(mean, cv, cs, np.array(STD_P))
    check("成果表曲线计算", y.shape == (len(STD_P),) and np.all(np.isfinite(y)))


def test_curve():
    print("== 频率曲线采样 ==")
    u, y = curve_sample_u(1000.0, 0.5, 1.0)
    check("采样形状", u.shape == y.shape and u.size == 400)
    # 单调性：u 越大（频率越大）设计值越小
    check("曲线单调", np.all(np.diff(y) < 0))
    # 中位数处 P=50%: u=0
    p50 = norm_cdf(0.0)
    v50 = 1000.0 * (1 + 0.5 * phi_pearson3(1.0, p50))
    idx = np.argmin(np.abs(u))
    check("曲线过中位数点", abs(y[idx] - v50) / v50 < 0.02)


if __name__ == "__main__":
    test_norm()
    test_gamma()
    test_phi()
    test_moments()
    test_fit()
    test_curve()
    print()
    if FAIL:
        print(f"共 {len(FAIL)} 项失败: {FAIL}")
        sys.exit(1)
    print("全部测试通过 ✔")
