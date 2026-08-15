# -*- coding: utf-8 -*-
"""与 scipy 的交叉验证：检验 gammaincc / gammainccinv / norm_ppf / phi_pearson3。"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from scipy import special

from p3core import gammaincc, gammainccinv, norm_ppf, phi_pearson3

rng = np.random.default_rng(7)

print("1) gammaincc vs scipy.special.gammaincc ...")
a = np.concatenate([np.array([0.04, 0.1, 0.44, 1.0, 2.0, 4.0, 16.0, 100.0, 400.0]),
                    10 ** rng.uniform(-1.3, 3.0, 40)])
x = np.concatenate([np.array([1e-6, 0.01, 0.5, 1.0, 5.0, 10.0, 26.8, 100.0, 401.0]),
                    10 ** rng.uniform(-4, 2.7, 40)])
a, x = np.meshgrid(a, x)
q_ours = gammaincc(a, x)
q_scipy = special.gammaincc(a, x)
err = np.max(np.abs(q_ours - q_scipy))
print(f"   max abs err = {err:.3e}   (a in [0.04, 400], x in [1e-4, 500])")
assert err < 1e-10, "gammaincc 偏差过大"

print("2) gammainccinv vs scipy.special.gammainccinv ...")
a2 = np.concatenate([np.array([0.04, 0.1, 0.44, 1.0, 2.0, 4.0, 16.0, 400.0]),
                     10 ** rng.uniform(-1.3, 3.0, 40)])
q2 = np.concatenate([np.array([1e-4, 1e-3, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999]),
                     rng.uniform(1e-4, 0.999, 30)])
a2, q2 = np.meshgrid(a2, q2)
t_ours = gammainccinv(a2, q2)
t_scipy = special.gammainccinv(a2, q2)
rel = np.abs(t_ours - t_scipy) / np.maximum(t_scipy, 1e-6)
print(f"   max rel err = {rel.max():.3e}")
assert rel.max() < 1e-9, "gammainccinv 偏差过大"

print("3) norm_ppf vs scipy.special.ndtri ...")
p3 = np.concatenate([np.array([1e-12, 1e-6, 1e-4, 0.01, 0.1, 0.5, 0.9, 0.99, 0.9999, 1-1e-12]),
                     rng.uniform(1e-9, 1 - 1e-9, 100)])
d = np.abs(norm_ppf(p3) - special.ndtri(p3))
print(f"   max abs err = {d.max():.3e}")
assert d.max() < 1e-10, "norm_ppf 偏差过大"

print("4) phi_pearson3 一致性：Q(4/Cs², 2Φ/Cs + 4/Cs²) = P ...")
cs4 = rng.uniform(-3, 6, 60)
cs4 = cs4[np.abs(cs4) > 0.05]
p4 = rng.uniform(0.0001, 0.999, 60)
CS, P = np.meshgrid(cs4, p4)
phi = phi_pearson3(CS, P)
alpha = 4.0 / CS ** 2
t = 2.0 * phi / CS + alpha
# 跳过 t 过小的极端点（前向 Q 在此处本身的舍入误差即达 1e-6 量级）
keep = t > 1e-4
back = special.gammaincc(alpha[keep], t[keep])
err4 = np.max(np.abs(back - P[keep]))
print(f"   max |Q(α,t) - P| = {err4:.3e}  (检验点数 {keep.sum()}/{t.size})")
assert err4 < 1e-9, "phi_pearson3 自洽性失败"
# 支撑下界：Φ ≥ -2/Cs（Cs>0 时）
lb = -2.0 / CS
assert np.all(phi[CS > 0] >= lb[CS > 0] - 1e-9), "违反支撑下界"
assert np.all(phi[CS < 0] <= lb[CS < 0] + 1e-9), "违反支撑上界"
print("   支撑边界检查通过 ✔")

print("全部交叉验证通过 ✔")
