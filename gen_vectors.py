# -*- coding: utf-8 -*-
"""用已验证的 Python p3core 生成测试向量 test_vectors.json，供 Node 版交叉验证。"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "P3FrequencyTool"))

import numpy as np
from p3core import (norm_ppf, gammaincc, gammainccinv, phi_pearson3,
                    fit_curve, sample_stats, empirical_frequency)

out = {}

# 1) norm_ppf
p_norm = [1e-12, 1e-8, 1e-4, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.9999, 1 - 1e-12]
out["norm_ppf"] = {"p": p_norm,
                   "x": [float(norm_ppf(p)) for p in p_norm]}

# 2) gammaincc
a_g = [0.04, 0.1, 0.44, 1.0, 2.0, 4.0, 16.0, 400.0]
x_g = [1e-4, 0.01, 0.5, 1.0, 5.0, 10.0, 26.8, 100.0, 401.0]
out["gammaincc"] = {"a": a_g, "x": x_g,
                    "q": [[float(gammaincc(a, x)) for x in x_g] for a in a_g]}

# 3) gammainccinv
q_i = [1e-4, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999]
out["gammainccinv"] = {"a": a_g, "q": q_i,
                       "t": [[float(gammainccinv(a, q)) for q in q_i] for a in a_g]}

# 4) phi_pearson3
cs_p = [-1.0, -0.5, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
p_p = [1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.99, 0.999]
out["phi"] = {"cs": cs_p, "p": p_p,
              "phi": [[float(phi_pearson3(c, p)) for p in p_p] for c in cs_p]}

# 5) 适线：合成数据（固定种子）
rng = np.random.default_rng(42)
mean0, cv0, cs0 = 1200.0, 0.45, 1.2
p_fit = np.arange(1, 41) / 41.0
phi_fit = phi_pearson3(cs0, p_fit)
x_fit = mean0 * (1 + cv0 * phi_fit) * (1 + rng.normal(0, 0.01, 40))

r = fit_curve(x_fit, p_fit)
out["fit_free"] = {"x": x_fit.tolist(), "p": p_fit.tolist(),
                   "cs": r["cs"], "mean": r["mean"], "cv": r["cv"],
                   "obj_norm": r["obj_norm"], "rmse": r["rmse"],
                   "max_rel": r["max_rel"], "r2": r["r2"]}
r2 = fit_curve(x_fit, p_fit, mean_fixed=mean0)
out["fit_fixed"] = {"cs": r2["cs"], "mean": r2["mean"], "cv": r2["cv"],
                    "obj_norm": r2["obj_norm"]}

# 6) 矩法统计
out["stats"] = {"x": [float(v) for v in x_fit],
                "mean": float(sample_stats(x_fit)["mean"]),
                "cv": float(sample_stats(x_fit)["cv"]),
                "cs": float(sample_stats(x_fit)["cs"])}

# 7) 经验频率
xs, pp = empirical_frequency(x_fit)
out["emp"] = {"xs": xs.tolist(), "p": pp.tolist()}

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "test_vectors.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("test_vectors.json 已生成（phi 网格 %d×%d，gammaincc %d×%d）"
      % (len(cs_p), len(p_p), len(a_g), len(x_g)))
