# -*- coding: utf-8 -*-
"""GUI 冒烟测试：实例化主程序、模拟参数调节、绘制并保存图片、检验成果表。"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import tkinter as tk

from app import P3App
from p3core import sample_stats


def main():
    root = tk.Tk()
    app = P3App(root)
    root.update()
    time.sleep(0.3)

    # 1) 默认示例数据已载入
    assert app.x is not None and app.x.size == 30, "示例数据未载入"
    s = sample_stats(app.x)
    print(f"示例数据 n={s['n']}  x̄={s['mean']:.2f}  Cv={s['cv']:.4f}  Cs={s['cs']:.4f}")

    # 2) 模拟人工适线：调节参数并重绘
    app.mean_var.set(2010.0)
    app.cv_var.set(0.45)
    app.cs_var.set(1.60)
    app._update_scale_ranges()
    app._sync_entries()
    app.redraw()
    root.update()
    time.sleep(0.2)

    # 3) 保存绘图结果（供人工查看）
    png = os.path.join(os.path.dirname(os.path.abspath(__file__)), "smoke_plot.png")
    app.fig.savefig(png, dpi=110, facecolor="white", bbox_inches="tight")
    print("绘图已保存:", png)

    # 4) 成果表
    rows = app.tree.get_children()
    assert len(rows) == len(__import__("p3core").STD_P), "成果表行数错误"
    print("成果表行数:", len(rows), "| 首行:", app.tree.item(rows[0], "values"))

    # 5) 优化适线
    t0 = time.time()
    app.auto_fit()
    dt = time.time() - t0
    root.update()
    print(f"优化适线耗时 {dt*1000:.0f} ms → Cv={app.cv_var.get():.4f} Cs={app.cs_var.get():.4f}")

    # 6) 拟合检验标签
    print("拟合标签:", app.lbl_err.cget("text")[:80], "...")

    # 7) 重新计算（手动输入数据）
    app.txt.delete("1.0", tk.END)
    app.txt.insert("1.0", "\n".join(f"{i} {v}" for i, v in
                                    enumerate(np.linspace(500, 900, 25), 1990)))
    app.recompute()
    root.update()
    assert app.x.size == 25, "重新计算失败"
    print("手动输入重算 OK (n=25)")

    root.destroy()
    print("GUI 冒烟测试全部通过 ✔")


if __name__ == "__main__":
    main()
