# -*- coding: utf-8 -*-
"""
P-Ⅲ型频率曲线适线分析软件（水文频率计算）
================================================
功能：
  1. 输入多年水文资料（年最大洪峰流量 / 年降水量 / 年径流量等），
     自动计算均值 x̄、变差系数 Cv、偏态系数 Cs（矩法初值）；
  2. 在海森机率格纸上绘制经验频率点与 P-Ⅲ 理论频率曲线；
  3. 可手动调节 x̄、Cv、Cs（滑杆 + 精确输入），实时适线；
  4. 一键优化适线（加权最小二乘，均值可固定为样本均值）；
  5. 输出设计成果表（频率—重现期—离均系数—设计值）、
     导出图片（PNG/SVG/PDF）、导出 CSV、生成计算成果报告；
  6. 工程保存/载入（JSON）。

运行方式：
  python app.py             （源码运行）
  或打包后的 P3频率计算软件.exe
"""

import ctypes
import datetime
import json
import os
import re
import sys

import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk

import numpy as np

import matplotlib
matplotlib.use("TkAgg")
from matplotlib import rcParams
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from p3core import (STD_P, PAPER_MAJOR, PAPER_MINOR, curve_sample_u,
                    empirical_frequency, fit_curve, frequency_curve,
                    norm_ppf, pct_to_x, phi_pearson3, sample_stats)

# ---------------------------------------------------------------------------
# 全局设置
# ---------------------------------------------------------------------------
APP_NAME = "P-Ⅲ型频率曲线适线分析软件"
APP_VERSION = "1.0.0"

rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                               "Noto Sans CJK SC", "Arial Unicode MS", "DejaVu Sans"]
rcParams["axes.unicode_minus"] = False

try:  # Windows 高 DPI 支持
    ctypes.windll.shcore.SetProcessDpiAwareness(1)
except Exception:
    pass

# 示例数据：某水文站 1991–2020 年最大洪峰流量（m³/s）
SAMPLE_YEARS = list(range(1991, 2021))
SAMPLE_DATA = [1870, 1560, 2140, 1320, 2480, 1780, 3020, 1650, 1980, 2250,
               1430, 2680, 1850, 2090, 1210, 2340, 1760, 2950, 1520, 2460,
               1890, 2150, 1390, 2580, 1710, 2230, 1940, 1620, 2810, 2010]


def set_dpi_aware():
    try:
        ctypes.windll.shcore.SetProcessDpiAwareness(1)
    except Exception:
        pass


# ---------------------------------------------------------------------------
# 数据解析
# ---------------------------------------------------------------------------
def parse_data_text(text):
    """解析输入文本：每行一个数值，或 '年份 数值' 两列。
    返回 (years, values)；years 为 None 表示没有年份信息。"""
    rows = []
    for ln in text.splitlines():
        if not ln.strip():
            continue
        toks = re.split(r"[\s,，;；\t]+", ln.strip().replace("－", "-"))
        nums = []
        for t in toks:
            try:
                nums.append(float(t))
            except ValueError:
                pass
        if nums:
            rows.append(nums)
    if not rows:
        raise ValueError("未识别到任何数值，请检查输入格式")
    two_col = all(len(r) == 2 for r in rows)
    if two_col and all(1700 <= r[0] <= 2200 and r[0].is_integer() for r in rows):
        years = [int(r[0]) for r in rows]
        vals = [r[1] for r in rows]
    else:
        years = None
        vals = [n for r in rows for n in r]
    return years, np.asarray(vals, dtype=float)


def load_data_file(path):
    """从 txt/csv/xlsx/xls 读取数据，返回 (years, values, label)。"""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xlsm"):
        import openpyxl
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        ws = wb.worksheets[0]
        rows = []
        for row in ws.iter_rows(values_only=True):
            nums = [float(c) for c in row if isinstance(c, (int, float))]
            if nums:
                rows.append(nums)
        wb.close()
        if not rows:
            raise ValueError("Excel 文件中未找到数值数据")
        two_col = all(len(r) == 2 for r in rows)
        if two_col and all(1700 <= r[0] <= 2200 and float(r[0]).is_integer() for r in rows):
            years = [int(r[0]) for r in rows]
            vals = [r[1] for r in rows]
        else:
            years = None
            vals = [n for r in rows for n in r]
        return years, np.asarray(vals, float), os.path.basename(path)
    if ext == ".xls":
        try:
            import xlrd
        except ImportError:
            raise ValueError("读取旧版 .xls 文件需要 xlrd 组件，请将文件另存为 "
                             ".xlsx 格式后再导入")
        wb = xlrd.open_workbook(path)
        ws = wb.sheet_by_index(0)
        rows = []
        for r in range(ws.nrows):
            nums = [c for c in ws.row_values(r)
                    if isinstance(c, (int, float)) and not isinstance(c, bool)]
            if nums:
                rows.append(nums)
        if not rows:
            raise ValueError("Excel 文件中未找到数值数据")
        two_col = all(len(r) == 2 for r in rows)
        if two_col and all(1700 <= r[0] <= 2200 and float(r[0]).is_integer() for r in rows):
            years = [int(r[0]) for r in rows]
            vals = [r[1] for r in rows]
        else:
            years = None
            vals = [n for r in rows for n in r]
        return years, np.asarray(vals, float), os.path.basename(path)
    # 文本文件
    raw = None
    for enc in ("utf-8-sig", "gbk", "utf-16", "latin-1"):
        try:
            with open(path, "r", encoding=enc) as f:
                raw = f.read()
            break
        except (UnicodeDecodeError, UnicodeError):
            continue
    if raw is None:
        raise ValueError("无法识别文件编码")
    years, vals = parse_data_text(raw)
    return years, vals, os.path.basename(path)


# ---------------------------------------------------------------------------
# 主程序
# ---------------------------------------------------------------------------
class P3App:
    def __init__(self, root):
        self.root = root
        self.x = None          # 数据数组
        self.years = None      # 年份列表或 None
        self.source = ""       # 数据来源说明
        self.var_name = "年最大洪峰流量"
        self.unit = "m³/s"
        self._syncing = False

        root.title(f"{APP_NAME}  v{APP_VERSION}")
        root.geometry("1320x860")
        root.minsize(1150, 720)

        self._build_style()
        self._build_menu()
        self._build_ui()
        self.load_sample_data()
        self.set_status("就绪。可直接输入数据，或点击【示例数据】【打开文件…】")

    # ------------------------------------------------------------------ UI
    def _build_style(self):
        style = ttk.Style()
        try:
            if "vista" in style.theme_names():
                style.theme_use("vista")
        except Exception:
            pass
        style.configure(".", font=("Microsoft YaHei UI", 9))
        style.configure("Treeview", font=("Microsoft YaHei UI", 9), rowheight=22)
        style.configure("Treeview.Heading", font=("Microsoft YaHei UI", 9, "bold"))
        style.configure("TLabelframe.Label", font=("Microsoft YaHei UI", 9, "bold"))

    def _build_menu(self):
        mb = tk.Menu(self.root)
        m_file = tk.Menu(mb, tearoff=0)
        m_file.add_command(label="打开数据文件… (Ctrl+O)", command=self.open_file)
        m_file.add_separator()
        m_file.add_command(label="保存工程… (Ctrl+S)", command=self.save_project)
        m_file.add_command(label="载入工程…", command=self.load_project)
        m_file.add_separator()
        m_file.add_command(label="退出", command=self.root.destroy)
        mb.add_cascade(label="文件", menu=m_file)

        m_data = tk.Menu(mb, tearoff=0)
        m_data.add_command(label="示例数据", command=self.load_sample_data)
        m_data.add_command(label="粘贴自剪贴板", command=self.paste_data)
        m_data.add_command(label="清空数据", command=self.clear_data)
        mb.add_cascade(label="数据", menu=m_data)

        m_ana = tk.Menu(mb, tearoff=0)
        m_ana.add_command(label="矩法估计参数", command=self.use_moment_params)
        m_ana.add_command(label="优化适线（最小二乘）", command=self.auto_fit)
        m_ana.add_command(label="重新计算 (F5)", command=self.recompute)
        mb.add_cascade(label="分析", menu=m_ana)

        m_exp = tk.Menu(mb, tearoff=0)
        m_exp.add_command(label="导出图片 (PNG, 300dpi)…", command=self.export_png)
        m_exp.add_command(label="导出图片 (SVG/PDF)…", command=self.export_vector)
        m_exp.add_command(label="导出成果表 (CSV)…", command=self.export_csv)
        m_exp.add_command(label="生成计算成果报告 (TXT)…", command=self.export_report)
        mb.add_cascade(label="导出", menu=m_exp)

        m_help = tk.Menu(mb, tearoff=0)
        m_help.add_command(label="使用说明", command=self.show_help)
        m_help.add_command(label="关于", command=self.show_about)
        mb.add_cascade(label="帮助", menu=m_help)
        self.root.config(menu=mb)

        self.root.bind("<Control-o>", lambda e: self.open_file())
        self.root.bind("<Control-s>", lambda e: self.save_project())
        self.root.bind("<F5>", lambda e: self.recompute())

    def _build_ui(self):
        paned = ttk.Panedwindow(self.root, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True)

        # ------------------------- 左侧控制面板 -------------------------
        left = ttk.Frame(paned, padding=(6, 6, 6, 6))
        paned.add(left, weight=0)

        # -- 数据输入 --
        f_data = ttk.LabelFrame(left, text="① 数据输入（每行一个年值，或“年份 数值”两列）", padding=4)
        f_data.pack(fill=tk.X, pady=(0, 6))
        self.txt = scrolledtext.ScrolledText(f_data, height=8, font=("Consolas", 9),
                                             wrap=tk.NONE)
        self.txt.pack(fill=tk.X)
        self.txt.bind("<Control-Return>", lambda e: self.recompute())
        btns = ttk.Frame(f_data)
        btns.pack(fill=tk.X, pady=(4, 0))
        ttk.Button(btns, text="示例数据", width=9, command=self.load_sample_data).pack(side=tk.LEFT)
        ttk.Button(btns, text="打开文件…", width=9, command=self.open_file).pack(side=tk.LEFT, padx=4)
        ttk.Button(btns, text="粘贴", width=7, command=self.paste_data).pack(side=tk.LEFT)
        ttk.Button(btns, text="清空", width=7, command=self.clear_data).pack(side=tk.LEFT, padx=4)
        ttk.Button(btns, text="重新计算 (F5)", command=self.recompute).pack(side=tk.RIGHT)
        self.lbl_source = ttk.Label(f_data, text="", foreground="#555555")
        self.lbl_source.pack(anchor=tk.W, pady=(3, 0))

        # -- 统计量 --
        f_stat = ttk.LabelFrame(left, text="② 样本统计量（矩法初值）", padding=4)
        f_stat.pack(fill=tk.X, pady=(0, 6))
        self.lbl_stat = ttk.Label(f_stat, text="尚无数据", justify=tk.LEFT)
        self.lbl_stat.pack(anchor=tk.W)

        # -- 参数调节 --
        f_par = ttk.LabelFrame(left, text="③ 频率曲线参数（适线：拖动滑杆或输入精确值）", padding=4)
        f_par.pack(fill=tk.X, pady=(0, 6))
        self.mean_var = tk.DoubleVar(value=1000.0)
        self.cv_var = tk.DoubleVar(value=0.3)
        self.cs_var = tk.DoubleVar(value=1.0)
        self.mean_entry_var = tk.StringVar()
        self.cv_entry_var = tk.StringVar()
        self.cs_entry_var = tk.StringVar()

        self._param_row(f_par, 0, "均值", self.mean_var, self.mean_entry_var, "mean")
        self._param_row(f_par, 1, "变差系数 Cv", self.cv_var, self.cv_entry_var, "cv")
        self._param_row(f_par, 2, "偏态系数 Cs", self.cs_var, self.cs_entry_var, "cs")

        self.chk_fix_mean = tk.BooleanVar(value=True)
        self.chk_points = tk.BooleanVar(value=True)
        self.chk_years = tk.BooleanVar(value=False)
        row = ttk.Frame(f_par)
        row.pack(fill=tk.X, pady=(4, 0))
        ttk.Checkbutton(row, text="均值固定为样本均值",
                        variable=self.chk_fix_mean,
                        command=self._on_fix_mean).pack(side=tk.LEFT)
        ttk.Checkbutton(row, text="显示经验点", variable=self.chk_points,
                        command=self.redraw).pack(side=tk.LEFT, padx=8)
        ttk.Checkbutton(row, text="显示年份标签", variable=self.chk_years,
                        command=self.redraw).pack(side=tk.LEFT, padx=8)
        row2 = ttk.Frame(f_par)
        row2.pack(fill=tk.X, pady=(5, 0))
        ttk.Button(row2, text="矩法估计", command=self.use_moment_params).pack(side=tk.LEFT)
        ttk.Button(row2, text="优化适线（最小二乘）", command=self.auto_fit).pack(side=tk.LEFT, padx=6)
        self.lbl_fit = ttk.Label(row2, text="", foreground="#0066aa")
        self.lbl_fit.pack(side=tk.LEFT)

        # -- 成果表 --
        f_tab = ttk.LabelFrame(left, text="④ 设计成果表（P-Ⅲ 理论频率曲线）", padding=4)
        f_tab.pack(fill=tk.BOTH, expand=True)
        cols = ("p", "t", "phi", "xp")
        self.tree = ttk.Treeview(f_tab, columns=cols, show="headings", height=11)
        heads = {"p": ("频率P(%)", 70, "center"),
                 "t": ("重现期T(年)", 90, "center"),
                 "phi": ("离均系数Φ", 90, "center"),
                 "xp": ("设计值", 120, "e")}
        for c, (txt, w, anc) in heads.items():
            self.tree.heading(c, text=txt)
            self.tree.column(c, width=w, anchor=anc)
        self.tree.tag_configure("design", background="#fff7d6")
        sb = ttk.Scrollbar(f_tab, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=sb.set)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.lbl_err = ttk.Label(left, text="", foreground="#aa3300")
        self.lbl_err.pack(anchor=tk.W, pady=(4, 0))

        # ------------------------- 右侧绘图区 -------------------------
        right = ttk.Frame(paned, padding=(0, 6, 6, 6))
        paned.add(right, weight=1)
        self.fig = Figure(figsize=(7.4, 6.2), dpi=100, facecolor="white")
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        tb = NavigationToolbar2Tk(self.canvas, right)
        tb.update()

        # ------------------------- 状态栏 -------------------------
        self.status = ttk.Label(self.root, text="", relief=tk.SUNKEN, anchor=tk.W,
                                padding=(6, 2))
        self.status.pack(fill=tk.X, side=tk.BOTTOM)

    def _param_row(self, parent, row, label, var, entry_var, key):
        f = ttk.Frame(parent)
        f.pack(fill=tk.X, pady=2)
        ttk.Label(f, text=label, width=11).pack(side=tk.LEFT)
        e = ttk.Entry(f, textvariable=entry_var, width=11)
        e.pack(side=tk.LEFT)
        e.bind("<Return>", lambda ev, k=key: self._entry_apply(k))
        e.bind("<FocusOut>", lambda ev, k=key: self._entry_apply(k))
        scale = tk.Scale(f, from_=0, to=1, orient=tk.HORIZONTAL, resolution=0.001,
                         showvalue=False, length=230,
                         command=lambda v, k=key: self._scale_move(k))
        scale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(8, 0))
        setattr(self, f"entry_{key}", e)
        setattr(self, f"scale_{key}", scale)

    # ------------------------------------------------------------- 行为
    def set_status(self, msg):
        self.status.config(text=msg)

    def _sync_entries(self):
        self._syncing = True
        self.mean_entry_var.set(f"{self.mean_var.get():.3g}")
        self.cv_entry_var.set(f"{self.cv_var.get():.4g}")
        self.cs_entry_var.set(f"{self.cs_var.get():.4g}")
        self._syncing = False

    def _scale_move(self, key):
        if self._syncing:
            return
        self._sync_entries()
        self.redraw()

    def _entry_apply(self, key):
        if self._syncing:
            return
        try:
            v = float({"mean": self.mean_entry_var,
                       "cv": self.cv_entry_var,
                       "cs": self.cs_entry_var}[key].get())
        except ValueError:
            self._sync_entries()
            return
        if key == "mean" and v <= 0:
            messagebox.showwarning("参数无效", "均值必须大于 0")
            self._sync_entries()
            return
        if key == "cv" and not (0 <= v < 50):
            messagebox.showwarning("参数无效", "Cv 应为 0 ~ 50 之间的非负数")
            self._sync_entries()
            return
        var = {"mean": self.mean_var, "cv": self.cv_var, "cs": self.cs_var}[key]
        var.set(v)
        self._update_scale_ranges()
        self.redraw()

    def _on_fix_mean(self):
        state = tk.DISABLED if self.chk_fix_mean.get() else tk.NORMAL
        self.entry_mean.config(state=state)
        self.scale_mean.config(state=state)
        self.redraw()

    def _update_scale_ranges(self):
        mean = self.mean_var.get()
        self.scale_mean.config(from_=max(mean * 0.2, mean - 3 * mean),
                               to=max(mean * 3.0, mean + 1e-9),
                               resolution=max(mean / 400.0, 1e-6))
        self.scale_cv.config(from_=0.0, to=3.0, resolution=0.001)
        self.scale_cs.config(from_=-3.0, to=10.0, resolution=0.01)
        # 保证滑块取值在范围内
        for var, sc in ((self.mean_var, self.scale_mean),
                        (self.cv_var, self.scale_cv),
                        (self.cs_var, self.scale_cs)):
            lo, hi = sc.cget("from"), sc.cget("to")
            if var.get() < lo:
                var.set(lo)
            if var.get() > hi:
                var.set(hi)

    # ------------------------------------------------------------- 数据
    def load_sample_data(self):
        self.years = SAMPLE_YEARS
        self.x = np.asarray(SAMPLE_DATA, float)
        self.source = "示例数据：某水文站 1991–2020 年最大洪峰流量"
        self.txt.delete("1.0", tk.END)
        lines = "\n".join(f"{y} {v}" for y, v in zip(SAMPLE_YEARS, SAMPLE_DATA))
        self.txt.insert("1.0", lines)
        self._on_new_data()

    def open_file(self):
        path = filedialog.askopenfilename(
            title="打开水文数据文件",
            filetypes=[("数据文件", "*.txt *.csv *.xlsx *.xls"),
                       ("Excel 工作簿", "*.xlsx *.xls"),
                       ("文本文件", "*.txt *.csv"), ("所有文件", "*.*")])
        if not path:
            return
        try:
            years, vals, label = load_data_file(path)
        except Exception as exc:
            messagebox.showerror("读取失败", f"无法读取文件：\n{exc}")
            return
        if vals.size < 3:
            messagebox.showerror("数据不足", "至少需要 3 个数据点")
            return
        self.years = years
        self.x = vals
        self.source = label
        self.txt.delete("1.0", tk.END)
        if years:
            self.txt.insert("1.0", "\n".join(f"{y} {v}" for y, v in zip(years, vals)))
        else:
            self.txt.insert("1.0", "\n".join(f"{v:g}" for v in vals))
        self._on_new_data()
        self.set_status(f"已读取数据文件：{label}（{vals.size} 个数据）")

    def paste_data(self):
        try:
            text = self.root.clipboard_get()
        except tk.TclError:
            messagebox.showinfo("剪贴板为空", "剪贴板中没有文本内容")
            return
        self.txt.delete("1.0", tk.END)
        self.txt.insert("1.0", text)
        self.recompute()

    def clear_data(self):
        self.txt.delete("1.0", tk.END)
        self.x = None
        self.years = None
        self.source = ""
        self.lbl_source.config(text="")
        self.lbl_stat.config(text="尚无数据")
        self.tree.delete(*self.tree.get_children())
        self.lbl_err.config(text="")
        self.lbl_fit.config(text="")
        self._draw_empty()
        self.set_status("数据已清空")

    def _on_new_data(self):
        n = self.x.size
        self.lbl_source.config(text=f"数据来源：{self.source}  （n = {n} 年）")
        self._update_stats_label()
        # 矩法初值
        s = sample_stats(self.x)
        self.mean_var.set(s["mean"])
        self.cv_var.set(s["cv"])
        self.cs_var.set(max(s["cs"], 0.05) if s["cs"] > 0 else min(s["cs"], -0.05))
        self._update_scale_ranges()
        self._on_fix_mean()
        self._sync_entries()
        self.redraw()
        self.set_status(f"已载入 {n} 个数据，矩法初值：均值={s['mean']:.2f}，"
                        f"Cv={s['cv']:.3f}，Cs={s['cs']:.3f}。"
                        f"可拖动 Cv/Cs 滑杆人工适线，或点击【优化适线】")

    def _update_stats_label(self):
        s = sample_stats(self.x)
        txt = (f"n = {s['n']} 年    均值 = {s['mean']:,.2f}    "
               f"标准差 σ = {s['std']:,.2f}\n"
               f"矩法：Cv = {s['cv']:.4f}    Cs = {s['cs']:.4f}\n"
               f"最大 = {s['xmax']:,.2f}（{self._year_of(s['xmax'])}）    "
               f"最小 = {s['xmin']:,.2f}（{self._year_of(s['xmin'])}）")
        self.lbl_stat.config(text=txt)

    def _year_of(self, val):
        if self.years is None or self.x is None:
            return "—"
        i = int(np.argmin(np.abs(self.x - val)))
        return str(self.years[i])

    # ------------------------------------------------------------- 分析
    def recompute(self):
        text = self.txt.get("1.0", tk.END)
        try:
            years, vals = parse_data_text(text)
        except ValueError as exc:
            messagebox.showwarning("数据格式错误", str(exc))
            return
        if vals.size < 3:
            messagebox.showwarning("数据不足", "至少需要 3 个数据点")
            return
        self.years = years
        self.x = vals
        self.source = "手动输入"
        self._on_new_data()
        self.set_status("已按输入数据重新计算")

    def use_moment_params(self):
        if self.x is None:
            messagebox.showinfo("提示", "请先输入数据")
            return
        s = sample_stats(self.x)
        self.mean_var.set(s["mean"])
        self.cv_var.set(s["cv"])
        self.cs_var.set(s["cs"])
        self._update_scale_ranges()
        self._sync_entries()
        self.redraw()
        self.set_status(f"已采用矩法估计：均值={s['mean']:.2f}，Cv={s['cv']:.4f}，Cs={s['cs']:.4f}")

    def auto_fit(self):
        if self.x is None:
            messagebox.showinfo("提示", "请先输入数据")
            return
        xs, p = empirical_frequency(self.x)
        if xs.std() < 1e-12:   # 常数序列的退化情况
            self.cv_var.set(0.0)
            self.cs_var.set(1.0)
            self._update_scale_ranges()
            self._sync_entries()
            self.redraw()
            self.set_status("数据为常数序列，Cv 已置 0（曲线为水平线）")
            return
        self.root.config(cursor="watch")
        self.root.update_idletasks()
        try:
            mean_fixed = self.mean_var.get() if self.chk_fix_mean.get() else None
            r = fit_curve(xs, p / 100.0, mean_fixed=mean_fixed)
        finally:
            self.root.config(cursor="")
        if r is None:
            messagebox.showerror("适线失败", "未能找到合适的参数组合。\n"
                                 "请确认数据均为正值且存在变幅。")
            return
        self.mean_var.set(r["mean"])
        self.cv_var.set(r["cv"])
        self.cs_var.set(r["cs"])
        self._update_scale_ranges()
        self._sync_entries()
        self.redraw()
        self.lbl_fit.config(
            text=f"  最优：Cv={r['cv']:.3f} Cs={r['cs']:.3f} "
                 f"相对SSE={r['obj_norm']:.2e}")
        self.set_status(f"优化适线完成：均值={r['mean']:.2f}，Cv={r['cv']:.4f}，"
                        f"Cs={r['cs']:.4f}，RMSE={r['rmse']:.2f}，"
                        f"最大相对误差={r['max_rel']*100:.2f}%")

    # ------------------------------------------------------------- 绘图
    def _draw_empty(self):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        ax.text(0.5, 0.5, "请先输入数据（每行一个年值），或点击【示例数据】【打开文件…】",
                ha="center", va="center", fontsize=12, color="#888888",
                transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        self.canvas.draw_idle()

    def redraw(self):
        if self.x is None:
            return
        mean = self.mean_var.get()
        cv = self.cv_var.get()
        cs = self.cs_var.get()
        xs, p_emp = empirical_frequency(self.x)

        self.fig.clear()
        ax = self.fig.add_subplot(111)

        # 海森机率格纸
        for p in PAPER_MINOR:
            ax.axvline(pct_to_x(p), color="#d9d9d9", lw=0.5, zorder=0)
        for p in PAPER_MAJOR:
            ax.axvline(pct_to_x(p), color="#9a9a9a", lw=0.8, zorder=0)
        xt = pct_to_x(np.asarray(PAPER_MAJOR, float))
        ax.set_xticks(xt)
        ax.set_xticklabels([f"{p:g}" for p in PAPER_MAJOR], fontsize=8)
        x_lo, x_hi = pct_to_x(0.01) - 0.12, pct_to_x(99.9) + 0.18
        ax.set_xlim(x_lo, x_hi)

        # 上方重现期轴
        t_ticks = [10000, 1000, 200, 100, 50, 20, 10, 5, 2]
        ax_t = ax.twiny()
        ax_t.set_xticks(pct_to_x(100.0 / np.asarray(t_ticks, float)))
        ax_t.set_xticklabels([f"{t:g}" for t in t_ticks], fontsize=8)
        ax_t.set_xlim(x_lo, x_hi)
        ax_t.set_xlabel("重现期 T（年）", fontsize=9)
        ax_t.grid(False)

        # 理论曲线
        u, y = curve_sample_u(mean, cv, cs)
        ax.plot(u, y, color="#d62728", lw=2.0, zorder=2,
                label=fr"P-Ⅲ理论曲线  $\bar{{x}}$={mean:.3g}  "
                      fr"Cv={cv:.4g}  Cs={cs:.4g}")

        # 经验点
        if self.chk_points.get():
            ax.scatter(pct_to_x(p_emp), xs, s=30, facecolors="none",
                       edgecolors="#1f77b4", linewidths=1.2, zorder=3,
                       label="经验点（P=m/(n+1)）")
            if self.chk_years.get() and self.years is not None:
                for yy, pp, vv in zip(self.years, p_emp, xs):
                    ax.annotate(str(yy), (pct_to_x(pp), vv),
                                textcoords="offset points", xytext=(0, 6),
                                fontsize=6.5, ha="center", color="#1f77b4",
                                zorder=4)

        # 坐标轴与图例
        ax.set_xlabel("频率 P（%）", fontsize=10)
        ax.set_ylabel(f"{self.var_name}  {self.unit}", fontsize=10)
        ax.set_title(f"{self.var_name} P-Ⅲ型频率曲线"
                     f"（n={self.x.size} 年，适线法）", fontsize=12)
        ax.legend(loc="lower right", fontsize=8.5, framealpha=0.9)
        ax.tick_params(axis="y", labelsize=9)

        # 纵轴范围
        mid = (u > pct_to_x(0.05)) & (u < pct_to_x(99.5))
        yc = y[mid]
        if yc.size:
            ylo = min(xs.min(), yc.min())
            yhi = max(xs.max(), yc.max())
        else:
            ylo, yhi = xs.min(), xs.max()
        span = (yhi - ylo) or abs(yhi) or 1.0
        ax.set_ylim(ylo - 0.1 * span, yhi + 0.06 * span)

        self.fig.subplots_adjust(left=0.10, right=0.965, top=0.89, bottom=0.105)
        self.canvas.draw_idle()
        self._update_table(mean, cv, cs)
        self._update_err(xs, p_emp, mean, cv, cs)

    def _update_table(self, mean, cv, cs):
        self.tree.delete(*self.tree.get_children())
        phi = phi_pearson3(cs, np.asarray(STD_P, float) / 100.0)
        xp = mean * (1.0 + cv * phi)
        for p, ph, x in zip(STD_P, phi, xp):
            t = 100.0 / p
            tag = "design" if p <= 1.0 else ""
            self.tree.insert("", tk.END, values=(
                f"{p:g}", f"{t:g}", f"{ph:.4f}", f"{x:,.2f}"), tags=(tag,))

    def _update_err(self, xs, p_emp, mean, cv, cs):
        phi = phi_pearson3(cs, p_emp / 100.0)
        xhat = mean * (1.0 + cv * phi)
        resid = xs - xhat
        rel = resid / xs
        rmse = float(np.sqrt((resid ** 2).mean()))
        sse = float((rel ** 2).sum())
        max_rel = float(np.abs(rel).max())
        r2 = float(np.corrcoef(xs, xhat)[0, 1] ** 2)
        self.lbl_err.config(
            text=f"拟合检验：相对残差平方和 Σ[(x-x̂)/x]² = {sse:.4f}    "
                 f"RMSE = {rmse:.2f}   最大相对误差 = {max_rel*100:.2f}%   "
                 f"R² = {r2:.4f}")

    # ------------------------------------------------------------- 导出
    def _save_figure(self, path, dpi=300):
        self.fig.savefig(path, dpi=dpi, facecolor="white", bbox_inches="tight")

    def export_png(self):
        if self.x is None:
            messagebox.showinfo("提示", "请先输入数据")
            return
        path = filedialog.asksaveasfilename(
            title="导出图片", defaultextension=".png",
            filetypes=[("PNG 图片", "*.png")],
            initialfile="P3频率曲线.png")
        if not path:
            return
        self._save_figure(path)
        self.set_status(f"图片已导出：{path}")

    def export_vector(self):
        if self.x is None:
            messagebox.showinfo("提示", "请先输入数据")
            return
        path = filedialog.asksaveasfilename(
            title="导出矢量图", defaultextension=".svg",
            filetypes=[("SVG 矢量图", "*.svg"), ("PDF 文档", "*.pdf")],
            initialfile="P3频率曲线.svg")
        if not path:
            return
        self._save_figure(path, dpi=300)
        self.set_status(f"图片已导出：{path}")

    def _result_rows(self):
        mean = self.mean_var.get()
        cv = self.cv_var.get()
        cs = self.cs_var.get()
        phi = phi_pearson3(cs, np.asarray(STD_P, float) / 100.0)
        xp = mean * (1.0 + cv * phi)
        rows = []
        for p, ph, x in zip(STD_P, phi, xp):
            rows.append((p, 100.0 / p, ph, x))
        return rows

    def export_csv(self):
        if self.x is None:
            messagebox.showinfo("提示", "请先输入数据")
            return
        path = filedialog.asksaveasfilename(
            title="导出成果表", defaultextension=".csv",
            filetypes=[("CSV 文件", "*.csv")], initialfile="P3成果表.csv")
        if not path:
            return
        mean, cv, cs = (self.mean_var.get(), self.cv_var.get(), self.cs_var.get())
        with open(path, "w", encoding="utf-8-sig", newline="") as f:
            f.write(f"变量,{self.var_name}\n")
            f.write(f"单位,{self.unit}\n")
            f.write(f"样本数,{self.x.size}\n")
            f.write(f"均值,{mean:.6g}\n")
            f.write(f"变差系数Cv,{cv:.6g}\n")
            f.write(f"偏态系数Cs,{cs:.6g}\n")
            f.write("频率P(%),重现期T(年),离均系数Φ,设计值\n")
            for p, t, ph, x in self._result_rows():
                f.write(f"{p:g},{t:g},{ph:.6f},{x:.4f}\n")
        self.set_status(f"成果表已导出：{path}")

    def export_report(self):
        if self.x is None:
            messagebox.showinfo("提示", "请先输入数据")
            return
        path = filedialog.asksaveasfilename(
            title="导出成果报告", defaultextension=".txt",
            filetypes=[("文本文件", "*.txt")], initialfile="P3计算成果报告.txt")
        if not path:
            return
        mean, cv, cs = (self.mean_var.get(), self.cv_var.get(), self.cs_var.get())
        s = sample_stats(self.x)
        xs, p_emp = empirical_frequency(self.x)
        phi = phi_pearson3(cs, p_emp / 100.0)
        xhat = mean * (1.0 + cv * phi)
        rel = (xs - xhat) / xs
        lines = []
        lines.append("=" * 56)
        lines.append("  P-Ⅲ型频率曲线适线计算成果报告")
        lines.append("=" * 56)
        lines.append(f"计算时间：{datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
        lines.append(f"数据来源：{self.source or '手动输入'}")
        lines.append(f"变量名称：{self.var_name}（{self.unit}）")
        lines.append(f"样本数 n = {s['n']}")
        lines.append("-" * 56)
        lines.append(f"均值        x̄ = {s['mean']:.4f}")
        lines.append(f"样本标准差  σ = {s['std']:.4f}")
        lines.append(f"变差系数   Cv = {cv:.4f}")
        lines.append(f"偏态系数   Cs = {cs:.4f}    （Cs/Cv = {cs/cv:.3f}）" if cv > 0
                     else f"偏态系数   Cs = {cs:.4f}")
        lines.append(f"最大值 {s['xmax']:.2f}，最小值 {s['xmin']:.2f}")
        lines.append("-" * 56)
        lines.append("经验频率公式：P = m/(n+1) × 100%（数学期望公式）")
        lines.append("理论曲线：x_p = x̄(1 + Cv·Φ_p)，Φ_p 由皮尔逊Ⅲ型分布反解")
        lines.append(f"拟合检验：相对残差平方和 = {float((rel**2).sum()):.4f}，"
                     f"RMSE = {float(np.sqrt((rel**2).mean()))*mean:.2f}，"
                     f"最大相对误差 = {float(np.abs(rel).max())*100:.2f}%")
        lines.append("-" * 56)
        lines.append(f"{'频率P(%)':>10}{'重现期T(年)':>12}{'离均系数Φ':>12}{'设计值':>14}")
        lines.append("-" * 56)
        for p, t, ph, x in self._result_rows():
            lines.append(f"{p:>10g}{t:>12g}{ph:>12.4f}{x:>14.2f}")
        lines.append("-" * 56)
        lines.append("注：P-Ⅲ 密度 f(x)=β^α/Γ(α)·(x-a0)^(α-1)·e^(-β(x-a0))，")
        lines.append("    α=4/Cs²，β=2/(x̄·Cv·Cs)，a0=x̄(1-2Cv/Cs)。")
        with open(path, "w", encoding="utf-8-sig") as f:
            f.write("\n".join(lines))
        self.set_status(f"成果报告已导出：{path}")

    # ------------------------------------------------------------- 工程
    def save_project(self):
        if self.x is None:
            messagebox.showinfo("提示", "请先输入数据")
            return
        path = filedialog.asksaveasfilename(
            title="保存工程", defaultextension=".json",
            filetypes=[("工程文件", "*.json")], initialfile="P3工程.json")
        if not path:
            return
        data = {
            "version": APP_VERSION,
            "var_name": self.var_name,
            "unit": self.unit,
            "years": self.years,
            "values": self.x.tolist(),
            "source": self.source,
            "mean": self.mean_var.get(),
            "cv": self.cv_var.get(),
            "cs": self.cs_var.get(),
            "fix_mean": self.chk_fix_mean.get(),
            "show_points": self.chk_points.get(),
            "show_years": self.chk_years.get(),
        }
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        self.set_status(f"工程已保存：{path}")

    def load_project(self):
        path = filedialog.askopenfilename(
            title="载入工程", filetypes=[("工程文件", "*.json"), ("所有文件", "*.*")])
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
            vals = np.asarray(data["values"], float)
            if vals.size < 3:
                raise ValueError("数据点不足")
            self.years = data.get("years")
            self.x = vals
            self.source = data.get("source", os.path.basename(path))
            self.var_name = data.get("var_name", self.var_name)
            self.unit = data.get("unit", self.unit)
            self.mean_var.set(float(data.get("mean", sample_stats(vals)["mean"])))
            self.cv_var.set(float(data.get("cv", 0.3)))
            self.cs_var.set(float(data.get("cs", 1.0)))
            self.chk_fix_mean.set(bool(data.get("fix_mean", True)))
            self.chk_points.set(bool(data.get("show_points", True)))
            self.chk_years.set(bool(data.get("show_years", False)))
            self.txt.delete("1.0", tk.END)
            if self.years:
                self.txt.insert("1.0", "\n".join(
                    f"{y} {v}" for y, v in zip(self.years, vals)))
            else:
                self.txt.insert("1.0", "\n".join(f"{v:g}" for v in vals))
            self._update_scale_ranges()
            self._on_fix_mean()
            self._sync_entries()
            self._update_stats_label()
            self.redraw()
            self.set_status(f"已载入工程：{os.path.basename(path)}")
        except Exception as exc:
            messagebox.showerror("载入失败", f"无法载入工程文件：\n{exc}")

    # ------------------------------------------------------------- 帮助
    def show_help(self):
        text = (
            "【P-Ⅲ型频率曲线适线分析软件】使用说明\n"
            "────────────────────────────────────\n"
            "1. 数据输入\n"
            "   每行输入一个年值；也支持“年份 数值”两列格式；\n"
            "   或点击【打开文件…】导入 txt/csv/Excel 文件。\n"
            "2. 自动计算\n"
            "   软件自动用矩法估计均值、变差系数 Cv、偏态系数 Cs，\n"            "   并按 P=m/(n+1) 计算经验频率，在海森机率格纸上绘制。\n"
            "3. 适线（两种方式）\n"
            "   · 人工适线：拖动 Cv、Cs 滑杆（或输入精确值），曲线实时更新；\n"
            "   · 优化适线：点击【优化适线】，自动求加权最小二乘最优参数。\n"
            "   勾选“均值固定为样本均值”时仅调节 Cv、Cs（推荐）。\n"
            "4. 成果输出\n"
            "   · 左侧成果表给出各设计频率的离均系数与设计值；\n"
            "   · 导出图片（PNG 300dpi / SVG / PDF）；\n"
            "   · 导出 CSV 成果表、生成 TXT 计算成果报告；\n"
            "   · 保存/载入工程（json），可随时继续调整。\n"
            "5. 建议\n"
            "   样本数 n ≥ 20 时适线结果较稳定；Cs 一般取 Cv 的 1~4 倍。\n"
            "   快捷键：F5 重新计算，Ctrl+O 打开文件，Ctrl+S 保存工程。")
        messagebox.showinfo("使用说明", text)

    def show_about(self):
        text = (f"{APP_NAME}  v{APP_VERSION}\n\n"
                "水文频率分析专业工具\n"
                "皮尔逊Ⅲ型（Pearson Type III）分布适线法\n\n"
                "计算方法：\n"
                "  · 经验频率：P = m/(n+1) × 100%\n"
                "  · 离均系数：由不完全伽马函数 Q(α,t) = P 反解\n"
                "  · 优化适线：加权最小二乘\n\n"
                "完全离线运行，数据仅保存在本机。")
        messagebox.showinfo("关于", text)


def main():
    set_dpi_aware()
    # 打包为 --windowed 时标准流为 None，重定向避免崩溃
    if getattr(sys, "frozen", False):
        try:
            if sys.stdout is None:
                sys.stdout = open(os.devnull, "w", encoding="utf-8")
            if sys.stderr is None:
                sys.stderr = open(os.devnull, "w", encoding="utf-8")
        except Exception:
            pass
    try:
        _run_gui()
        return 0
    except Exception:
        import traceback
        log = os.path.join(os.environ.get("TEMP", "."), "P3频率计算软件_错误日志.txt")
        try:
            with open(log, "w", encoding="utf-8") as f:
                f.write(traceback.format_exc())
        except Exception:
            pass
        try:
            import tkinter.messagebox as mb
            mb.showerror("程序错误", f"程序发生错误：\n{traceback.format_exc(limit=3)}\n"
                                    f"详细信息已写入：\n{log}")
        except Exception:
            pass
        return 1


def _run_gui():
    root = tk.Tk()
    try:
        import base64
        from app_icon import ICON_PNG_B64
        _icon = tk.PhotoImage(data=base64.b64decode(ICON_PNG_B64))
        root.iconphoto(True, _icon)
    except Exception:
        pass
    app = P3App(root)
    root.mainloop()


if __name__ == "__main__":
    sys.exit(main())
