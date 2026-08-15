# -*- coding: utf-8 -*-
"""构建网页版单文件：把 p3core.js 内联进 p3.html。

输出两个文件：
  * index.html               —— GitHub Pages 部署名（仓库根目录跟踪此文件）
  * P3频率曲线网页版.html     —— 本地交付名（git 忽略，双击即用）
"""
import os

here = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(here, "p3core.js"), "r", encoding="utf-8") as f:
    core = f.read()
with open(os.path.join(here, "p3.html"), "r", encoding="utf-8") as f:
    html = f.read()

MARK = '<script src="p3core.js"></script>'
assert MARK in html, "p3.html 中未找到脚本引用占位符"

html = html.replace(MARK, "<script>\n/* ===== p3core.js（内联）===== */\n" + core + "\n</script>")

for name in ("index.html", "P3频率曲线网页版.html"):
    out = os.path.join(here, name)
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"已生成：{out}（{os.path.getsize(out)/1024:.1f} KB）")
