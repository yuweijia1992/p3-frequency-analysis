# -*- coding: utf-8 -*-
"""构建单文件网页版：把 p3core.js 内联进 p3.html，生成 P3频率曲线网页版.html。"""
import os

here = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(here, "p3core.js"), "r", encoding="utf-8") as f:
    core = f.read()
with open(os.path.join(here, "p3.html"), "r", encoding="utf-8") as f:
    html = f.read()

MARK = '<script src="p3core.js"></script>'
assert MARK in html, "p3.html 中未找到脚本引用占位符"

html = html.replace(MARK, "<script>\n/* ===== p3core.js（内联）===== */\n" + core + "\n</script>")

out = os.path.join(here, "P3频率曲线网页版.html")
with open(out, "w", encoding="utf-8") as f:
    f.write(html)

print(f"已生成单文件：{out}")
print(f"大小：{os.path.getsize(out)/1024:.1f} KB（p3core.js 内联 {len(core)/1024:.1f} KB）")
