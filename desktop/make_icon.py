# -*- coding: utf-8 -*-
"""生成软件图标：icon.ico（exe 用）+ app_icon.py（窗口图标，base64 PNG）。"""
import base64
import io
import os

import numpy as np
from PIL import Image, ImageDraw, ImageFont

S = 256


def draw_icon(size):
    img = Image.new("RGB", (size, size), "#eaf3fb")
    d = ImageDraw.Draw(img)
    s = size / 256.0
    # 圆角边框
    d.rounded_rectangle([3*s, 3*s, size-3*s, size-3*s], radius=18*s,
                        outline="#2b6cb0", width=max(4, int(5*s)))
    # 频率格纸网格
    for i in range(9):
        x = 22*s + i*(size-44*s)/8
        d.line([x, 20*s, x, size-44*s], fill="#bcd7ee", width=max(1, int(1.5*s)))
    for i in range(5):
        y = 20*s + i*(size-64*s)/4
        d.line([22*s, y, size-22*s, y], fill="#bcd7ee", width=max(1, int(1.5*s)))
    # 坐标轴
    d.line([22*s, 20*s, 22*s, size-44*s], fill="#4a7fb5", width=max(2, int(3*s)))
    d.line([22*s, size-44*s, size-22*s, size-44*s], fill="#4a7fb5", width=max(2, int(3*s)))
    # P-Ⅲ 理论曲线（红色，先陡后缓的 S 形）
    u = np.linspace(-3.4, 3.0, 300)
    p = 0.5 * (1 + np.tanh(u / 2.2))            # 近似正态 CDF
    yv = (size - 44*s) - 18*s - (size - 82*s) * p ** 1.35   # 左高右低、左侧偏陡
    pts = [(22*s + (u + 3.4) / 6.4 * (size - 44*s), yv) for u, yv in zip(u, yv)]
    d.line(pts, fill="#d62728", width=max(3, int(4*s)), joint="curve")
    # 经验点（蓝色空心圆）
    for uu in (-2.4, -1.5, -0.7, 0.0, 0.9, 1.9, 2.6):
        x = 22*s + (uu + 3.4) / 6.4 * (size - 44*s)
        pp = 0.5 * (1 + np.tanh(uu / 2.2))
        y = (size - 44*s) - 18*s - (size - 82*s) * pp ** 1.35
        r = 4.5*s
        d.ellipse([x-r, y-r, x+r, y+r], outline="#1f77b4", width=max(2, int(2.5*s)))
    # 底部文字 "P-Ⅲ"
    try:
        font = ImageFont.truetype("C:/Windows/Fonts/msyhbd.ttc", int(34*s))
    except Exception:
        font = ImageFont.load_default()
    text = "P-Ⅲ"
    bb = d.textbbox((0, 0), text, font=font)
    tw, th = bb[2]-bb[0], bb[3]-bb[1]
    d.text(((size-tw)/2 - bb[0], size-40*s - bb[1]), text, fill="#1a365d", font=font)
    return img


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    sizes = [16, 24, 32, 48, 64, 128, 256]
    imgs = [draw_icon(sz) for sz in sizes]
    imgs[-1].save(os.path.join(here, "icon.ico"), format="ICO",
                  sizes=[(sz, sz) for sz in sizes], append_images=imgs[1:])
    # 窗口图标 PNG（64px）→ app_icon.py
    buf = io.BytesIO()
    draw_icon(64).save(buf, format="PNG")
    b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    with open(os.path.join(here, "app_icon.py"), "w", encoding="utf-8") as f:
        f.write("# -*- coding: utf-8 -*-\n# 窗口图标（64x64 PNG，base64 内嵌，便于单文件打包）\n")
        f.write(f'ICON_PNG_B64 = "{b64}"\n')
    print("已生成 icon.ico 与 app_icon.py")


if __name__ == "__main__":
    main()
