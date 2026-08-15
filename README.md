# P-Ⅲ型频率曲线适线分析（网页版）

轻量网页版水文频率计算工具：**粘贴数据 → 适线 → 绘图**，完全离线、无任何外部依赖。

## 使用

- 直接双击打开 **`P3频率曲线网页版.html`**（单文件，Chrome/Edge/Firefox 均可，离线可用）；
- 把多年数据粘贴到左侧文本框（每行一个年值，或“年份 数值”两列），按 `F5` 或点击【计算】；
- 拖动 **Cv / Cs / 均值** 滑杆（或输入精确值）人工适线，曲线实时更新；
- 点击【优化适线（最小二乘）】自动求最优参数（勾选"均值固定"则只调 Cv、Cs）；
- 左侧给出设计成果表（频率—重现期—离均系数Φ—设计值）与拟合检验指标；
- 鼠标悬停经验点可查看年份/实测值/理论值；右上角可下载 PNG 图片与 CSV 成果表。

## 文件结构

| 文件 | 说明 |
| --- | --- |
| `P3频率曲线网页版.html` | **交付文件**（单文件，内联全部算法与界面） |
| `p3core.js` | 核心算法（纯 JS：不完全伽马函数及其逆、离均系数、矩法、优化适线） |
| `p3.html` | 界面源文件（引用 p3core.js） |
| `build_single.py` | 构建脚本：把 p3core.js 内联进 p3.html 生成单文件版 |
| `gen_vectors.py` | 用已验证的 Python 版（P3FrequencyTool 目录）生成测试向量 |
| `test_p3core.js` | 核心算法 Node 交叉验证（与 Python 版比对，容差 1e-8） |
| `test_ui.js` | 界面冒烟测试（jsdom，需 `npm install jsdom`） |
| `test_geometry.js` | 绘图几何测试（记录型 canvas 桩检查坐标与刻度） |

## 开发与测试

```bat
python build_single.py        :: 修改 p3core.js / p3.html 后重新生成单文件版
npm install jsdom             :: 安装 UI 测试依赖（仅测试需要）
node test_p3core.js           :: 算法测试
node test_ui.js               :: UI 测试
node test_geometry.js         :: 绘图几何测试
```

## 数值精度

与桌面版（Python）同一套算法实现，离均系数 Φ 经与 SciPy 交叉验证，
两种语言实现间最大差异 < 1e-8（深尾处迭代终止噪声），满足水文计算要求。

## 版本

v1.0.0
