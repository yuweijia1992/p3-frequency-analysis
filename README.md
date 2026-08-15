# P-Ⅲ型频率曲线适线分析（水文频率计算）

水文学专业的皮尔逊Ⅲ型（Pearson Type III）频率曲线适线分析工具，
包含 **网页版** 与 **Windows 桌面版** 两个版本，算法完全一致、结果可互相印证。

- 输入多年水文资料（年最大洪峰流量、年降水量、年径流量等），自动按矩法
  估计均值、变差系数 Cv、偏态系数 Cs；
- 在海森机率格纸上绘制经验点（P=m/(n+1)）与 P-Ⅲ 理论频率曲线；
- 支持人工适线（实时调节参数）与优化适线（加权最小二乘）；
- 输出设计成果表（频率—重现期—离均系数Φ—设计值）与拟合检验指标。

## 目录结构

```
├── P3频率曲线网页版.html    ← 网页版（单文件，双击即用，完全离线）
├── p3core.js / p3.html       ← 网页版核心算法 / 界面源码
├── build_single.py           ← 网页版单文件构建脚本
├── test_*.js / gen_vectors.py ← 网页版测试与验证
├── desktop/                  ← Windows 桌面版（Python + tkinter + matplotlib）
│   ├── dist/P3频率计算软件.exe  ← 免安装可执行程序（38MB）
│   ├── app.py / p3core.py      ← 桌面版主程序 / 核心算法
│   ├── build_exe.bat           ← 一键打包 exe 脚本
│   └── ...
└── README.md
```

## 使用方法

**网页版**：双击 `P3频率曲线网页版.html` → 粘贴数据 → 按 F5 计算 → 拖动滑杆适线。

**桌面版**：运行 `desktop/dist/P3频率计算软件.exe`（免安装）；
或 `pip install -r desktop/requirements.txt` 后 `python desktop/app.py` 从源码运行。

详见各版本目录内 README 与《使用说明.txt》。

## 计算方法

- 经验频率：P = m/(n+1) × 100%（数学期望公式）
- 理论曲线：x_p = x̄(1 + Cv·Φ_p)，离均系数 Φ_p 由 P-Ⅲ 分布与不完全伽马函数
  的关系 Q(4/Cs², 2Φ/Cs + 4/Cs²) = P 反解（α=4/Cs²，β=2/(x̄·Cv·Cs)）
- 矩法初值：Cv=√(Σ(Ki-1)²/n)，Cs=Σ(Ki-1)³/(n·Cv³)
- 优化适线：Cs 网格搜索 + 加权最小二乘（均值可固定为样本均值）

数值实现（Python 与 JavaScript 双语言）经与 SciPy 交叉验证，
精度达 1e-13 量级，不依赖任何外部计算服务。

## 测试

```bat
:: 网页版（需 Node.js）
npm install jsdom
node test_p3core.js      :: 核心算法交叉验证（与 Python 版比对）
node test_ui.js          :: 界面冒烟测试
node test_geometry.js    :: 绘图几何测试

:: 桌面版（需 Python）
python desktop/test_p3core.py      :: 核心模块单元测试
python desktop/test_vs_scipy.py    :: 与 SciPy 交叉验证
python desktop/smoke_test.py       :: GUI 冒烟测试
```

## 版本

v1.0.0
