@echo off
rem ============================================
rem  P-Ⅲ频率计算软件 一键打包脚本
rem  需要: python + pip install pyinstaller
rem ============================================
cd /d "%~dp0"
python -m PyInstaller --noconfirm --clean --onefile --windowed ^
  --name "P3频率计算软件" ^
  --icon icon.ico ^
  --hidden-import openpyxl ^
  --hidden-import xlrd ^
  --hidden-import matplotlib.backends.backend_tkagg ^
  app.py
echo.
echo 打包完成，程序位于 dist\P3频率计算软件.exe
pause
