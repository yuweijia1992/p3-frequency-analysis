/* ============================================================================
 * test_ui.js —— 网页版 UI 冒烟测试（jsdom + canvas 桩）
 * 运行：node test_ui.js
 * ==========================================================================*/
"use strict";

const fs = require("fs");
const path = require("path");
const { JSDOM } = require("jsdom");

const html = fs.readFileSync(path.join(__dirname, "P3频率曲线网页版.html"), "utf8");

/* Canvas 2D 桩：吸收所有绘图调用 */
const noop = () => {};
const ctxStub = new Proxy({}, {
    get: (t, p) => (p === "canvas" ? null
                    : p === "measureText" ? () => ({ width: 10 })
                    : noop),
    set: () => true,
});

const errors = [];
const dom = new JSDOM(html, {
    runScripts: "dangerously",
    pretendToBeVisual: true,
    beforeParse(window) {
        window.HTMLCanvasElement.prototype.getContext = () => ctxStub;
        window.HTMLCanvasElement.prototype.toDataURL = () => "data:image/png;base64,AA==";
        window.alert = msg => { console.log("  [alert]", msg); };
        window.addEventListener("error", ev => errors.push(ev.message));
    },
});

const doc = dom.window.document;
let failed = 0;
function check(name, cond, detail) {
    if (cond) { console.log(`  [OK] ${name}`); }
    else { failed++; console.log(`  [FAIL] ${name}  ${detail || ""}`); }
}

setTimeout(() => {
    console.log("== 启动与示例数据 ==");
    const stats = doc.getElementById("stats").textContent;
    check("示例数据已载入统计", stats.includes("n = 30"), stats.slice(0, 60));
    check("成果表 13 行", doc.querySelectorAll("#tbody tr").length === 13);
    check("无脚本错误", errors.length === 0, errors.join(";"));

    console.log("== 手动输入并计算 ==");
    const ta = doc.getElementById("data");
    ta.value = "1990 100\n1991 120\n1992 90\n1993 150\n1994 110\n1995 95\n1996 130\n1997 105\n1998 85\n1999 140";
    doc.getElementById("btnCalc").click();
    check("统计更新 n=10", doc.getElementById("stats").textContent.includes("n = 10"));
    check("成果表仍 13 行", doc.querySelectorAll("#tbody tr").length === 13);
    check("来源信息", doc.getElementById("srcInfo").textContent.includes("10 个数据"));

    console.log("== 参数调节（滑杆联动） ==");
    const cvNum = doc.getElementById("cvNum"), cvRng = doc.getElementById("cvRng");
    cvRng.value = "0.5";
    cvRng.dispatchEvent(new dom.window.Event("input", { bubbles: true }));
    check("Cv 滑杆→数字框", cvNum.value === "0.5", `got ${cvNum.value}`);
    cvNum.value = "0.42";
    cvNum.dispatchEvent(new dom.window.Event("input", { bubbles: true }));
    check("Cv 数字框→滑杆", cvRng.value === "0.42", `got ${cvRng.value}`);
    const meanNum = doc.getElementById("meanNum");
    check("均值固定时输入禁用", meanNum.disabled === true);
    doc.getElementById("chkFix").click();  // 取消固定
    check("取消固定后均值可调", meanNum.disabled === false);

    console.log("== 矩法估计与优化适线 ==");
    doc.getElementById("btnMoment").click();
    doc.getElementById("btnFit").click();
    setTimeout(() => {
        const fi = doc.getElementById("fitInfo").textContent;
        check("适线结果给出 Cv/Cs", /Cv=\d/.test(fi), fi.slice(0, 60));
        check("状态栏提示适线完成", doc.getElementById("status").textContent.includes("优化适线完成"));
        check("拟合检验标签", doc.getElementById("errInfo").textContent.includes("拟合检验"));
        check("全程无脚本错误", errors.length === 0, errors.join(";"));

        console.log("== 异常路径 ==");
        doc.getElementById("btnClear").click();
        check("清空后无成果表", doc.querySelectorAll("#tbody tr").length === 0);
        doc.getElementById("btnCalc").click();   // 空数据 → alert
        ta.value = "1\n2";                        // 数据不足
        doc.getElementById("btnCalc").click();   // alert
        ta.value = "abc\ndef";                    // 非数值
        doc.getElementById("btnCalc").click();   // alert
        check("异常路径无脚本错误", errors.length === 0, errors.join(";"));

        console.log("== 导出（桩） ==");
        doc.getElementById("btnPng").click();
        doc.getElementById("btnCsv").click();
        check("导出无异常", errors.length === 0, errors.join(";"));

        console.log();
        if (failed || errors.length) {
            console.log(`失败 ${failed} 项，脚本错误 ${errors.length} 条`);
            process.exit(1);
        }
        console.log("网页版 UI 测试全部通过 ✔");
        process.exit(0);
    }, 150);
}, 400);
