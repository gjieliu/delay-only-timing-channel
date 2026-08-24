# -*- coding: utf-8 -*-
"""端到端复现：多分辨率主结果、走廊损失下端、相关宿主表、清零对照表、单字母窗口、GMI 费用与全部论文图表。

缺省快速档：随机数值脚本以 --quick 运行（分钟级冒烟，写 *_quick.json，不覆盖正式 JSON）。
确定性的 multires_bound 与 make_figs 仍重建同名正式基线和图件。
`python3 run_all.py --full` 跑正式档（覆盖正式 JSON；相关宿主表约 30–60 分钟）。
校验：随后运行 `python3 check_consistency.py`（只认正式 JSON）。
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    full = "--full" in sys.argv
    import lambda_delta_minus
    import make_figs
    import multires_bound

    print("[1/9] multires_bound")
    multires_bound.main()
    print("\n[2/9] certify_delta_plus")
    subprocess.check_call([sys.executable, os.path.join(HERE, "certify_delta_plus.py")])
    print("\n[3/9] lambda_delta_minus")
    lambda_delta_minus.main(quick=not full)
    print("\n[4/9] make_figs")
    make_figs.main()
    flag = [] if full else ["--quick"]
    for i, script in enumerate(
            ("correlation_full.py", "delta_minus_width_bounds.py", "gmi_demo.py",
             "gmi_finite.py", "run_ding.py"), start=5):
        print("\n[%d/9] %s %s" % (i, script, " ".join(flag)))
        subprocess.check_call([sys.executable, os.path.join(HERE, script)] + flag)
    print("\nrun_all done" + ("" if full else "（快速档；正式结果请用 --full 重跑）"))
