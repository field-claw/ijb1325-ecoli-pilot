# -*- coding: utf-8 -*-
"""ecoli_darkmatter_eval —— E. coli ecGEM「注释 / kcat 双缺口」基准（eval_hub 第 6 个条目）。

来源（对比学习：钟博子韬 Bozitao Zhong, Flash Talk "Virtual E. coli", Westlake AIVC Week 2026）：
    其 Venus 生态（VenusX / VenusRXN / MAPred / VenusFactory）从蛋白语言模型侧补齐酶功能
    「暗物质」。要判断这条路对 VitaMind 的 GECKO 酶约束栈有多大价值，必须先量化：
        E. coli ecGEM 里到底有多少反应是「有基因但无 EC / 无 kcat」？
    本模块把该量化固化为**可复现的正式 benchmark**，与 Yeast8 生长基准并列。

指标（全部为**可数指标**，不依赖求解器、不需 GPU，确定性可复现）：
    - kcat_coverage / ec_coverage（双口径：原始 ID 计数 vs base 去重）
    - dark_matter_full（完全暗物质：有 GPR · 无 EC · 无 kcat）
    - ec_missing_only（有 kcat 无 EC：注释可零成本回填）
    - inherit_potential（变体继承 base EC 可修复的反应数）
    - true_targets_in_top50（归并 base 后仍无 EC 的高酶池占用靶标数）

诚实边界（已写进 SPEC，登记即随记录落盘）：
    - 完全暗物质 706–2923 是**口径差异**（同工酶变体合并与否），非计算错误，须注明口径。
    - **补 EC ≠ 补 kcat**：GECKO 要 kcat_MW，Venus 生态无公开 kcat 预测能力。
    - base 级覆盖率已隐含继承 → 继承在 base 级增益为 0，收益只在反应级；
      不得把「口径选择」宣传为「注释修复增益」。
    - 本基准**不预测、不评测模型精度**，只量化缺口——是「缺口审计」而非「性能基准」。

用法：
    python vitamind_core/ecoli_darkmatter_eval.py            # 算指标 + 登记 eval_hub
    from vitamind_core.ecoli_darkmatter_eval import run, ECOLI_DARKMATTER_SPEC
    res = run(register=False)
"""
from __future__ import annotations

import json
import os
import sys
from typing import Dict, Optional

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GAP_EVIDENCE = os.path.join(ROOT, "data", "ecoli_darkmatter_gap.json")
TARGET_EVIDENCE = os.path.join(ROOT, "data", "ecoli_ec_annotation_targets.json")
EVAL_HUB_DIR = os.path.join(ROOT, "data", "eval_hub")

BENCHMARK_ID = "ecoli_darkmatter_gap"

# —— 基准规格（制度化） —— #
ECOLI_DARKMATTER_SPEC: Dict[str, object] = {
    "name": "E. coli ecGEM annotation / kcat dual-gap audit",
    "source": ("VitaMind 自建（对比学习触发：钟博子韬 Bozitao Zhong, Flash Talk "
               "'Virtual E. coli', Westlake AIVC Week 2026）"),
    "model_under_audit": "eciML1515（真实 E. coli GECKO ecGEM，5883 rxn / 1516 gene / 2399 kcat_MW）",
    "task": ("审计：给定 ecGEM，统计「有 GPR 但无 EC 号 / 无 kcat」的反应规模，"
             "作为酶功能注释路线（PLM/结构/反应多模态）的靶标规模"),
    "benchmark_type": "缺口审计（可数指标），**非**预测精度基准",
    "metrics": [
        "kcat_coverage_raw", "kcat_coverage_dedup",
        "ec_coverage_raw", "ec_coverage_dedup",
        "dark_matter_full_raw", "dark_matter_full_dedup",
        "ec_missing_only_dedup", "inherit_potential", "true_targets_in_top50",
    ],
    "why_two_calibers": ("eciML1515 是 ecModels 格式，66.4%（3895/5870）反应 ID 带 "
                         "_numN/_reverse 同工酶拆分后缀；按 ID 直数会重复计同一生化反应，"
                         "故必须同时给「原始 ID 计数」与「base 去重」两个口径"),
    "honest_boundaries": [
        "完全暗物质 = 706–2923 条（占酶促 31%–54%），是口径定义差异而非计算错误；引用须注明口径",
        "补 EC ≠ 补 kcat：GECKO 需要 kcat_MW（s⁻¹·MW⁻¹），Venus 生态无公开 kcat 预测能力",
        "base 级覆盖率（48.32%）已隐含继承（任一变体有 EC 即算覆盖）→ 变体继承在 base 级增益为 0",
        "继承的真实收益只在反应级（+296 条，27.44%→32.93%），不得把口径选择宣传为注释修复增益",
        "本基准不评测任何模型的预测精度，只量化缺口规模；不构成对钟博子韬报告内容的描述",
        "kcat_MW 为 ecModels 原始文献自带参数，非本机实测，未逐条复核",
    ],
    "reproduce": ("python vitamind_core/ecoli_darkmatter_gap.py && "
                  "python vitamind_core/ecoli_darkmatter_ec_targets.py && "
                  "python vitamind_core/ecoli_darkmatter_eval.py"),
}


def _load(path: str) -> Optional[dict]:
    if not os.path.exists(path):
        return None
    return json.load(open(path, encoding="utf-8"))


def run(register: bool = True, hub_root: str = EVAL_HUB_DIR) -> Dict[str, object]:
    """读预计算证据，算指标，登记 eval_hub。"""
    gap = _load(GAP_EVIDENCE)
    tgt = _load(TARGET_EVIDENCE)
    if gap is None:
        return {"status": "missing_evidence",
                "missing": GAP_EVIDENCE,
                "fix": "先跑 python vitamind_core/ecoli_darkmatter_gap.py"}

    c = gap.get("counts", {})
    cov = gap.get("coverage", {})
    dedup = gap.get("dedup", {})
    dark = gap.get("dark_matter_full", {})
    ecgap = gap.get("ec_gap_only", {})
    orph = gap.get("orphan_reactions", {})
    top50 = gap.get("top50_bottleneck", {})

    full = (tgt or {}).get("full_gem_inherit_potential", {})
    metrics: Dict[str, object] = {
        "n_reactions_backbone": c.get("n_total"),
        "n_enzymatic_raw": c.get("n_with_gpr"),
        "n_enzymatic_dedup": dedup.get("n_enzymatic"),
        "kcat_coverage_raw": cov.get("kcat_coverage_of_enzymatic"),
        "kcat_coverage_dedup": dedup.get("kcat_coverage_of_enzymatic"),
        "ec_coverage_raw": cov.get("ec_coverage_of_enzymatic"),
        "ec_coverage_dedup": dedup.get("ec_coverage_of_enzymatic"),
        "dark_matter_full_raw": dark.get("n"),
        "dark_matter_full_dedup": dedup.get("dark_full"),
        "ec_missing_only_dedup": dedup.get("ec_missing_only"),
        "orphan_reactions": orph.get("n"),
        "top50_no_ec_count": top50.get("n_without_ec"),
        "top50_no_ec_pool_pct": top50.get("pool_pct_without_ec"),
        "inherit_potential": full.get("n_inheritable_enzymatic"),
        "ec_coverage_raw_after_inherit": full.get("reaction_level_ec_coverage_after_inherit"),
        "true_targets_in_top50": (tgt or {}).get("n_true_targets"),
        "true_targets_pool_pct": (tgt or {}).get("pool_pct_true_targets"),
    }

    out: Dict[str, object] = {"status": "ok", "benchmark": BENCHMARK_ID, "metrics": metrics}

    if register:
        from vitamind_core.eval_hub import EvalHub, record

        hub = EvalHub(root=hub_root)
        recs = [record(
            id=f"{BENCHMARK_ID}:spec",
            question=("E. coli ecGEM 注释/kcat 双缺口基准：有 GPR 但无 EC / 无 kcat 的反应规模？"
                      "（酶功能注释路线的靶标规模审计）"),
            data=("source: eciML1515 真实 GECKO ecGEM | trigger: 钟博子韬 Flash Talk "
                  "'Virtual E. coli' (Westlake AIVC Week 2026)"),
            answer=json.dumps(ECOLI_DARKMATTER_SPEC, ensure_ascii=False),
            has_file=False,
            metadata={"type": "gap_audit_spec", "role": "third_party_benchmark",
                      "recommended_by": "钟博子韬 (Westlake AIVC Week 2026, Flash Talk)"},
        )]

        # 指标记录
        recs.append(record(
            id=f"{BENCHMARK_ID}:coverage",
            question="eciML1515 的 kcat 与 EC 覆盖率是多少（双口径）？",
            data="eciML1515 | 原始 ID 计数 vs base 去重（合并 _numN/_reverse 同工酶变体）",
            answer=round(float(metrics["kcat_coverage_dedup"] or 0), 4),
            has_file=False,
            metadata={"kcat_coverage_raw": metrics["kcat_coverage_raw"],
                      "kcat_coverage_dedup": metrics["kcat_coverage_dedup"],
                      "ec_coverage_raw": metrics["ec_coverage_raw"],
                      "ec_coverage_dedup": metrics["ec_coverage_dedup"],
                      "metric_source": "ecoli_darkmatter_gap.json"},
        ))
        recs.append(record(
            id=f"{BENCHMARK_ID}:dark_matter",
            question="完全暗物质反应（有 GPR · 无 EC · 无 kcat）有多少条？",
            data="eciML1515 | 双口径区间",
            answer=float(metrics["dark_matter_full_dedup"] or 0),
            has_file=False,
            metadata={"dark_matter_raw": metrics["dark_matter_full_raw"],
                      "dark_matter_dedup": metrics["dark_matter_full_dedup"],
                      "caliber_range": "706–2923（占酶促 31%–54%）",
                      "metric_source": "ecoli_darkmatter_gap.json"},
        ))
        recs.append(record(
            id=f"{BENCHMARK_ID}:inherit_potential",
            question="变体继承 base EC 可修复多少条反应（零 AI 成本）？",
            data="eciML1515 | 自身无 EC 但 base 家族有 EC",
            answer=float(metrics["inherit_potential"] or 0),
            has_file=False,
            metadata={"inherit_potential": metrics["inherit_potential"],
                      "ec_coverage_raw_before": metrics["ec_coverage_raw"],
                      "ec_coverage_raw_after": metrics["ec_coverage_raw_after_inherit"],
                      "base_level_gain": 0,
                      "caliber_warning": (full.get("caliber_warning")
                                          or "base 级口径已隐含继承，继承增益为 0")},
        ))
        recs.append(record(
            id=f"{BENCHMARK_ID}:top50_targets",
            question="Top50 瓶颈酶中归并 base 后仍无 EC 的真靶标有几条、占酶池多少？",
            data="eciML1515 Top50 瓶颈酶（占酶池 62.81%）| 归并 base 再判定",
            answer=float(metrics["true_targets_in_top50"] or 0),
            has_file=False,
            metadata={"n_input_no_ec": metrics["top50_no_ec_count"],
                      "pool_pct_no_ec_listed": metrics["top50_no_ec_pool_pct"],
                      "n_true_targets": metrics["true_targets_in_top50"],
                      "true_targets_pool_pct": metrics["true_targets_pool_pct"],
                      "targets": [t["reaction"] for t in
                                  (tgt or {}).get("true_targets", [])],
                      "metric_source": "ecoli_ec_annotation_targets.json"},
        ))

        path = hub.register(
            BENCHMARK_ID, recs,
            source_note=("E. coli ecGEM 注释/kcat 双缺口审计（对比学习触发：钟博子韬 Virtual E. coli）；"
                         "可数指标、双口径（原始 ID vs base 去重）；"
                         "补 EC ≠ 补 kcat，不评测模型精度"))
        out["registered_path"] = path
        out["n_records"] = len(recs)
    return out


if __name__ == "__main__":
    res = run(register=True)
    print("=== E. coli ecGEM 注释 / kcat 双缺口基准 ===")
    print(f"status: {res.get('status')}")
    if res.get("status") == "ok":
        m = res["metrics"]
        print(f"  主干反应 {m['n_reactions_backbone']} ｜ 酶促 原始 {m['n_enzymatic_raw']} / 去重 {m['n_enzymatic_dedup']}")
        print(f"  kcat 覆盖 : 原始 {m['kcat_coverage_raw']:.2%} ｜ 去重 {m['kcat_coverage_dedup']:.2%}")
        print(f"  EC   覆盖 : 原始 {m['ec_coverage_raw']:.2%} ｜ 去重 {m['ec_coverage_dedup']:.2%}")
        print(f"  完全暗物质: 原始 {m['dark_matter_full_raw']} ｜ 去重 {m['dark_matter_full_dedup']}"
              f"   (区间 706–2923)")
        print(f"  有kcat无EC（注释可回填）: {m['ec_missing_only_dedup']}")
        print(f"  变体继承可修复: {m['inherit_potential']} 条 "
              f"(反应级 EC {m['ec_coverage_raw']:.2%} → {m['ec_coverage_raw_after_inherit']:.2%})")
        print(f"  Top50 无EC {m['top50_no_ec_count']} 条（酶池 {m['top50_no_ec_pool_pct']}%）→ "
              f"归并base后真靶标 {m['true_targets_in_top50']} 条（酶池 {m['true_targets_pool_pct']}%）")
        if "registered_path" in res:
            print(f"\n已登记 eval_hub: {res['registered_path']} ({res['n_records']} 条记录)")
