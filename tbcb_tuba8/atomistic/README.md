# TBCB–TUBA8 原子级中间态物理计算（2026-11-04 起）

执行 `protocol/LOCAL_EXECUTION_PROTOCOL.md` 及其补充要求中的全部原子级计算：折叠中间态（REST2、伞形采样 PMF、310 K 无偏 MD、MSM、高温展开）、残基替换的折叠自由能 ΔΔG_fold（炼金 FEP）、TBCB/TBCE 对 α 尾巴与中间态的结合（相对/绝对结合自由能、Kd、τRAMD 解离、PMF 推导的 kon/koff、MM/GBSA–PBSA 排序、集合对接）、TBCE–α 与 α/β 异二聚体界面。**不含**任何动力学网络模型拟合或 Purkinje 反事实（按要求只做原子级）。

## 体系与力场（已完成）
- 55 个起始结构，全部保留 α N 位点 GTP·Mg²⁺：14 个 α 单体（TUBA1A/1B/4A/8、I114S、S114I、双向 loop 交换、A40K、K40A、G43I、I43G、loop+114 双突变）、AF3 TBCB–α / TBCE–α（5 个 seed）、9M1L TBCE–α、4 个 α/β 异二聚体（9M1M）、10 个 CAP-Gly–尾巴体系（TBCB WT/K188A × Y/F 尾巴、TBCE 槽）。
- 力场：ff14SB（ISOLDE 转换，与 OpenMM 自带 ff14SB 逐项比对：键/角/非键完全一致，扭转差 0.09/4970 kcal/mol）+ Meagher 2003 GTP + Li–Merz IOD Mg²⁺（电荷 +2；ISOLDE 原文件 Mg 电荷写成 1.0，已改）+ TIP3P/JC 离子。Mg²⁺ 内壳层从 1.8 Å 晶体 4I4T 嫁接（β/γ 磷酸氧 + 4 个水），最小化后保持八面体。
- 能量最小化与残余力：`results/minimization_forces.csv`（每个结构：初始/限制性最小化后/自由最小化后的势能、最大受力及其原子、RMS 受力；复合物另有真空原样受力、界面接触数与界面 Cα RMSD）。
- ISOLDE 的两个参数文件因许可证不在仓库中：`ff/fetch_isolde_ff.sh` 按固定 commit 下载并重建 `gtp_meagher.xml`（md5 见下）。

## 自动流程（已配置并在集群上排队）
- Sherlock：8 个分析作业（fep, rest2, md310, md450, pmf, usphi, tail, mmgbsa）各自以 `afterany` 依赖其计算作业；`report_tbcb` 依赖全部计算与分析作业，运行 `a21_compile.py` + `a24_report.py` 生成 `master.csv/json` 与中文报告 `report_zh.md`，并以 Slurm 邮件（END,FAIL）通知。作业号见 `results/manifests/`。
- Marlowe：wave-2 的 4 个分析作业 + `report_tbcb_marlowe`（同上，邮件 END,FAIL）。
- 两个集群之间无法由作业互相取数，**跨集群合并、推送到本仓库需要我的会话执行**（集群上不存放 GitHub 凭据）。
- 结果目录：Sherlock `$SCRATCH/tbcb_atomistic/{runs,results}`（90 天清除）与 `$HOME/tbcb_atomistic/results_summary`（不清除，小文件）；Marlowe `/scratch/m000251/zongminl/tbcb_mech/atomistic/{runs,results}`。

## 物理量的类别
每个数都标注：实算物理量 / 由实算量经物理模型推导 / 排序分数（MM/GBSA、HADDOCK，非自由能，不换算 Kd）/ 未计算（附原因）。未通过协议规定验证门（λ 窗重叠 < 0.03、半程漂移超出重复离散、MSM 未过 CK 检验、PMF 未到体相）的数一律不报。

## 明确不计算的
全长 TBCB–α 的绝对结合自由能（~5000 Å² 界面）；TBCB→TBCE 转移的反应路径/能垒/速率（不存在 TBCB–TBCE–α 原子结构）；K40A/A40K 的炼金 ΔΔG（改变净电荷，PME 有限尺寸误差大于信号，改用 MD/REST2 布居）；G43I/I43G 的炼金 ΔΔG（甘氨酸主链参数不同，改用 φ43 伞形采样 PMF）。
