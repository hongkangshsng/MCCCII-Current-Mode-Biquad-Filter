# MCCCII 主動式電流模式二階濾波器

Current-mode biquad filters using a single multi-output current-controlled conveyor.

洪鏮展｜輔仁大學電機工程學系工程作品集。由原始專題報告出發，重建 LPF、HPF、BPF 電路方程、修正可確認的公式錯誤，並以 Python 建立可重現的理想頻率響應與元件調參實驗。

**Evidence boundary / 證據界線：** 新增曲線是 Rx=0 理想模型的數值模擬；原報告圖表是歷史來源，尚無可供重新分析的原始實測資料。本專案不宣稱已完成硬體重新量測或電晶體層模擬。

## 開始使用 / Quick start

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python run_experiments.py
```

執行目錄為 repository 根目錄。輸出包含 21 組頻率響應 CSV、9 張增益／相位 SVG 曲線、掃描表與環境版本。Run from the repository root; no external dataset is required.

## 技術文件 / Technical documentation

- [電路拓撲、完整數學推導與原報告勘誤 / Derivation and audit](docs/derivation.md)
- [參數配置、調參方法、全部曲線與驗證 / Experiments and validation](docs/experiments.md)
- [掃描結果表 / Generated results](results/summary.md)
- [中英文專題摘要 / Bilingual project summary](docs/project-summary.md)
- [參考文獻 / References](references/references.md)

| Baseline | f0 (Hz) | Q |
|---|---:|---:|
| LPF | 159.154943 | 0.500000 |
| HPF | 159.154943 | 0.500000 |
| BPF | 91.888149 | 0.346410 |

![LPF Q tuning](results/LPF_q.svg)

## 工程能力 / Engineering scope

本次補強涵蓋節點分析、二階系統係數比較、自然頻率與 Q、R/C/電導調參、數值輸出及獨立求解器測試。程式 `src/filters.py` 提供三種模型；`run_experiments.py` 自動輸出曲線與資料；`tests/test_filters.py` 驗證節點解、解析式、穩定性與頻寬。CI 設定在 `.github/workflows/verify.yml`。

This project reconstructs three admittance configurations from the author's report, derives transfer functions using KCL, corrects documented inconsistencies, and supplies reproducible parameter sweeps with independent nodal-solver checks. Seven local test methods passed. Finite Rx, frequency-dependent tracking, noise, physical implementation and measured-data comparison remain future work.

All new plots and tables represent ideal transfer-function calculations, not new physical measurements. The natural frequency f0 is distinct from the LPF/HPF −3 dB cutoff. See the audit for topology-label ambiguities and the corrected BPF damping coefficient.
