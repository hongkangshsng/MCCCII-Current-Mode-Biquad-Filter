# 專題摘要 / Project summary

以單一 MCCCII 的導納網路為基礎，研究電流模式二階低通、高通及帶通濾波器。由原始報告電路圖重建節點方程，釐清導納標籤與帶通公式的不一致，推導 H(s)、自然頻率及品質因數。以 Python 產生 21 組理想頻率響應、9 張增益／相位圖及元件掃描表，並以獨立節點求解、SciPy 及解析性質驗證程式。本次七項本機測試通過；新增結果屬理論模型數值模擬，無新增硬體量測。原報告中的實驗圖未作為本次程式輸出或量化誤差證據。

Reconstructed single-MCCCII LPF, HPF and BPF admittance networks from the author's report; derived transfer functions and audited formula and labeling errors. Implemented 21 reproducible ideal-model cases with gain/phase plots and component sweep tables. Validated the implementation using independent nodal solutions, SciPy frequency responses, pole stability and analytical bandwidth checks. Seven local test methods passed. New outputs are ideal-model numerical simulations; no new physical measurements are claimed.

[推導 / Derivation](derivation.md) · [實驗 / Experiments](experiments.md)
