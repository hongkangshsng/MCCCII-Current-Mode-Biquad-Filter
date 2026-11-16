# MCCCII 電流模式二階濾波器：電路架構與靈敏度分析筆記

* **專題題目**：Realization of Current-Mode Biquad using Single Multi-output Current Controlled Conveyor (MCCCII)
* **所屬單位**：天主教輔仁大學 電機工程學系（指導老師：鄞永昌）
* **線上互動實驗台**：[MCCCII Interactive Lab (Streamlit)](https://www.google.com/search?q=https://mcccii-current-mode-biquad-filter-clf2bbkf4wffc2dcqvhmvq.streamlit.app/)

---

## 1. MCCCII 主動元件埠口特性

多輸出電流控制輸送器（Multi-output Current Controlled Conveyor, MCCCII）為本專題的核心主動元件，其理想與非理想埠口關係如下：

| 模型狀態 | 電壓與電流埠口關係 | 參數說明 |
| :--- | :--- | :--- |
| **理想模型** | $V_x = V_y + I_x R_x$，$I_y = 0$，$I_{Z+} = I_x$，$I_{Z-} = -I_x$ | $Y$ 端高輸入阻抗（$I_y = 0$）；理想分析常假設 $R_x = 0$ |
| **非理想模型** | $V_x = \alpha V_y + I_x R_x$，$I_{Z+} = \beta I_x$ | $\alpha = 1 - \varepsilon_1$（電壓追蹤誤差，$\lvert\varepsilon_1\rvert \ll 1$）；$\beta = 1 - \varepsilon_2$（電流追蹤誤差，$\lvert\varepsilon_2\rvert \ll 1$） |

---

## 2. 通用單一 MCCCII 架構與三種濾波器實現

透過單顆 MCCCII 搭配四組導納網路（$Y_1, Y_2, Y_3, Y_4$），可推導出通用電流模式二階轉移函數：

$$
\frac{I_o}{I_{in}} = \frac{Y_2 Y_4}{Y_1 Y_2 + Y_1 Y_3 + Y_2 Y_3 + Y_2 Y_4}
$$

依據不同的被動元件（$G_i, C_i$）配置，可由同一主動架構實現低通（LPF）、高通（HPF）與帶通（BPF）濾波器（對應專案 `src/filters.py` 實作）：

| 濾波器型態 | 導納配置 ($Y_1, Y_2, Y_3, Y_4$) | 轉移函數 $\frac{I_o}{I_{in}}$ | 中心頻率 $\omega_0$ 與品質因數 $Q$ |
| :--- | :--- | :--- | :--- |
| **LPF (低通)** | $Y_1=sC_1,\; Y_2=G_2,\; Y_3=sC_3,\; Y_4=G_4$ | $\frac{G_2 G_4}{s^2 C_1 C_3 + s G_2(C_1 + C_3) + G_2 G_4}$ | $\omega_0 = \sqrt{\frac{G_2 G_4}{C_1 C_3}}$，$Q = \frac{1}{C_1+C_3}\sqrt{\frac{C_1 C_3 G_4}{G_2}}$ |
| **HPF (高通)** | $Y_1=G_1,\; Y_2=sC_2,\; Y_3=G_3,\; Y_4=sC_4$ | $\frac{s^2 C_2 C_4}{s^2 C_2 C_4 + s C_2(G_1 + G_3) + G_1 G_3}$ | $\omega_0 = \sqrt{\frac{G_1 G_3}{C_2 C_4}}$，$Q = \frac{1}{G_1+G_3}\sqrt{\frac{C_4 G_1 G_3}{C_2}}$ |
| **BPF (帶通)** | $Y_1=sC_1+G_1,\; Y_2=sC_2,\; Y_3=sC_3+G_3,\; Y_4=G_4$ | $\frac{s C_2 G_4}{\Delta s^2 + \Delta_1 s + G_1 G_3}$ | $\omega_0 = \sqrt{\frac{G_1 G_3}{\Delta}}$，$Q = \frac{\sqrt{G_1 G_3 \Delta}}{\Delta_1}$ |

其中 BPF 分母輔助係數為：
* $\Delta = C_1 C_2 + C_1 C_3 + C_2 C_3$
* $\Delta_1 = C_2 G_1 + C_1 G_3 + C_3 G_1 + C_2 G_3 + C_2 G_4$

### 架構擴充：單顆 vs. 雙顆 MCCCII 比較

| 架構配置 | 可實現濾波器型態 | 應用特點 |
| :--- | :--- | :--- |
| **單顆 MCCCII**（本 Repo 驗證核心） | LPF、HPF、BPF（3 種） | 電路結構精簡、主動元件數最少，適合低功耗設計 |
| **雙顆 MCCCII**（延伸架構） | LPF、HPF、BPF、BSF（帶阻）、APF（全通）（5 種） | 可進一步實現特定頻段抑制（BSF）與相位延遲補償（APF） |

---

## 3. 非理想追蹤誤差與電路靈敏度分析

考慮主動元件非理想電壓/電流追蹤誤差（$\alpha, \beta$）時，通用轉移函數修正為：

$$
\frac{I_o}{I_{in}} = \frac{Y_2 Y_4 + Y_3 Y_4 - \alpha Y_3 Y_4}{\beta(Y_1 Y_2 + Y_1 Y_3 + Y_2 Y_3 + Y_2 Y_4 + Y_3 Y_4 - \alpha Y_3 Y_4)}
$$

靈敏度（Sensitivity）定義為參數 $X_i$ 的相對變動對指標 $F$（如 $\omega_0$ 或 $Q$）造成的相對影響：

$$
S_{X_i}^{F} = \frac{X_i}{F} \frac{\partial F}{\partial X_i}
$$

### 3.1 低通濾波器（LPF）靈敏度
* $S_{G_2, G_4}^{\omega_0} = S_{G_4}^{Q} = \frac{1}{2}$，$S_{C_1, C_3}^{\omega_0} = S_{G_2}^{Q} = -\frac{1}{2}$
* $S_{C_1}^{Q} = \frac{1}{2} - \frac{C_1}{C_1 + C_3}$，$S_{C_3}^{Q} = \frac{1}{2} - \frac{C_3}{C_1 + C_3}$
* 主動追蹤靈敏度：$S_{\alpha}^{\omega_0} = S_{\beta}^{\omega_0} = S_{\beta}^{Q} = 0$，$S_{\alpha}^{Q} = \frac{\alpha C_3 G_4}{C_1 G_2 + C_3 G_2 + C_3 G_4 - \alpha C_3 G_4}$

### 3.2 高通濾波器（HPF）靈敏度
* $S_{G_1, G_3}^{\omega_0} = S_{C_4}^{Q} = \frac{1}{2}$，$S_{C_2, C_4}^{\omega_0} = S_{C_2}^{Q} = -\frac{1}{2}$
* $S_{G_1}^{Q} = \frac{1}{2} - \frac{G_1}{G_1 + G_3}$，$S_{G_3}^{Q} = \frac{1}{2} - \frac{G_3}{G_1 + G_3}$
* 主動追蹤靈敏度：$S_{\alpha}^{\omega_0} = S_{\beta}^{\omega_0} = S_{\beta}^{Q} = 0$，$S_{\alpha}^{Q} = \frac{\alpha C_4 G_3}{C_2 G_3 + C_2 G_1 + C_4 G_3 - \alpha C_4 G_3}$

### 3.3 帶通濾波器（BPF）靈敏度
* $S_{G_1, G_3}^{\omega_0} = \frac{1}{2}$，$S_{C_1}^{\omega_0} = -\frac{1}{2}\frac{C_1(C_2+C_3)}{\Delta}$，$S_{C_2}^{\omega_0} = -\frac{1}{2}\frac{C_2(C_1+C_3)}{\Delta}$，$S_{C_3}^{\omega_0} = -\frac{1}{2}\frac{C_3(C_1+C_2)}{\Delta}$
* $S_{G_4}^{Q} = -1$，$S_{G_1}^{Q} = 1 - \frac{G_1(C_2+C_3)}{\Delta_1}$，$S_{G_3}^{Q} = 1 - \frac{G_3(C_1+C_2)}{\Delta_1}$
* $S_{C_1}^{Q} = \frac{C_1(C_2+C_3)}{\Delta} - \frac{\Delta G_3}{\Delta_1}$，$S_{C_3}^{Q} = \frac{C_3(C_1+C_2)}{\Delta} - \frac{\Delta G_1}{\Delta_1}$，$S_{C_2}^{Q} = \frac{C_2(C_1+C_3)}{\Delta} - \frac{\Delta(G_1+G_3+G_4)}{\Delta_1}$
* 主動追蹤靈敏度：$S_{\beta}^{\omega_0} = S_{\beta}^{Q} = 0$，$S_{\alpha}^{\omega_0} = \frac{1}{2}$，$S_{\alpha}^{Q} = \frac{\alpha C_3 G_4}{\Delta_1 - \alpha C_3 G_4 + C_3 G_4}$

> **工程結論**：三種架構的被動靈敏度絕對值多數落於 $\le 0.5$（或 $\le 1$）的低靈敏度區間，且 $\omega_0$ 與 $Q$ 對電流追蹤誤差 $\beta$ 均呈零靈敏度（$S_\beta^{\omega_0} = S_\beta^Q = 0$）。實務元件公差造成的響應偏移可於本專案 **Streamlit 互動實驗台（Tab 2: Monte Carlo 公差模擬）** 中即時驗證。

---

## 4. 軟體驗證對照與證據邊界

* **已自動化驗證項目**：理想模型（$R_x = 0, \alpha = \beta = 1$）已於 `src/filters.py` 實作，並透過 `tests/test_filters.py` 的 7 項單元測試與 3×3 KCL 節點矩陣交叉比對（誤差 `< 2e-9`）。
* **未來擴充方向**：將有限寄生電阻 $R_x$、高頻追蹤誤差 $(\alpha, \beta)$ 與電晶體層級 Spice 模擬納入數值模型。
