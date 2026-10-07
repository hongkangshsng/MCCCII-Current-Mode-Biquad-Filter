# 主動式電流模式二階濾波器設計與實現
### Current-Mode Biquad Filter using MCCCII

> 輔仁大學電機工程相關專題作品集｜類比電路設計、電流模式濾波器、數學分析與實驗驗證

本專題以**單一多輸出電流控制傳輸器（Multi-output Current Controlled Conveyor, MCCCII）**搭配被動元件，研究電流模式二階濾波器的設計與實現。透過不同導納配置，可由相同的主動元件架構實現**低通（LPF）、高通（HPF）與帶通（BPF）**等濾波功能。

**核心技術：** `Analog Circuit` · `MCCCII` · `Current-Mode Circuit` · `LPF` · `HPF` · `BPF` · `Transfer Function` · `Sensitivity Analysis` · `MATLAB`

---

## 1. 專題背景與研究目的

電流模式（Current-Mode）電路具有適合類比訊號處理的特性。本專題以 MCCCII 作為主要主動元件，探討如何利用單一 MCCCII 與 R、C 等被動元件建立二階濾波器，並從電路方程式、傳輸函數、中心頻率、品質因數、靈敏度以及實驗結果等面向進行分析。

研究重點包括：

- 使用單一 MCCCII 實現多種二階濾波功能
- 分析 Low-Pass、High-Pass、Band-Pass Filter
- 推導 Transfer Function H(s)
- 分析中心頻率 ω0 與品質因數 Q
- 探討 MCCCII 非理想特性與 Current Tracking Error
- 進行主動／被動 Sensitivity Analysis
- 比較理論分析與實驗／模擬結果

---

## 2. MCCCII 工作原理

理想 MCCCII 的基本關係可表示為：

```text
Vx = Vy + IxRx
Iy = 0
Iz+ = Ix
Iz- = -Ix
```

透過 X、Y、Z+、Z− 等端點的電壓／電流關係，以及外部導納元件的不同配置，可形成不同的電流模式濾波器。

---

## 3. 二階濾波器電路架構

本專題研究的單一 MCCCII 架構可透過不同的被動導納 Y1～Y4，合成不同的二階濾波器功能。主要研究三種：

| 濾波器 | 英文 | 功能 |
|---|---|---|
| 低通濾波器 | Low-Pass Filter (LPF) | 保留低頻訊號並衰減高頻訊號 |
| 高通濾波器 | High-Pass Filter (HPF) | 衰減低頻並保留高頻訊號 |
| 帶通濾波器 | Band-Pass Filter (BPF) | 允許指定頻帶內的訊號通過 |

---

## 4. 低通濾波器（LPF）

其中一組低通濾波器導納配置為：

```text
Y1 = sC1
Y2 = G2
Y3 = sC3
Y4 = G4
```

其傳輸函數可表示為：

```text
Io/Iin = G2G4 /
         [s²C1C3 + s(C1G2 + C3G2) + G2G4]
```

中心頻率為：

```text
ω0 = √(G2G4 / C1C3)
```

藉由電容與電導參數，可進一步分析濾波器的中心頻率與品質因數。

---

## 5. 高通與帶通濾波器（HPF / BPF）

透過改變 MCCCII 周圍的被動元件與導納配置，同一類型架構亦可實現 High-Pass 與 Band-Pass 二階濾波器。本專題分別分析其傳輸函數、中心頻率、品質因數，以及 Gain / Phase Response。

此設計概念的重點在於：**以較少的主動元件實現多種濾波功能，使電路架構更加簡化與彈性化。**

---

## 6. 非理想模型與靈敏度分析

考慮 MCCCII 非理想特性後，可引入追蹤係數 α、β：

```text
Vx = αVy + IxRx
Iz+ = βIx
```

接著分析電路元件與主動元件非理想效應對下列參數的影響：

- 中心頻率 ω0
- 品質因數 Q
- Passive Sensitivity
- Active Sensitivity
- Current Tracking Error

專題分析顯示，所研究的濾波器具有較低的被動靈敏度；中心頻率與品質因數對 MCCCII 電流追蹤誤差亦具有一定穩健性。

---

## 7. 實驗與模擬驗證

專題針對不同濾波器進行 Gain Response 與 Phase Response 分析，包括：

- LPF Gain / Phase Response
- HPF Gain / Phase Response
- BPF Gain / Phase Response

透過理論結果與實驗／模擬數據進行交叉比較，結果僅呈現小幅差異，可用以驗證所提出電路架構與理論分析的可行性。

---

## 8. 本專題展現的工程能力

透過本專題完成並練習以下能力：

- **類比電路分析（Analog Circuit Analysis）**
- **主動式濾波器設計（Active Filter Design）**
- **電流模式電路（Current-Mode Circuit）**
- **二階系統與傳輸函數推導**
- **中心頻率 ω0 / 品質因數 Q 分析**
- **頻率響應、Gain、Phase Analysis**
- **Sensitivity Analysis**
- **MATLAB 訊號／頻率響應分析**
- **工程論文閱讀與文獻整理**
- **理論與實驗結果比較**

---

## 9. 專題成果

本專題研究顯示，使用**單一 MCCCII** 即可建立電流模式二階低通、高通與帶通濾波器。透過數學推導、靈敏度分析與實驗結果驗證，進一步確認此類架構在類比訊號處理與主動式濾波器設計上的可行性。

---

## 10. Repository Structure

```text
MCCCII-Current-Mode-Biquad-Filter/
├── README.md                 # 中文作品集首頁
├── docs/
│   └── project-summary.md    # 專題摘要
└── references/
    └── references.md        # 主要參考文獻
```

---

## 11. 參考文獻

本專題參考 Current Conveyor、Current-Mode Filter、CFCCII 與 MCCCII 等相關研究。完整整理請參閱：

➡️ [References](references/references.md)

---

## 作者 / Author

**洪鏮展**  
輔仁大學 電機工程學系  
Electrical Engineering, Fu Jen Catholic University