# 電路核對與數學推導 / Topology audit and derivation

## 來源、假設與拓撲 / Source, assumptions and topology

依據洪鏮展原始中文報告（16 頁），PDF 頁碼 2–4 的 Fig.2–5 與式 (1)、(2)、(8)–(10)。頁碼為 PDF 頁序。原檔 SHA256：`0137f2bcf18123d7c51235491b6ecbf6b403564f1ffecd9cad65fd9a0f88acec`。原 PDF 未新增公開副本，避免擴散內含學號與第三方圖片；以下以文字重建連接關係，並列出歧義。

Source: the author's 16-page Chinese report, PDF pages 2–4. The original PDF is not republished. This document reconstructs the connectivity and explicitly records conflicts.

令輸入節點為 A；Y1：A 對地，Y2：Y 對地，Y3：A 到 Y，Y4：A 到 X。輸入電流流入 A，Ix 定義流入 X，Io 為 Z+ 輸出電流。Y 端不吸電流，Z− 未使用；輸出視為理想電流輸出，不含負載效應。

Input node A has shunt Y1; Y2 shunts terminal Y; Y3 connects A to Y; Y4 connects A to X. Ix enters X and Io leaves Z+. Terminal Y draws no current; Z− is unused.

**模型界線：** 原式 (1) 對應 `Vx=Vy`、`Io=Ix`，亦即忽略 Rx 的理想極限。原元件關係其實寫 `Vx=Vy+Ix Rx`；若 Rx 非零，以下 H(s) 不再完整。不能以這些曲線宣稱已驗證有限 Rx、電晶體頻寬、雜訊或硬體。

The implemented model uses Rx=0, unity tracking, and ideal unloaded current output. Finite Rx, transistor bandwidth, noise and hardware validation remain outside its scope.

## 節點推導 / Nodal derivation

由 Y 端 KCL、X 支路關係及 A 節點 KCL：

$$ (Y_2+Y_3)V_y=Y_3V_a,\qquad I_x=Y_4(V_a-V_y), $$
$$ I_{in}=Y_1V_a+Y_3(V_a-V_y)+I_x. $$

消去電壓且令 $I_o=I_x$ 得：

$$H(s)=\frac{I_o}{I_{in}}=\frac{Y_2Y_4}{Y_1Y_2+Y_1Y_3+Y_2Y_3+Y_2Y_4}.$$

These three equations are also solved independently in `nodal_response`; this checks topology against the polynomial implementation rather than comparing only two evaluators of the same polynomial.

## 三種配置 / Three admittance configurations

| Filter | Y1 (A-ground) | Y2 (Y-ground) | Y3 (A-Y) | Y4 (A-X) |
|---|---|---|---|---|
| LPF | sC1 | G2 | sC3 | G4 |
| HPF | G1 | sC2 | G3 | sC4 |
| BPF | G1+sC1 | sC2 | G3+sC3 | G4 |

$G=1/R$，單位 S；C 單位 F，s 單位 rad/s。下列寫 $H=N/(as^2+bs+c)$；$ω_0=\sqrt{c/a}$、$f_0=ω_0/(2\pi)$、$Q=\sqrt{ac}/b$。

### LPF

逐項代入分母為 $(sC_1)G_2+(sC_1)(sC_3)+G_2(sC_3)+G_2G_4$。

$$N=G_2G_4,\quad a=C_1C_3,\quad b=G_2(C_1+C_3),\quad c=G_2G_4,$$
$$ω_0=\sqrt{\frac{G_2G_4}{C_1C_3}},\quad Q=\frac{\sqrt{C_1C_3G_2G_4}}{G_2(C_1+C_3)}.$$

### HPF

逐項代入為 $G_1(sC_2)+G_1G_3+(sC_2)G_3+(sC_2)(sC_4)$。

$$N=s^2C_2C_4,\quad a=C_2C_4,\quad b=C_2(G_1+G_3),\quad c=G_1G_3,$$
$$ω_0=\sqrt{\frac{G_1G_3}{C_2C_4}},\quad Q=\frac{\sqrt{C_2C_4G_1G_3}}{C_2(G_1+G_3)}.$$

### BPF

逐項代入為 $(G_1+sC_1)sC_2+(G_1+sC_1)(G_3+sC_3)+sC_2(G_3+sC_3)+sC_2G_4$。收集次方後：

$$N=sC_2G_4,\quad a=k=C_1C_2+C_1C_3+C_2C_3,$$
$$b=C_2G_1+C_1G_3+C_3G_1+C_2G_3+C_2G_4,\quad c=G_1G_3,$$
$$ω_0=\sqrt{\frac{G_1G_3}{k}},\quad Q=\frac{\sqrt{kG_1G_3}}{b}.$$

BPF 峰值在 $ω_0$，$|H(jω_0)|=C_2G_4/b$，不保證單位增益。半功率頻寬為 $BW_ω=b/a$，$Q=ω_0/BW_ω$。

## 可確認修正 / Confirmed corrections and unresolved labels

| Location | 原文問題 / Issue | 處理 / Resolution |
|---|---|---|
| p2 Fig.2 | A-X 與 A-ground 同標 Y1 | A-X 採 Y4，與式 (1)、Fig.3/4 一致；literal Fig.2 is inconsistent |
| p3 Fig.4 下段 | HP 圖下面直接接 BP 敘述 | HP 依 Fig.4 另推；BP 段落依其導納配置分類 |
| p4 Fig.5 | 上支路 G2、Y 對地 C1；與 BP 文字 G4、C2 不一致 | 採 p3 BP 明確配置及 p12 參數：上支路 G4、Y-ground C2；此為有註記的重建 |
| p3 eq.(8) | BPF 阻尼少 C2G3 | 展開節點方程補上；分母共有五個一次項 |
| p3 eq.(10) | Q 分子未開根號 | 改為 sqrt(kG1G3)，使 Q 無因次 |
| p5 | 20log10(0.5) 寫成 −3 dB | 實為 −6.0206 dB；半功率幅度 1/sqrt(2) 為 −3.0103 dB |
| p8 | α=ε1≪1，卻又把 α 當電壓追蹤倍率 | 使用 Vx=αVy 時，理想 α=1；小誤差可定義 α=1−εv |
| p8 eq.(11) | β 在分母，與 Io=βIx 定義不合 | 依輸出定義 β 應乘在分子；未沿用原非理想式作模擬 |
| p9 LPF sensitivities | 電容 Q 靈敏度含 1−(1/2)Ci/(C1+C3) | 應為 1/2−Ci/(C1+C3) |
| p10 conclusion | 「無論追蹤誤差多大」皆不影響 | 僅常數輸出倍率 β 在此理想模型不移動極點；頻率相依追蹤或 Rx 不在此結論內 |

For finite Rx the branch relation is $I_x=Y_4(V_a-\alpha V_y)/(1+R_xY_4)$ and $I_o=\beta I_x$. This provides an extension path, not a tested hardware model. At Rx=0 a consistent expression is

$$H=\frac{\beta Y_4[Y_2+(1-\alpha)Y_3]}{Y_1Y_2+Y_1Y_3+Y_2Y_3+Y_4[Y_2+(1-\alpha)Y_3]}.$$

## 靈敏度 / Sensitivity

定義 $S_x^F=\partial\ln F/\partial\ln x$。由係數可直接求得：

$$S_x^{ω_0}=\tfrac12(S_x^c-S_x^a),\qquad S_x^Q=\tfrac12(S_x^a+S_x^c)-S_x^b.$$

LPF 的 G2/G4 對 $ω_0$ 均為 +1/2，C1/C3 為 −1/2。對 Q：G2=−1/2、G4=+1/2、$S_{C_i}^Q=1/2-C_i/(C_1+C_3)$。HPF：$S_{G_i}^Q=1/2-G_i/(G_1+G_3)$、C4=+1/2、C2=−1/2。BPF 可由上列 a,b,c 套通式；不能無條件宣稱所有非理想效應均低靈敏。
