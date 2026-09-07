# 外力驅動的非線性阻尼振盪

[ 中文 | [English](./README.en.md) ]

## 系統示意圖

<img src="./images/schematic_of_system_3.png" width="200">

### 物理模型設定
一個磁針（以 $\mu$ 表示）受到兩個分量組成的合成磁場驅動：
- $x$ 軸磁場：一個固定不變的靜磁場，記為 $B_1$（常數）。
- $y$ 軸磁場：一個隨時間週期變化的振盪磁場 $B_2 \cos(\omega t)$，角頻率為 $\omega$。

### 動力學與運動

當磁偶極矩 $\mu$ 與這些磁場交互作用時，會與 $x$ 軸形成夾角 $\theta$。

由於 $y$ 軸分量會隨 $\cos(\omega t)$ 振盪，合成磁場的方向也會隨時間週期性改變。

這種時變的驅動力會對磁針施加週期性力矩，使其繞平衡位置產生轉動振盪，如圖中彎曲的雙向箭頭所示。


## 無因次化後的運動方程式


$\frac{d^2\theta}{dT^2}=-\gamma\frac{d\theta}{dT}-b_1\sin\theta+b_2\cos\theta\cos (2\pi T)$


本專案採用 $\gamma = 6.0, b_1 = 36.0$。



---

## 特定初始條件下的分岔圖

<img src="./images/BifEuThAveIniTh0Btwo80(90)_0.1_120.png" >


---

## $b_2=95.0$ （2種 週期-1）

### $\theta-t$ 圖

<img src="./images/Btwo095.00(2period-1)/ThTBtwo95.00.png" height="300">

---

### 兩種振盪態的動畫

黑色箭頭為總磁場。
紅色箭頭為動畫中的週期- $1^+$ 振盪。
藍色箭頭為動畫中的週期- $1^-$ 振盪。

<img src="./images/Btwo095.00(2period-1)/Eu_Btwo95.00IniTh.25_Video.gif" height="300"> 

<img src="./images/Btwo095.00(2period-1)/Eu_Btwo95.00IniTh0_Video.gif" height="300">

---

### 相空間的二維投影

<img src="./images/Btwo095.00(2period-1)/OmeThBtwo95.00.png" height="300">

---

### 吸引域（Basins of attraction）

<img src="./images/Btwo095.00(2period-1)/OmeZeThZeBtwo95.00.png" height="300">

---

## $b_2=103.2$ （週期-3, 5共存）

### $\theta-t$ 圖

<img src="./images/Btwo103.20(period-3,5)/ThTBtwo103.20_RK4_4.png" height="300">

---

### 兩種振盪態的動畫

黑色箭頭為總磁場。
紅色箭頭為週期-3 振盪。
藍色箭頭為週期-5 振盪。

<img src="./images/Btwo103.20(period-3,5)/Eu_Btwo103.20IniTh.278_Video.gif" height="300"> 

<img src="./images/Btwo103.20(period-3,5)/Eu_Btwo103.20IniTh0_Video.gif" height="300">

---

### 相空間的二維投影

<img src="./images/Btwo103.20(period-3,5)/OmeThBtwo103.20_7.png" height="300">

---

### 吸引域（Basins of attraction）

<img src="./images/Btwo103.20(period-3,5)/OmeZeThZeBtwo103.20_RK4_2.png" height="300">

---
