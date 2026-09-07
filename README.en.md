# Nonlinear Driven Damped Oscillation

[ [中文](./README.md) | English ]

## Schematic of the system

<img src="./images/schematic_of_system_3.png" width="200">

### Physical Setup
A magnetic needle (represented by $\mu$ ) is driven by a combined magnetic field consisting of two components:  
- $x$ -axis magnetic field: A constant (static) magnetic field denoted as $B_1$ (const.).
- $y$ -axis magnetic field: An oscillating magnetic field given by $B_2 \cos(\omega t)$ , which varies periodically in time with an angular frequency $\omega$ .

### Dynamics and Motion

As the magnetic dipole $\mu$ interacts with these fields, it forms an angle $\theta$ with the $x$ -axis.  

Because the $y$-axis component oscillates ( $\cos(\omega t)$ ), the net magnetic field direction changes periodically.  

This time-dependent driving force exerts a periodic torque on the magnetic needle, causing it to undergo rotational oscillation around its equilibrium position, as indicated by the curved double-headed arrows.


## Rescaled equation of motion


$\frac{d^2\theta}{dT^2}=-\gamma\frac{d\theta}{dT}-b_1\sin\theta+b_2\cos\theta\cos (2\pi T)$


This project uses $\gamma = 6.0, b_1 = 36.0$ .



---

## Bifurcation diagram for specific initial conditions

<img src="./images/BifEuThAveIniTh0Btwo80(90)_0.1_120.png" >


---

## $b_2=95.0$ (2period-1)

### $\theta-t$ plot

<img src="./images/Btwo095.00(2period-1)/ThTBtwo95.00.png" height="300">

---

### Animation of two oscillation states

Black arrow is total magnetic field.  
Red arrow is period- $1^+$ oscillation in animation.  
Blue arrows is period- $1^-$ oscillation in animation.

<img src="./images/Btwo095.00(2period-1)/Eu_Btwo95.00IniTh.25_Video.gif" height="300"> 

<img src="./images/Btwo095.00(2period-1)/Eu_Btwo95.00IniTh0_Video.gif" height="300">

---

### 2D projection of the phase space

<img src="./images/Btwo095.00(2period-1)/OmeThBtwo95.00.png" height="300">

---

### Basins of attraction

<img src="./images/Btwo095.00(2period-1)/OmeZeThZeBtwo95.00.png" height="300">

---

## $b_2=103.2$ (period-3,5)

### $\theta-t$ plot

<img src="./images/Btwo103.20(period-3,5)/ThTBtwo103.20_RK4_4.png" height="300">

---

### Animation of two oscillation states

Black arrow is total magnetic field.  
Red arrow is period-3 oscillation in animation.  
Blue arrows is period-5 oscillation in animation.

<img src="./images/Btwo103.20(period-3,5)/Eu_Btwo103.20IniTh.278_Video.gif" height="300"> 

<img src="./images/Btwo103.20(period-3,5)/Eu_Btwo103.20IniTh0_Video.gif" height="300">

---

### 2D projection of the phase space

<img src="./images/Btwo103.20(period-3,5)/OmeThBtwo103.20_7.png" height="300">

---

### Basins of attraction

<img src="./images/Btwo103.20(period-3,5)/OmeZeThZeBtwo103.20_RK4_2.png" height="300">

---