# hybrid-sigma 座標

## 基礎方程式

地球半径 $a$ は大気の厚さより十分大きいので、鉛直座標 $z$ 一定の球面微分を

$$
\bm\nabla_z=\bm{e}_\lambda\frac{1}{a\cos\varphi}\frac{\partial}{\partial\lambda}+\bm{e}_\varphi\frac{1}{a}\frac{\partial}{\partial\varphi}
$$

として、Navier-Stokes 方程式と連続の式を

$$
\frac{D\bm{u}}{Dt}+f\bm{k}\times\bm{u}=-\frac{1}{\rho}\bm\nabla_zp+\bm{F},\\
\frac{Dw}{Dt}=-\frac{1}{\rho}\frac{\partial p}{\partial z}-g\\
\frac{\partial\rho}{\partial t}+\bm\nabla_z\cdot(\rho\bm{u})+\frac{\partial(\rho w)}{\partial z}=0
$$

と近似する。さらに静水圧近似として、$Dw/Dt=0$ を課す。これにより

$$
\frac{\partial p}{\partial z}=-\rho g
$$

となるので $dp=-\rho gdz$ が成り立つ。

## hybrid-sigma 座標

$z$ の代わりに地表面気圧 $p_s$ と関数 $A(\eta),B(\eta)$ により

$$
p(\eta)=A(\eta)+B(\eta)p_s
$$

と定義される $\eta$ を鉛直座標に用いる。ただし $0\le\eta\le1$ で、$A(1)=0,B(1)=1$ で $\eta=1$ が地表になるようにする。また $B(0)=0$ とし、上端気圧 $p_T=A(0)$ は時間依存しないものとする。

### 連続の式
$d\eta$ の中にある質量は

$$
dm=\rho |dxdydz|=\rho\left|\frac{\partial z}{\partial p} dxdydp\right|=\frac{1}{g}|dxdydp|=\frac{1}{g}\frac{\partial p}{\partial\eta}|dxdyd\eta|
$$

となるので、

$$
\mu=\frac{\partial p}{\partial\eta}=A'(\eta)+B'(\eta)p_s
$$

と置けば連続の式は

$$
\frac{\partial\mu}{\partial t}+\bm\nabla_\eta\cdot(\mu\bm{u})+\frac{\partial(\mu\dot\eta)}{\partial\eta}=0
$$

となる。ただし

$$
\dot\eta=\frac{D\eta}{Dt}
$$

である。鉛直質量フラックスを $M=\mu\dot\eta$ と置く。

### 地表面気圧

連続の式を $\eta=0$ から $1$ まで積分すると、上端と地表で空気が突き抜けないことから $\dot\eta(0)=\dot\eta(1)=0$ なので

$$
\frac{\partial p_s}{\partial t}+\bm\nabla_\eta\cdot\left(\int_0^1\mu\bm{u}d\eta\right)=0
$$

つまり、

$$
\frac{\partial p_s}{\partial t}=-\bm\nabla_\eta\cdot\left(\int_{p_T}^{p_s}\bm{u}dp\right)
$$

となる。

### 鉛直質量フラックス

連続の式を $\eta=0$ から $\eta$ まで積分した場合は、

$$
\frac{\partial p(\eta)}{\partial t}+\bm\nabla_\eta\cdot\left(\int_{p_T}^{p(\eta)}\bm{u}dp\right)+M(\eta)=0
$$

となるので、hybrid-sigma 座標の式を代入して

$$
M(\eta)=-B(\eta)\frac{\partial p_s}{\partial t}-\bm\nabla\cdot\left(\int_{p_T}^{p(\eta)}\bm{u}dp\right)
$$

が得られる。

物質微分

$$
\frac{D}{Dt}=\frac{\partial}{\partial t}+\bm{u}\cdot\bm\nabla_\eta+\dot\eta\frac{\partial}{\partial\eta}
$$

の中の鉛直移流成分は

$$
\dot\eta\frac{\partial}{\partial\eta}=\frac{M}{\mu}\frac{\partial}{\partial\eta}
$$

となるので $M$ から得られる。

### ジオポテンシャル

ジオポテンシャル

$$
\Phi=gz
$$

を導入すると、

$$
\frac{\partial\Phi}{\partial p}=g\frac{\partial z}{\partial p}=-\frac{1}{\rho}
$$

となる。理想気体を仮定すると右辺は $-RT/p$ となるので、

$$
\frac{\partial\Phi}{\partial\ln p}=-RT
$$

を得る。地表面の高さを $z_s$ とし、$\Phi_s=gz_s$ からこれを積分して

$$
\Phi=\Phi_s+\int_{\ln p(\eta)}^{\ln p_s}RT d\ln p
$$

としてジオポテンシャルが得られる。

### 水平運動方程式

偏微分を計算すると

$$
-\frac{1}{\rho}\bm\nabla_zp
=-\bm\nabla_\eta\Phi-\frac{1}{\rho}\bm\nabla_\eta p
$$

と変形できるので、水平方向の運動方程式は

$$
\frac{D\bm{u}}{Dt}+f\bm{k}\times\bm{u}=-\bm\nabla_\eta\Phi-\frac{1}{\rho}\bm\nabla_\eta p+\bm{F}
$$

となる。物質微分を展開すれば

$$
\frac{\partial\bm{u}}{\partial t}=-(\bm{u}\cdot\bm\nabla_\eta)\bm{u}-\frac{M}{\mu}\frac{\partial\bm{u}}{\partial\eta}-f\bm{k}\times\bm{u}-\bm\nabla_\eta\Phi-\frac{RT}{p}\bm\nabla_\eta p+\bm{F}
$$

を得る。このとき運動エネルギー $K=\frac{1}{2}(u^2+v^2)$ を使うと

$$
(\bm{u}\cdot\bm\nabla_\eta)\bm{u}=\bm\nabla_\eta K+\zeta\bm{k}\times\bm{u}
$$

となるので、

$$
\frac{\partial\bm{u}}{\partial t}=-\bm\nabla_\eta(K+\Phi)-\frac{M}{\mu}\frac{\partial\bm{u}}{\partial\eta}-(\zeta+f)\bm{k}\times\bm{u}-\frac{RT}{p}\bm\nabla_\eta p+\bm{F}
$$

となる。

### 渦度

水平運動方程式に $\bm{k}\cdot\bm\nabla_\eta\times$ を作用させて、渦度の方程式

$$
\frac{\partial\zeta}{\partial t}=-\bm\nabla_\eta\cdot[(\zeta+f)\bm{u}]+\bm{k}\cdot\bm\nabla_\eta\times\left[-\frac{M}{\mu}\frac{\partial\bm{u}}{\partial\eta}-RT\bm\nabla_\eta\ln p+\bm{F}\right]
$$

を得る。

### 発散

水平運動方程式に $\bm\nabla\cdot$ を作用させて、発散の方程式

$$
\frac{\partial\delta}{\partial t}=-\nabla_\eta^2(K+\Phi)+\bm{k}\cdot\bm\nabla_\eta\times[(\zeta+f)\bm{u}]+\bm\nabla_\eta\cdot\left[-\frac{M}{\mu}\frac{\partial\bm{u}}{\partial\eta}-RT\bm\nabla_\eta\ln p+\bm{F}\right]
$$

となる。

### 温度

熱力学第一法則により、乾燥断熱の場合は

$$
\frac{DT}{Dt}=\kappa\frac{T}{p}\frac{Dp}{Dt}
$$

となる。ただし $\kappa=2/7$ である。左辺の物質微分を展開すれば

$$
\frac{\partial T}{\partial t}=-\bm{u}\cdot\bm\nabla_\eta T-\frac{M}{\mu}\frac{\partial T}{\partial\eta}+\kappa T\frac{D\ln p}{Dt}
$$

が得られる。

## A,B

$p_0=1000\,\mathrm{hPa}$ を基準とし、上端気圧を $p_{1/2}=10\,\mathrm{hPa}$ とする。$N=10$ では $A,B$ は

$$
\begin{aligned}
A_{k+1/2}&=1000,5000,10000,8000,8000,10000,12000,10000,7000,3000,0\,\mathrm{Pa}\\
B_{k+1/2}&=0,0,0,0.1,0.2,0.3,0.4,0.55,0.7,0.85,1
\end{aligned}
$$

と置く。離散化した各面の $\eta$ は

$$
\eta_{k+1/2}=\frac{A_{k+1/2}}{p_0}+B_{k+1/2}
$$

と対応づけ、$p_s=p_0$ のとき $p_{k+1/2}=\eta_{k+1/2}p_0$ となるようにする。上端は $\eta_{1/2}=p_T/p_0=0.01$ である。予報変数を置く full level では $\eta_k=(\eta_{k-1/2}+\eta_{k+1/2})/2$ とする。
