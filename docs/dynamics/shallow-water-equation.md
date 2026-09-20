# 浅水方程式

地球半径 $a=6371\,\mathrm{km}$ の球面上の浅水波を考える。座標、基底、球面上の微分演算子の定義は [球面上の微分演算子](./spherical-derivation.md)に従う。一般の流れは流線関数 $\psi$ と速度ポテンシャル $\chi$ によって

$$
\boldsymbol{u}
=\boldsymbol{k}\times\boldsymbol{\nabla}_s\psi
+\boldsymbol{\nabla}_s\chi
$$

と表せる。相対渦度 $\zeta$ と発散 $\delta$ は

$$
\begin{aligned}
\zeta&=\boldsymbol{k}\cdot(\boldsymbol{\nabla}_s\times\boldsymbol{u})=\nabla_s^2\psi,\\
\delta&=\boldsymbol{\nabla}_s\cdot\boldsymbol{u}=\nabla_s^2\chi
\end{aligned}
$$

と置く。

自転角速度を $\Omega=7.2921159\times10^{-5}\,\mathrm{rad/s}$ と置き、惑星渦度を

$$
f=2\Omega\sin\varphi
$$

とする。$h$ を波の高さとし、$g=9.8\,\mathrm{m/s^2}$ を重力加速度とする。また、運動エネルギーを

$$
K=\frac{1}{2}(u^2+v^2)
$$

と置く。

浅水方程式は

$$
\begin{aligned}
\frac{\partial\zeta}{\partial t}&=-\boldsymbol{\nabla}_s\cdot((\zeta+f)\boldsymbol{u}),\\
\frac{\partial\delta}{\partial t}&=-\boldsymbol{\nabla}_s\cdot((\zeta+f)\boldsymbol{k}\times\boldsymbol{u})-\nabla_s^2(gh+K),\\
\frac{\partial h}{\partial t}&=-\boldsymbol{\nabla}_s\cdot(h\boldsymbol{u})
\end{aligned}
$$

と書ける。

## 重力波

波の高さを平均 $H=10\,\mathrm{km}$ と偏差 $\eta$ に分けて $h=H+\eta$ とする。このとき浅水方程式は非常に高速な重力波のモードを含む。変数をまとめて

$$
X=\begin{pmatrix}\zeta\\\delta\\\eta\end{pmatrix}
$$

と書くと重力波に関係する部分を

$$
GX=\begin{pmatrix}0\\-g\nabla_s^2\eta\\-H\delta\end{pmatrix}\\
$$

として方程式は

$$
\frac{\partial X}{\partial t}=N(X)+GX
$$

という形になる。非線形項は

$$
N(X)=\begin{pmatrix}
-\bm\nabla_s\cdot((\zeta+f)\bm{u})\\
-\bm\nabla_s\cdot((\zeta+f)\bm{k}\times\bm{u})-\nabla_s^2K\\
-\bm\nabla_s\cdot(\eta\bm{u})
\end{pmatrix}
$$

である。

LeapFrog法では、重力波は陰的に解く。陰的にする度合いを $\beta=0.5$ として

$$
\frac{X^{k+1}_*-\overline X^{k-1}}{2\Delta t}=N(X^k)+(1-\beta)G\overline{X}^{k-1}+\beta GX^{k+1}_*
$$

とする。つまり

$$
(1-2\beta\Delta t\,G)X^{k+1}_*=\overline{X}^{k-1}+2\Delta t[N(X^k)+(1-\beta)G\overline{X}^{k-1}]
$$

として $X_*^{k+1}$ を得る。$\zeta$ に重力波は関係しないので、左辺は実質的に $2\times2$ 行列がかかっていることになり、この方程式は解析的に解ける。

時間ステップは $\Delta t=1200\,\mathrm{s}$ とする。

## 超粘性

$\zeta,\delta$ にだけ超粘性を入れる。モード $(n,m)$ の減衰率を

$$
\kappa_n=\frac{1}{\tau}\left(\frac{n(n+1)}{T(T+1)}\right)^p
$$

とし、$\tau$ を $6$ 時間、$p=4$ とする。超粘性は

$$
\frac{\zeta^{k+1}_{**}-\zeta^{k+1}_*}{2\Delta t}=-\kappa_n\zeta^{k+1}_{**}
$$

のように入れ、

$$
\zeta^{k+1}_{**}=\frac{\zeta^{k+1}_*}{1+2\Delta t\,\kappa_n}
$$

とする。$\delta$ についても同様である。$\eta$ は単に

$$
\eta^{k+1}_{**}=\eta^{k+1}_*
$$

とする。

## RAW フィルター

RAW フィルターは順圧渦度方程式と同様に

$$
d^k=\frac{\epsilon}{2}(\overline X^{k-1}-2X^k+X^{k+1}_{**})
$$

と定義し、

$$
\begin{aligned}
\overline X^k&=X^k+\alpha d^k,\\
X^{k+1}&=X^{k+1}_{**}-(1-\alpha)d^k
\end{aligned}
$$

とする。ここでは

$$
\epsilon=0.1,\quad\alpha=0.53
$$

と置く。

このようにして更新した後、$\zeta^0_0=\delta^0_0=\eta^0_0=0$ と手で置く。

## 初期化

上の手順をまとめれば結局、

$$
X^{k+1}=F(\Delta t,\overline{X}^{k-1},X^k)
$$

としていることになる。初期化の際には、初期値 $X^0$ から

$$
\begin{aligned}
\overline{X}^0&=X^0,\\
X^{1/2}&=F_\text{noRAW}\left(\frac{\Delta t}{4},X^0,X^0\right),\\
X^1&=F\left(\frac{\Delta t}{2},X^0,X^{1/2}\right)
\end{aligned}
$$

として LeapFrog を始める。ただし $X^{1/2}$ だけは RAW フィルターを使わない。
