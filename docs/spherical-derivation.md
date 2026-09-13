# 球面上の微分演算子

半径 $a$ の球面上で、経度を $\lambda$、緯度を $\varphi$ とする。東向き・北向きの単位ベクトルをそれぞれ $\boldsymbol{e}_\lambda,\boldsymbol{e}_\varphi$、外向き単位法線を $\boldsymbol{k}$ とし、接ベクトルの成分を

$$
\boldsymbol{A}=A_\lambda\boldsymbol{e}_\lambda+A_\varphi\boldsymbol{e}_\varphi
$$

と書く。基底の向きは

$$
\boldsymbol{k}\times\boldsymbol{e}_\lambda=\boldsymbol{e}_\varphi,
\qquad
\boldsymbol{k}\times\boldsymbol{e}_\varphi=-\boldsymbol{e}_\lambda
$$

とする。

## 計量とスケール因子

球面上の線素と面積要素は

$$
ds^2=a^2\cos^2\varphi\,d\lambda^2+a^2\,d\varphi^2,
\qquad
dS=a^2\cos\varphi\,d\lambda\,d\varphi
$$

である。したがって直交曲線座標 $(\lambda,\varphi)$ のスケール因子は

$$
h_\lambda=a\cos\varphi,
\qquad
h_\varphi=a
$$

となる。以下の各演算子は、このスケール因子を直交曲線座標の公式へ代入すれば得られる。

## 勾配と90度回転

スカラー場 $A$ の球面勾配は

$$
\boldsymbol{\nabla}_s A
=\boldsymbol{e}_\lambda\frac{1}{a\cos\varphi}\frac{\partial A}{\partial\lambda}
+\boldsymbol{e}_\varphi\frac{1}{a}\frac{\partial A}{\partial\varphi}
$$

である。これを接平面内で反時計回りに90度回転すると

$$
\boldsymbol{k}\times\boldsymbol{\nabla}_s A
=-\boldsymbol{e}_\lambda\frac{1}{a}\frac{\partial A}{\partial\varphi}
+\boldsymbol{e}_\varphi\frac{1}{a\cos\varphi}\frac{\partial A}{\partial\lambda}
$$

となる。

## 発散・回転・ラプラシアン

接ベクトル場 $\boldsymbol{A}$ の発散は

$$
\boldsymbol{\nabla}_s\cdot\boldsymbol{A}
=\frac{1}{a\cos\varphi}
\left[
\frac{\partial A_\lambda}{\partial\lambda}
+\frac{\partial(A_\varphi\cos\varphi)}{\partial\varphi}
\right]
$$

であり、回転の球面に垂直な成分は

$$
\boldsymbol{k}\cdot(\boldsymbol{\nabla}_s\times\boldsymbol{A})
=\frac{1}{a\cos\varphi}
\left[
\frac{\partial A_\varphi}{\partial\lambda}
-\frac{\partial(A_\lambda\cos\varphi)}{\partial\varphi}
\right]
$$

である。スカラー場に作用するラプラシアンは、発散と勾配の合成として

$$
\begin{aligned}
\nabla_s^2 A
&=\boldsymbol{\nabla}_s\cdot\boldsymbol{\nabla}_s A\\
&=\frac{1}{a^2\cos^2\varphi}\frac{\partial^2 A}{\partial\lambda^2}
+\frac{1}{a^2\cos\varphi}\frac{\partial}{\partial\varphi}
\left(\cos\varphi\frac{\partial A}{\partial\varphi}\right)
\end{aligned}
$$

となる。

## Helmholtz 分解

球面上の滑らかな接ベクトル場 $\boldsymbol{u}=u\boldsymbol{e}_\lambda+v\boldsymbol{e}_\varphi$ は、流線関数 $\psi$ と速度ポテンシャル $\chi$ を使って

$$
\boldsymbol{u}
=\boldsymbol{k}\times\boldsymbol{\nabla}_s\psi
+\boldsymbol{\nabla}_s\chi
$$

と分解できる。成分表示では

$$
\begin{aligned}
u&=-\frac{1}{a}\frac{\partial\psi}{\partial\varphi}
+\frac{1}{a\cos\varphi}\frac{\partial\chi}{\partial\lambda},\\
v&=\frac{1}{a\cos\varphi}\frac{\partial\psi}{\partial\lambda}
+\frac{1}{a}\frac{\partial\chi}{\partial\varphi}
\end{aligned}
$$

である。相対渦度 $\zeta$ と発散 $\delta$ はそれぞれ

$$
\begin{aligned}
\zeta
&=\boldsymbol{k}\cdot(\boldsymbol{\nabla}_s\times\boldsymbol{u})
=\nabla_s^2\psi,\\
\delta
&=\boldsymbol{\nabla}_s\cdot\boldsymbol{u}
=\nabla_s^2\chi
\end{aligned}
$$

となる。特に非発散流では $\chi=0$ としてよい。球面では $H^1(S^2)=0$ なので調和 $1$-形式の項はなく、$\psi,\chi$ は加法定数を除いて定まる。また閉じた球面上では $\zeta$ と $\delta$ の面積積分はいずれもゼロである。
