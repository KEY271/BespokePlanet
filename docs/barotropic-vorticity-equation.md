# 順圧渦度方程式

## 方程式

地球半径 $a=6371\,\mathrm{km}$ の球面上の水平流が鉛直方向に変化しないと仮定する。経度を $\lambda$、緯度を $\varphi$、東向き・北向き風速をそれぞれ $u,v$ とすると、非圧縮流なら $u,v$ は流線関数 $\psi$ によって

$$
u=-\frac{1}{a}\frac{\partial\psi}{\partial\varphi},\qquad
v=\frac{1}{a\cos\varphi}\frac{\partial\psi}{\partial\lambda}
$$

と表せる。これは非発散流は閉 $1$-形式とみなせるが、$H^1(S^2)=0$ だからである。相対渦度 $\zeta$ を

$$
\begin{aligned}
\zeta=\frac{1}{a\cos\varphi}\left[\frac{\partial v}{\partial\lambda}-\frac{\partial(u\cos\varphi)}{\partial\varphi}\right]=\nabla_s^2\psi
\end{aligned}
$$

と定義する。ただし球面ラプラシアンは

$$
\nabla_s^2 A
=\frac{1}{a^2\cos^2\varphi}\frac{\partial^2 A}{\partial\lambda^2}
+\frac{1}{a^2\cos\varphi}\frac{\partial}{\partial\varphi}
\left(\cos\varphi\frac{\partial A}{\partial\varphi}\right)
$$

である。自転角速度を $\Omega=7.2921159\times10^{-5}\,\mathrm{rad/s}$ とし、惑星渦度を

$$
f=2\Omega\sin\varphi
$$

と定義し、絶対渦度を $q=\zeta+f$ と置く。外力も散逸もない順圧渦度方程式は

$$
\frac{\partial\zeta}{\partial t}+J(\psi,q)=0
$$

となる。ただし球面上のヤコビアンを

$$
J(A,B)=\frac{1}{a^2\cos\varphi}
\left(
\frac{\partial A}{\partial\lambda}\frac{\partial B}{\partial\varphi}
-\frac{\partial A}{\partial\varphi}\frac{\partial B}{\partial\lambda}
\right)
$$

とした。したがって方程式は絶対渦度の物質保存則

$$
\frac{Dq}{Dt}=0,
\qquad
\frac{D}{Dt}
=\frac{\partial}{\partial t}
+\frac{u}{a\cos\varphi}\frac{\partial}{\partial\lambda}
+\frac{v}{a}\frac{\partial}{\partial\varphi}
$$

に等しい。

## 格子解法

格子上では次のステップで解く。

1. 相対渦度 $\zeta$ が既知であるとする。
2. 流線関数 $\psi$ を $\zeta=\nabla_s^2\psi$ で求める。
3. $u,v$ を $\psi$ から求める。
4. 格子に戻して $J(\psi,q)$ を求める。
5. 順圧渦度方程式で $\zeta$ を時間発展させる。

## スペクトル法

流線関数 $\psi$ を $\zeta$ から求めるのは球面調和関数の展開係数では非常に簡単であり、単に

$$
\psi_n^m=-\frac{a^2}{n(n+1)}\zeta_n^m,\qquad n\geq1
$$

とすればよい。

$n=0$ は空間的に一定な流線関数であり、速度に寄与しない。そこでゲージとして $\psi_0^0=0$ と置く。また、閉じた球面上では相対渦度の面積積分はゼロであるため、$\zeta_0^0$ もゼロでなければならない。丸め誤差で生じた $\zeta_0^0$ は各ステップでゼロに戻す。

$u,v$ は [Octahedral Gaussian Grid のスペクトル微分](./octahedral-gaussian-grid.md#スペクトル微分)で求まる。

2次の非線形項 $J(\psi,q)$ については、非発散性を利用すると

$$
J(\psi,q)
=\nabla_s\cdot(\boldsymbol{v}q)
=\frac{1}{a\cos\varphi}\left[\frac{\partial(uq)}{\partial\lambda}+\frac{\partial(vq\cos\varphi)}{\partial\varphi}\right]
$$

と書ける。$uq,vq\cos\varphi$ を格子で求めてからスペクトルに戻して微分すればよい。

## Leapfrog 法

空間離散化後は

$$
\frac{d\boldsymbol{\zeta}}{dt}=R(\boldsymbol{\zeta})
$$

という常微分方程式になる。時間積分は Robert--Asselin--Williams（RAW）フィルター付き Leapfrog 法を用いる。RAW フィルターはフィルターされていない値 $\zeta^k$ とフィルターされた値 $\overline\zeta^{k-1}$ を用いて $\overline\zeta^k,\zeta^{k+1}$ を決める。

### 超粘性

格子スケール付近の渦が大きくなっていかないように、人工的に超粘性を加える。モード $(n,m)$ の減衰率を

$$
\kappa_n
=\nu_{2p}
\left(\frac{n(n+1)}{a^2}\right)^p
$$

とする。LeapFrog の通常ステップに超粘性を加えたものは

$$
(\zeta_n^m)^{k+1}
=(\overline\zeta_n^m)^{k-1}
+2\Delta t
\left[
R_n^m(\zeta^k)
-\kappa_n\left(\zeta_n^m\right)^{k+1}
\right]
$$

なので、これを解いて

$$
\left(\zeta_n^m\right)^{k+1}_*
=\frac{(\overline\zeta_n^m)^{k-1}+2\Delta t\,R_n^m(\zeta^k)}{1+2\Delta t\,\kappa_n}
$$

とする。実際にこれを計算すると $\kappa_n$ が大きくなってしまうので正規化する。$n_\mathrm{max}=T$ として

$$
\left(\zeta_n^m\right)^{k+1}_*
=\frac{(\overline\zeta_n^m)^{k-1}+2\Delta t\,R_n^m(\zeta^k)}{1+2\dfrac{\Delta t}{\tau}\left(\dfrac{n(n+1)}{n_\mathrm{max}(n_\mathrm{max}+1)}\right)^p}
$$

とする。時定数 $\tau$ は $6$ 時間とし、$p=4$ とする。

### RAW フィルター

Leapfrog 法には物理解と独立な、時間ステップごとに符号が反転する計算モードがある。このモードを抑えるため RAW フィルターを適用する。まず

$$
d^k=\frac{\epsilon}{2}
\left(
\overline{\zeta}^{\,k-1}
-2\zeta^k
+\zeta_*^{k+1}
\right)
$$

を計算し、

$$
\begin{aligned}
\overline{\zeta}^{\,k}
&=\zeta^k+\alpha d^k,\\
\zeta^{k+1}
&=\zeta_*^{k+1}-(1-\alpha)d^k
\end{aligned}
$$

と補正する。$\epsilon$ が Robert フィルター係数、$\alpha$ が Williams フィルター係数である。ここでは

$$
\epsilon=0.1,\qquad \alpha=0.53
$$

とする。

### 初期化

LeapFrog は二時刻を必要とするので、初期値 $\zeta^0$ から次の二段階で $\zeta^1$ を作る。

$$
\begin{aligned}
(\zeta^m_n)^{1/2}
&=(\zeta^m_n)^0+\frac{\Delta t}{2}\left[R((\zeta^m_n)^0)-\kappa_n(\zeta^m_n)^{1/2}\right],\\
(\zeta^m_n)^1
&=(\zeta^m_n)^0+\Delta t\left[R((\zeta^m_n)^{1/2})-\kappa_n(\zeta^m_n)^1\right].
\end{aligned}
$$

一段目は半ステップの Euler 法、二段目は中点の傾向を使う最初の Leapfrog 更新である。これを解いて

$$
\begin{aligned}
(\overline\zeta^m_n)^0&=(\zeta^m_n)^0,\\
(\zeta^m_n)^{1/2}&=\frac{(\zeta^m_n)^0+\dfrac{\Delta t}{2}R((\zeta^m_n)^0)}{1+\kappa_n\dfrac{\Delta t}{2}},\\
(\zeta^m_n)^1&=\frac{(\zeta^m_n)^0+\Delta t\,R((\zeta^m_n)^{1/2})}{1+\kappa_n\Delta t}
\end{aligned}
$$

として LeapFrog を始める。

### 時間ステップ

最大速度 $U_\mathrm{max}=50\,\mathrm{m}/\mathrm{s}$ の流れを考える。移流 CFL 条件は

$$
\Delta t\lesssim C\frac{a}{U_\mathrm{max}\sqrt{T(T+1)}}
$$

なので $T=63$ で $C=0.5$ とすると、$\Delta t=900\,\mathrm{s}$ くらいで良いことになる。各ステップで

$$
C=\frac{U_\mathrm{max}\Delta t}{a}\sqrt{T(T+1)}
$$

を計算して $0.5$ 程度以下かを監視する。
