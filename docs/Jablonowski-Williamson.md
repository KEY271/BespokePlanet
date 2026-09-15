# Jablonowski–Williamson の初期状態

Jablonowski–Williamson の傾圧不安定テストでは、乾燥大気を東西対称な平衡基本場で初期化する。定常性を調べる場合は基本場だけを使い、傾圧波を発達させる場合は東西風に局所的な擾乱を加える。

経度を $\lambda$、緯度を $\varphi$ とする。地表面気圧は一様に

$$
p_s(\lambda,\varphi)=p_0=1000\,\mathrm{hPa}=10^5\,\mathrm{Pa}
$$

とし、地表ジオポテンシャル $\Phi_s(\varphi)$ を[後述](#地表ジオポテンシャル)のように与える。

初期化に使う定数は

$$
\begin{aligned}
a&=6.371229\times10^6\,\mathrm{m},
&\Omega&=7.29212\times10^{-5}\,\mathrm{s^{-1}},\\
g&=9.80616\,\mathrm{m\,s^{-2}},
&R_d&=287.0\,\mathrm{J\,kg^{-1}\,K^{-1}},\\
u_0&=35\,\mathrm{m\,s^{-1}},
&\eta_0&=0.252,\\
T_0&=288\,\mathrm{K},
&\Gamma&=0.005\,\mathrm{K\,m^{-1}},\\
\eta_t&=0.2,
&\Delta T&=4.8\times10^5\,\mathrm{K}
\end{aligned}
$$

である。まず鉛直方向の補助変数を

$$
\eta_v=\frac{\pi}{2}(\eta-\eta_0)
$$

と置く。

## 鉛直座標

鉛直座標 $\eta$ は [hybrid-sigma 座標](./hybrid-sigma.md)の係数 $A_{k+1/2},B_{k+1/2}$ と

$$
\eta_{k+1/2}=\frac{A_{k+1/2}}{p_0}+B_{k+1/2},\qquad k=0,1,\dots,N
$$

で対応づける。これにより $p_s=p_0$ のとき各面で $p=\eta p_0$ となり、以下の式の $\eta$ は初期状態の $p/p_0$ に一致する。地表は $\eta_{N+1/2}=1$、モデル上端は $\eta_{1/2}=p_T/p_0$ である。なお、ここでの $A_{k+1/2},B_{k+1/2}$ は温度の節で使う緯度の関数 $A(\varphi),B(\varphi)$ とは別物である。

予報変数を置く full level では

$$
\eta_k=\frac{\eta_{k-1/2}+\eta_{k+1/2}}{2}
$$

とし、速度と温度はこの $\eta_k$ で評価する。$N=10$ の係数では

$$
\begin{aligned}
\eta_{k+1/2}&=0.01,0.05,0.10,0.18,0.28,0.40,0.52,0.65,0.77,0.88,1,\\
\eta_k&=0.03,0.075,0.14,0.23,0.34,0.46,0.585,0.71,0.825,0.94
\end{aligned}
$$

となる。$\eta_k$ を $(k-1/2)/N$ のように係数と無関係な等間隔に取ると、速度・温度を評価した $\eta$ とモデルの気圧が食い違い、基本場が釣り合わなくなる。

## 速度

東西対称な基本場は、東向き風速 $u$ と北向き風速 $v$ を

$$
\begin{aligned}
u_{\mathrm{base}}(\lambda,\varphi,\eta)
  &=u_0\cos^{3/2}(\eta_v)\sin^2(2\varphi),\\
v_{\mathrm{base}}(\lambda,\varphi,\eta)&=0
\end{aligned}
$$

として初期化する。この基本場は水平発散がゼロで、両半球の緯度 $45^\circ$、$\eta=\eta_0$ に最大風速 $u_0$ のジェットを持つ。

傾圧波を発達させる実験では、中心

$$
(\lambda_c,\varphi_c)=\left(\frac{\pi}{9},\frac{2\pi}{9}\right)
=(20^\circ\mathrm{E},40^\circ\mathrm{N})
$$

からの大円距離を

$$
\begin{aligned}
X&=\sin\varphi_c\sin\varphi
  +\cos\varphi_c\cos\varphi\cos(\lambda-\lambda_c),\\
r&=a\arccos X
\end{aligned}
$$

とし、全ての鉛直レベルで東西風に

$$
u'(\lambda,\varphi,\eta)
=u_p\exp\left[-\left(\frac{r}{R}\right)^2\right],
\qquad
u_p=1\,\mathrm{m\,s^{-1}},\qquad R=\frac{a}{10}
$$

を加える。したがって最終的な初期速度は

$$
u=u_{\mathrm{base}}+u',\qquad v=0
$$

である。平衡基本場の維持だけを調べる場合は $u'=0$ とする。浮動小数点の丸めで $X$ が定義域を外れないよう、実装では $\arccos$ の前に $X$ を $[-1,1]$ に収める。

実際の初期値は、この $u,v$ から $\zeta,\delta$ をスペクトル法で求めることで定める。

## 温度

温度の水平平均を

$$
\overline{T}(\eta)=
\begin{cases}
T_0\eta^{R_d\Gamma/g}, & \eta_t\le\eta\le1,\\[2mm]
T_0\eta^{R_d\Gamma/g}+\Delta T(\eta_t-\eta)^5,
  & 0\le\eta<\eta_t
\end{cases}
$$

とする。緯度依存部分を簡潔に書くため

$$
\begin{aligned}
A(\varphi)
  &=-2\sin^6\varphi\left(\cos^2\varphi+\frac{1}{3}\right)
    +\frac{10}{63},\\
B(\varphi)
  &=\frac{8}{5}\cos^3\varphi
    \left(\sin^2\varphi+\frac{2}{3}\right)-\frac{\pi}{4}
\end{aligned}
$$

と置けば、初期温度は

$$
\begin{aligned}
T(\lambda,\varphi,\eta)
=\overline{T}(\eta)
&+\frac{3}{4}\frac{\eta\pi u_0}{R_d}
  \sin\eta_v\cos^{1/2}\eta_v\\
&\quad\times\left[
  2u_0\cos^{3/2}\eta_v\,A(\varphi)
  +a\Omega B(\varphi)
  \right]
\end{aligned}
$$

である。温度は経度に依存せず、傾圧波を発生させる場合も温度擾乱は加えない。

## 地表ジオポテンシャル

基本場は地表 $\eta=1$ でも $u_0\cos^{3/2}\eta_{v,s}\sin^2(2\varphi)$（最大約 $8.4\,\mathrm{m\,s^{-1}}$）の東西風を持つ。地表面気圧を一様にしたままこれを釣り合わせるため、地表ジオポテンシャルを

$$
\Phi_s(\varphi)
=u_0\cos^{3/2}\eta_{v,s}
\left[u_0\cos^{3/2}\eta_{v,s}\,A(\varphi)+a\Omega B(\varphi)\right],
\qquad
\eta_{v,s}=\frac{\pi}{2}(1-\eta_0)
$$

と与える。これは基本場のジオポテンシャル

$$
\Phi(\varphi,\eta)=\overline\Phi(\eta)
+u_0\cos^{3/2}\eta_v
\left[u_0\cos^{3/2}\eta_v\,A(\varphi)+a\Omega B(\varphi)\right]
$$

の $\eta=1$ での値である。ただし水平平均部分 $\overline\Phi(\eta)$ は $\overline\Phi(1)=0$ を満たす。上の温度は静水圧の関係

$$
T=-\frac{\eta}{R_d}\frac{\partial\Phi}{\partial\eta}
$$

からこの $\Phi$ に対応して得られるものであり、$\Phi_s$ を省くと全層のジオポテンシャルから緯度依存部分が失われて地衡風バランスが崩れる。

$A(\varphi),B(\varphi)$ の全球平均はどちらも $0$ なので、$\Phi_s$ の全球平均も $0$ である。値は赤道で約 $1106\,\mathrm{m^2\,s^{-2}}$（高さ約 $113\,\mathrm{m}$）、極で約 $-3094\,\mathrm{m^2\,s^{-2}}$（約 $-316\,\mathrm{m}$）となる。

$\Phi_s$ は予報変数ではなく、時間変化しない下端の境界条件として固定する。傾圧波を発達させる実験でも同じ $\Phi_s$ を使う。実装では格子上で $\Phi_s$ を計算してからスペクトル変換し、他の変数と同じ切断波数で切った値を用いる。

## 参考文献

- C. Jablonowski and D. L. Williamson, “A baroclinic instability test case for atmospheric model dynamical cores,” *Quarterly Journal of the Royal Meteorological Society*, 132, 2943–2975 (2006), [doi:10.1256/qj.06.12](https://doi.org/10.1256/qj.06.12). 初期速度と温度は印刷頁 2946–2947、風速擾乱は印刷頁 2948–2949 を参照。
