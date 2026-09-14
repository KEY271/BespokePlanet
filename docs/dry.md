# 乾燥大気

## 予報変数と基礎方程式

各ステップで時間発展させる変数は、相対渦度 $\zeta$ と発散 $\delta$、温度 $T$、地表面気圧 $\ln p_s$ である。以下は [hybrid-sigma 座標](./hybrid-sigma.md)の記号に従う。

基礎方程式は

$$
\begin{aligned}
\frac{\partial\ln p_s}{\partial t}&=-\frac{1}{p_s}\bm\nabla_\eta\cdot\left(\int_{p_T}^{p_s}\bm{u}dp\right)\\
\frac{\partial\zeta}{\partial t}&=-\bm\nabla_\eta\cdot[(\zeta+f)\bm{u}]+\bm{k}\cdot\bm\nabla_\eta\times\left[-\frac{M}{\mu}\frac{\partial\bm{u}}{\partial\eta}-RT\bm\nabla_\eta\ln p+\bm{F}\right]\\
\frac{\partial\delta}{\partial t}&=-\nabla_\eta^2(K+\Phi)+\bm{k}\cdot\bm\nabla_\eta\times[(\zeta+f)\bm{u}]+\bm\nabla_\eta\cdot\left[-\frac{M}{\mu}\frac{\partial\bm{u}}{\partial\eta}-RT\bm\nabla_\eta\ln p+\bm{F}\right]\\
\frac{\partial T}{\partial t}&=-\bm{u}\cdot\bm\nabla_\eta T-\frac{M}{\mu}\frac{\partial T}{\partial\eta}+\kappa T\frac{D\ln p}{Dt}
\end{aligned}
$$

である。右辺に出てくる量は

$$
\begin{aligned}
p(\eta)&=A(\eta)+B(\eta)p_s,\\
\mu&=\frac{\partial p}{\partial\eta},\\
M&=-Bp_s\frac{\partial\ln p_s}{\partial t}-\bm\nabla_\eta\cdot\left(\int_{p_T}^{p(\eta)}\bm{u}dp\right),\\
\Phi&=\Phi_s+\int_{\ln p(\eta)}^{\ln p_s}RTd\ln p
\end{aligned}
$$

として求まる。

## 鉛直方向の離散化

普通の変数 $\zeta,\delta,T$ は $k=1,2,\dots,N$ に対して $\zeta_k,\delta_k,T_k$ と離散化する。

### 圧力

圧力を

$$
p_{k+1/2}=A_{k+1/2}+B_{k+1/2}p_s,\quad k=0,1,\dots,N
$$

と置き、

$$
A_{1/2}=p_T,\quad B_{1/2}=0,\quad A_{N+1/2}=0,\quad B_{N+1/2}=1
$$

とする。つまり $p_{1/2}=p_T,p_{N+1/2}=p_s$ と固定する。そして $k=1,2,\dots,N$ に対して

$$
L_k=\ln\frac{p_{k+1/2}}{p_{k-1/2}},\quad\Delta p_k=p_{k+1/2}-p_{k-1/2},\quad\alpha_k=1-\frac{p_{k-1/2}}{\Delta p_k}L_k
$$

と置く。

### ジオポテンシャル

ジオポテンシャルは

$$
\Phi_{N+1/2}=\Phi_s
$$

から始めて、

$$
\Phi_{k-1/2}=\Phi_{k+1/2}+RT_kL_k
$$

と置く。そして層内で平均を取って

$$
\Phi_k=\Phi_{k+1/2}+\alpha_kRT_k
$$

とする。

### 圧力傾度力

次に、方程式の右辺に出てくる $\bm\nabla_\eta\ln p$ については

$$
(\bm\nabla_\eta\ln p)_k=\bm{G}_k=\frac{1}{\Delta p_k}[L_k\bm\nabla_\eta p_{k-1/2}+\alpha_k\bm\nabla_\eta\Delta p_k]
$$

と計算する。これは層内で圧力を線形補間して積分すると得られる。hybrid-sigma 座標を代入すると、

$$
\bm{G}_k=\frac{p_s}{\Delta p_k}(B_{k-1/2}L_k+\alpha_k\Delta B_k)\bm\nabla_\eta \ln p_s
$$

となる。

### 地表面気圧と鉛直質量フラックス

まず

$$
\begin{aligned}
F_k&=\bm\nabla_\eta\cdot(\bm{u}_k\Delta p_k)\\
&=\Delta p_k\delta_k+p_s\Delta B_k\bm{u}_k\cdot\bm\nabla_\eta\ln p_s
\end{aligned}
$$

と置く。これの累積を

$$
S_k=\sum_{j=1}^kF_j,\quad S_0=0
$$

と置くと地表面気圧の時間微分は

$$
\frac{\partial\ln p_s}{\partial t}=-\frac{S_N}{p_s}
$$

と与えられる。鉛直質量フラックスは

$$
\begin{aligned}
M_{k+1/2}&=-B_{k+1/2}p_s\frac{\partial\ln p_s}{\partial t}-S_k\\
&=B_{k+1/2}S_N-S_k
\end{aligned}
$$

とする。$F_k$ の定義により、

$$
M_{1/2}=M_{N+1/2}=0
$$

が成立する。

### 鉛直移流

鉛直方向の移流から来る

$$
\frac{M}{\mu}\frac{\partial X}{\partial\eta}
$$

という形の項については、

$$
W_k(X)=\frac{1}{2\Delta p_k}[M_{k+1/2}(X_{k+1}-X_k)+M_{k-1/2}(X_k-X_{k-1})]
$$

と離散化する。上下端は

$$
W_1(X)=\frac{M_{3/2}(X_2-X_1)}{2\Delta p_1},\quad W_N(X)=\frac{M_{N-1/2}(X_N-X_{N-1})}{2\Delta p_N}
$$

となる。

### 温度の方程式

温度の方程式に出てくる $D(\ln p)/Dt$ については連続理論で

$$
\begin{aligned}
\frac{D\ln p}{Dt}&=\frac{\partial \ln p}{\partial t}+\bm{u}\cdot\bm\nabla_\eta\ln p+\dot\eta\frac{\partial\ln p}{\partial\eta}\\
&=\frac{1}{p}\left[B\frac{\partial p_s}{\partial t}+B\bm{u}\cdot\bm\nabla_\eta p_s+M\right]\\
&=\frac{1}{p}\left[B\bm{u}\cdot\bm\nabla_\eta p_s-\bm\nabla_\eta\cdot\left(\int_{p_T}^p\bm{u}dp\right)\right]
\end{aligned}
$$

となるので、

$$
\left(\frac{D\ln p}{Dt}\right)_k=\bm{u}_k\cdot\bm{G}_k-\frac{L_kS_{k-1}+\alpha_kF_k}{\Delta p_k}
$$

と離散化する。

### 離散化された方程式

結局最終的に、離散化された方程式は

$$
\begin{aligned}
\frac{\partial\ln p_s}{\partial t}&=-\frac{S_N}{p_s},\\
\frac{\partial\zeta_k}{\partial t}&=-\bm\nabla_\eta\cdot[(\zeta_k+f)\bm{u}_k]+\bm{k}\cdot\bm\nabla_\eta\times[-W_k(\bm{u})-RT_k\bm{G}_k+\bm{F}_k],\\
\frac{\partial\delta_k}{\partial t}&=-\nabla_\eta^2(K_k+\Phi_k)+\bm{k}\cdot\bm\nabla_\eta\times[(\zeta_k+f)\bm{u}_k]+\bm\nabla_\eta\cdot[-W_k(\bm{u})-RT_k\bm{G}_k+\bm{F}_k],\\
\frac{\partial T_k}{\partial t}&=-\bm{u}_k\cdot\bm\nabla_\eta T_k-W_k(T)+\kappa T_k Q_k
\end{aligned}
$$

となる。外力 $\bm{F}_k$ は一応書いてあるが、実際にはとりあえず $0$ にする。ただし

$$
\begin{aligned}
p_{k+1/2}&=A_{k+1/2}+B_{k+1/2}p_s,\\
\Delta p_k&=p_{k+1/2}-p_{k-1/2},\\
L_k&=\ln\frac{p_{k+1/2}}{p_{k-1/2}},\\
\alpha_k&=1-\frac{p_{k-1/2}}{\Delta p_k}L_k,\\
F_k&=\Delta p_k\delta_k+p_s\Delta B_k\bm{u}_k\cdot\bm\nabla_\eta\ln p_s,\\
S_k&=\sum_{j=1}^kF_j,\\
M_{k+1/2}&=B_{k+1/2}S_N-S_k,\\
\bm{G}_k&=\frac{p_s}{\Delta p_k}(B_{k-1/2}L_k+\alpha_k\Delta B_k)\bm\nabla_\eta\ln p_s,\\
Q_k&=\bm{u}_k\cdot\bm{G}_k-\frac{L_kS_{k-1}+\alpha_kF_k}{\Delta p_k}
\end{aligned}
$$

であり、$\Phi_k$ は

$$
\begin{aligned}
\Phi_{N+1/2}&=\Phi_s,\\
\Phi_{k-1/2}&=\Phi_{k+1/2}+RT_kL_k,\\
\Phi_k&=\Phi_{k+1/2}+\alpha_kRT_k
\end{aligned}
$$

と計算し、$W_k(X)$ は

$$
W_k(X)=\frac{1}{2\Delta p_k}[M_{k+1/2}(X_{k+1}-X_k)+M_{k-1/2}(X_k-X_{k-1})]
$$

とする。

## 水平方向の離散化

水平方向の離散化は[浅水方程式](./shallow-water-equation.md)と同じように行う。つまり $\zeta,\delta$ から流線関数 $\psi$ と速度ポテンシャル $\chi$ を出し、そこから $u,v$ を出す。

## 重力波

基準大気からのずれを重力波として分離して陰的に解く。基準大気は [Jablonowski–Williamson の初期状態](./Jablonowski-Williamson.md)の水平平均温度を $\widetilde T_k$ とし、地表面気圧は $p_s=p_0=10^5\,\mathrm{Pa}$、$\zeta=\delta=0$ とする。基準大気の量にはチルダをつけて表し、基準大気からのずれはプライムをつけて表すと、線形化した方程式は

$$
\begin{aligned}
\frac{\partial(\ln p_s)'}{\partial t}&=-\frac{1}{p_0}\sum_{j=1}^N\widetilde{\Delta p_j}\delta'_j\\
\frac{\partial\delta'_k}{\partial t}&=-\nabla_\eta^2\left[R\widetilde\alpha_kT'_k+R\sum_{j=k+1}^N\widetilde L_jT'_j+RP_k(\ln p_s)'\right]\\
\frac{\partial T'_k}{\partial t}&=-\frac{1}{2\widetilde{\Delta p_k}}\left[M'_{k+1/2}(\widetilde T_{k+1}-\widetilde T_k)+M'_{k-1/2}(\widetilde T_k-\widetilde T_{k-1})\right]+\kappa\widetilde T_k\left(-\widetilde\alpha_k\delta'_k-\frac{\widetilde L_k}{\widetilde{\Delta p_k}}S'_{k-1}\right)
\end{aligned}
$$

となる。ただし

$$
\begin{aligned}
\Lambda_k&=p_0\left(\frac{B_{k+1/2}}{\widetilde p_{k+1/2}}-\frac{B_{k-1/2}}{\widetilde p_{k-1/2}}\right),\\
\Gamma_k&=-\frac{p_0B_{k-1/2}}{\widetilde{\Delta p_k}}\widetilde L_k+\frac{\widetilde p_{k-1/2}p_0\Delta B_k}{\widetilde{\Delta p_k}^2}\widetilde{L}_k-\frac{\widetilde p_{k-1/2}}{\widetilde{\Delta p_k}}\Lambda_k,\\
C_k&=\frac{p_0}{\widetilde{\Delta p_k}}(B_{k-1/2}\widetilde L_k+\widetilde\alpha_k\Delta B_k)
\end{aligned}
$$

は事前に計算できる量で、

$$
\begin{aligned}
S'_k&=\sum_{j=1}^{k}\widetilde{\Delta p_j}\delta'_j,\\
M'_{k+1/2}&=B_{k+1/2}S'_N-S'_k,\\
P_k&=\widetilde T_k\Gamma_k+\sum_{j=k+1}^N\widetilde T_j\Lambda_j+\widetilde T_kC_k
\end{aligned}
$$

とした。この重力波の方程式の右辺を $GX$ と置き、元の方程式の右辺を $\mathcal{R}(X)$ と置いたとき、陰的にする度合いを $\beta=0.5$ として

$$
\frac{X_*^{n+1}-\overline{X}^{n-1}}{2\Delta t}=\mathcal{R}(X^n)-GX^n+(1-\beta)G\overline{X}^{n-1}+\beta GX^{n+1}_*
$$

とする。つまり

$$
(1-2\beta\Delta tG)X^{k+1}_*=\overline{X}^{k-1}+2\Delta t[\mathcal{R}(X^k)-G(X^k)+(1-\beta)G\overline{X}^{k-1}]
$$

として $X^{n+1}_*$ を得る。時間ステップは $\Delta t=900\,\mathrm{s}$ とする。

## 超粘性

超粘性は相対渦度 $\zeta_k$、発散 $\delta_k$、温度 $T_k$ にだけ入れ、地表面気圧 $\ln p_s$ には入れない。球面調和モード $(n,m)$ に対し、変数 $Y\in\{\zeta,\delta,T\}$ の減衰率を

$$
\kappa_n^{(Y)}=\frac{1}{\tau_Y}\left(\frac{n(n+1)}{n_\mathrm{max}(n_\mathrm{max}+1)}\right)^p
$$

とする。ここで $n_\mathrm{max}$ は切断波数、$p=4$ とし、時定数は

$$
\tau_\delta=1\,\mathrm{hour},\qquad
\tau_\zeta=\tau_T=4\,\mathrm{hours}
$$

とする。重力波を陰的に解いて得た $X_*^{q+1}$ に対し、各鉛直層の各スペクトル係数を

$$
Y_{**}^{q+1}=\frac{Y_*^{q+1}}{1+2\Delta t\,\kappa_n^{(Y)}},\qquad
Y\in\{\zeta,\delta,T\}
$$

と更新する。一方、地表面気圧は

$$
(\ln p_s)_{**}^{q+1}=(\ln p_s)_*^{q+1}
$$

とする。

## RAW フィルター

全予報変数をまとめて

$$
X=(\zeta_1,\dots,\zeta_N,\delta_1,\dots,\delta_N,T_1,\dots,T_N,\ln p_s)^\mathsf{T}
$$

と書く。RAW フィルターは浅水方程式と同様に

$$
d^q=\frac{\epsilon}{2}(\overline X^{q-1}-2X^q+X_{**}^{q+1})
$$

と定義し、

$$
\begin{aligned}
\overline X^q&=X^q+\alpha d^q,\\
X^{q+1}&=X_{**}^{q+1}-(1-\alpha)d^q
\end{aligned}
$$

とする。ここでは

$$
\epsilon=0.1,\qquad\alpha=0.53
$$

と置く。この更新後、各鉛直層で $\zeta_0^0=\delta_0^0=0$ とする。

## 初期化

重力波の陰的更新、超粘性、RAW フィルターをまとめて

$$
X^{q+1}=F(\Delta t,\overline X^{q-1},X^q)
$$

と書く。初期値 $X^0$ には [Jablonowski--Williamson の初期状態](./Jablonowski-Williamson.md)を用いる。最初の二時刻は

$$
\begin{aligned}
\overline X^0&=X^0,\\
X^{1/2}&=F_\mathrm{noRAW}\left(\frac{\Delta t}{4},X^0,X^0\right),\\
X^1&=F\left(\frac{\Delta t}{2},X^0,X^{1/2}\right)
\end{aligned}
$$

として作り、LeapFrog を始める。ただし $X^{1/2}$ の計算にだけは RAW フィルターを使わない。
