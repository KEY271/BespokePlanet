# 湿潤大気

[乾燥大気](./dry.md)に比湿 $q$ を予報変数として加え、水蒸気が力学に及ぼす影響を仮想温度（仮温度）$T_v$ を通して取り入れる。水蒸気は地球大気と同程度に希薄（$q\ll1$）であると近似し、次の扱いをする。

- 水蒸気の質量が地表面気圧に与える影響は無視する。すなわち $p_s$ は乾燥空気の質量で決まるものとみなし、蒸発と降水による $p_s$ の変化は考えない。連続の式と $\ln p_s$ の方程式は変更しない。
- 水蒸気による密度の変化は $T_v$ を通して静水圧の関係、圧力傾度力、断熱昇温項に取り入れる。一方、水蒸気による比熱の変化は無視し、$c_p$ と $\kappa$ は乾燥空気の値のままとする。
- 雲は扱わない。凝結した水は大気中に滞留せず、直ちに降水として地表へ落ちる。落下中の再蒸発や、降水が運ぶエンタルピーも無視する。
- 水の相は水蒸気と液体の水だけとし、氷は扱わない。凝結の潜熱 $L$ は一定とする。

水蒸気に関わる定数と飽和比湿の式は[飽和比湿](../tendency/saturation-specific-humidity.md)にまとめる。

## 予報変数と基礎方程式

各ステップで時間発展させる変数は、相対渦度 $\zeta$ と発散 $\delta$、温度 $T$、比湿 $q$、地表面気圧 $\ln p_s$ である。比湿 $q$ は湿潤空気の単位質量に含まれる水蒸気の質量で、無次元である。地表ジオポテンシャル $\Phi_s$ は乾燥大気と同じく固定する。記号は [hybrid-sigma 座標](./hybrid-sigma.md)と[乾燥大気](./dry.md)に従う。

### 仮想温度

乾燥空気の気体定数を $R$、水蒸気の気体定数を $R_v$ とすると、比湿 $q$ の湿潤空気の気体定数は

$$
R_m=(1-q)R+qR_v=R(1+\delta_vq),\qquad
\delta_v=\frac{R_v}{R}-1=\frac{1}{\varepsilon}-1
$$

である。ここで $\varepsilon=R/R_v$ である。状態方程式 $p=\rho R_mT$ を $p=\rho RT_v$ と書けるように、仮想温度を

$$
T_v=(1+\delta_vq)T
$$

と定義する。したがって乾燥大気の式のうち密度 $1/\rho=RT/p$ に由来する $RT$ を $RT_v$ に置き換えれば、湿潤空気の式になる。

### 方程式

基礎方程式は

$$
\begin{aligned}
\frac{\partial\ln p_s}{\partial t}&=-\frac{1}{p_s}\bm\nabla_\eta\cdot\left(\int_{p_T}^{p_s}\bm{u}dp\right)\\
\frac{\partial\zeta}{\partial t}&=-\bm\nabla_\eta\cdot[(\zeta+f)\bm{u}]+\bm{k}\cdot\bm\nabla_\eta\times\left[-\frac{M}{\mu}\frac{\partial\bm{u}}{\partial\eta}-RT_v\bm\nabla_\eta\ln p+\bm{F}\right]\\
\frac{\partial\delta}{\partial t}&=-\nabla_\eta^2(K+\Phi)+\bm{k}\cdot\bm\nabla_\eta\times[(\zeta+f)\bm{u}]+\bm\nabla_\eta\cdot\left[-\frac{M}{\mu}\frac{\partial\bm{u}}{\partial\eta}-RT_v\bm\nabla_\eta\ln p+\bm{F}\right]\\
\frac{\partial T}{\partial t}&=-\bm{u}\cdot\bm\nabla_\eta T-\frac{M}{\mu}\frac{\partial T}{\partial\eta}+\kappa T_v\frac{D\ln p}{Dt}+S_T\\
\frac{\partial q}{\partial t}&=-\bm{u}\cdot\bm\nabla_\eta q-\frac{M}{\mu}\frac{\partial q}{\partial\eta}+S_q
\end{aligned}
$$

である。右辺に出てくる $p(\eta)$、$\mu$、$M$ は乾燥大気と同じで、ジオポテンシャルだけが

$$
\Phi=\Phi_s+\int_{\ln p(\eta)}^{\ln p_s}RT_vd\ln p
$$

と変わる。

温度の式は、熱力学第一法則 $c_pDT/Dt=(1/\rho)Dp/Dt+\dot{Q}$ に $1/\rho=RT_v/p$ を代入して得たものである。$\kappa=R/c_p$ は乾燥空気の値を使う。$S_T$ は潜熱加熱や放射などの物理過程による加熱、$S_q$ は蒸発・凝結・対流による比湿の変化で、いずれも後述の物理過程 $\mathcal{P}$ に含める。水蒸気には保存則以外の力学的な生成・消滅はないので、$S_q=0$ のとき $q$ は物質微分が 0 の単なるトレーサーである。

## 鉛直方向の離散化

$\zeta_k,\delta_k,T_k,q_k$ を $k=1,2,\dots,N$ の full level に置く。圧力 $p_{k+1/2}$、$L_k$、$\Delta p_k$、$\alpha_k$、鉛直質量フラックス $M_{k+1/2}$、$F_k$、$S_k$、圧力傾度 $\bm{G}_k$、鉛直移流 $W_k(X)$、および $Q_k=(D\ln p/Dt)_k$ の離散化は[乾燥大気](./dry.md)と同じである。

### 仮想温度

各層の仮想温度は

$$
T_{v,k}=(1+\delta_vq_k^+)T_k,\qquad q_k^+=\max(q_k,0)
$$

とする。$q_k^+$ を使う理由は後述の[負の比湿](#負の比湿)のとおりである。

### ジオポテンシャル

ジオポテンシャルは $T_k$ の代わりに $T_{v,k}$ を使って

$$
\begin{aligned}
\Phi_{N+1/2}&=\Phi_s,\\
\Phi_{k-1/2}&=\Phi_{k+1/2}+RT_{v,k}L_k,\\
\Phi_k&=\Phi_{k+1/2}+\alpha_kRT_{v,k}
\end{aligned}
$$

と計算する。

### 比湿の移流

比湿の移流は温度と同じ形で

$$
-\bm{u}_k\cdot\bm\nabla_\eta q_k-W_k(q)
$$

と離散化する。$W_k(q)$ は乾燥大気の $W_k(X)$ に $X=q$ を代入したものである。

### 離散化された方程式

結局、離散化された方程式は

$$
\begin{aligned}
\frac{\partial\ln p_s}{\partial t}&=-\frac{S_N}{p_s},\\
\frac{\partial\zeta_k}{\partial t}&=-\bm\nabla_\eta\cdot[(\zeta_k+f)\bm{u}_k]+\bm{k}\cdot\bm\nabla_\eta\times[-W_k(\bm{u})-RT_{v,k}\bm{G}_k+\bm{F}_k],\\
\frac{\partial\delta_k}{\partial t}&=-\nabla_\eta^2(K_k+\Phi_k)+\bm{k}\cdot\bm\nabla_\eta\times[(\zeta_k+f)\bm{u}_k]+\bm\nabla_\eta\cdot[-W_k(\bm{u})-RT_{v,k}\bm{G}_k+\bm{F}_k],\\
\frac{\partial T_k}{\partial t}&=-\bm{u}_k\cdot\bm\nabla_\eta T_k-W_k(T)+\kappa T_{v,k}Q_k+S_{T,k},\\
\frac{\partial q_k}{\partial t}&=-\bm{u}_k\cdot\bm\nabla_\eta q_k-W_k(q)+S_{q,k}
\end{aligned}
$$

となる。乾燥大気の式との違いは、$RT_k$ が $RT_{v,k}$ に置き換わったこと、$\Phi_k$ を $T_{v,k}$ から計算すること、および $q_k$ の式が加わったことだけである。$q_k=0$ のとき乾燥大気の式に一致する。

## 水平方向の離散化

水平方向の離散化は[乾燥大気](./dry.md)と同じである。$q$ は $T$ と同じく球面調和関数のスペクトル係数として保持し、移流項の非線形積を評価するときに格子へ変換する。$\bm\nabla_\eta q_k$ はスペクトル係数から直接求める。$q_k$ は各層でスカラーなので、変換の回数は $T_k$ と同じだけ増える。

## 負の比湿

スペクトル法の切断誤差により、格子に変換した $q_k$ は水蒸気の少ない領域で負になりうる。負の値は、仮想温度 $T_{v,k}$ と物理過程（蒸発、凝結、対流調節、飽和の判定）の評価に用いる前に $q_k^+=\max(q_k,0)$ で置き換える。一方、移流項 $-\bm{u}_k\cdot\bm\nabla_\eta q_k-W_k(q)$ と超粘性、RAW フィルターには置き換えない値をそのまま使う。スペクトル係数自体は修正せず、負の値を隣接層から補う穴埋めのような質量保存の補正も行わない。したがって格子で見た比湿は局所的に負になりうる。

物理過程には $q^+$ を見せるが、そこから得た傾向はそのまま $q$ に加える。したがって蒸発 $E$ と降水 $P_{\mathrm{conv}}+P_{\mathrm{ls}}$ は、符号付きの比湿で定義した気柱水蒸気量 $\sum_kq_k\Delta p_k/g$ の変化として閉じており、$q^+$ で定義した可降水量の変化とは負の領域の分だけずれる。また移流を移流形 $-\bm{u}_k\cdot\bm\nabla_\eta q_k-W_k(q)$ で離散化しているので、超粘性は各層の全球平均 $q$ を変えないものの、移流の離散化誤差と $\Delta p_k$ の空間変化のために $q$ の質量重み付き全球積分は厳密には保存されない。正値性を保つ移流スキームや全球の水分補正（moisture fixer）は導入せず、その代わりに[湿潤ケース](../cases/moist.md#出力)の出力で符号付きの気柱水蒸気量と負の部分の量を別に診断し、全球の $dW/dt=E-P$ がどの程度閉じるかを長期積分で確認する。

## 重力波

重力波を陰的に解く基準大気は乾燥大気と同じ $\widetilde T_k=300\,\mathrm{K}$、$p_s=p_0$、$\zeta=\delta=0$ とし、比湿は $\widetilde q_k=0$ とする。このとき $\widetilde T_{v,k}=\widetilde T_k$ なので、線形化した重力波の方程式と演算子 $G$、およびそこから作る陰的に解く行列は[乾燥大気](./dry.md#重力波)から変わらない。仮想温度と温度の差 $T_{v,k}-T_k=\delta_vq_k^+T_k$ に由来する項はすべて陽解法部分 $\mathcal{D}(X)$ に含める。$q$ には重力波の項がないので、$q$ の更新は

$$
q_*^{n+1}=\overline q^{n-1}+2\Delta t\left[\mathcal{D}_q(X^n)+\mathcal{P}_q(\overline X^{n-1})\right]
$$

と陽的に行う。全予報変数をまとめた更新式

$$
(1-2\beta\Delta tG)X^{n+1}_*=\overline{X}^{n-1}+2\Delta t[\mathcal{D}(X^n)+\mathcal{P}(\overline X^{n-1})-G(X^n)+(1-\beta)G\overline{X}^{n-1}]
$$

の形と $\beta=0.5$、$\Delta t=1200\,\mathrm{s}$ は乾燥大気のままである。

## 超粘性

超粘性は $\zeta_k,\delta_k,T_k$ に加えて $q_k$ にも入れ、$\ln p_s$ には入れない。比湿の減衰率は温度と同じ

$$
\kappa_n^{(q)}=\frac{1}{\tau_q}\left(\frac{n(n+1)}{n_\mathrm{max}(n_\mathrm{max}+1)}\right)^p,\qquad
\tau_q=\tau_T=4\,\mathrm{hours},\quad p=4
$$

とし、

$$
q_{**}^{n+1}=\frac{q_*^{n+1}}{1+2\Delta t\,\kappa_n^{(q)}}
$$

と更新する。$\eta$ 面が傾いていることによる超粘性の補正は行わない。

## RAW フィルター

全予報変数を

$$
X=(\zeta_1,\dots,\zeta_N,\delta_1,\dots,\delta_N,T_1,\dots,T_N,q_1,\dots,q_N,\ln p_s)^\mathsf{T}
$$

とまとめ、乾燥大気と同じ $\epsilon=0.1$、$\alpha=0.53$ の RAW フィルターを $q$ にも適用する。更新後に各鉛直層で $\zeta_0^0=\delta_0^0=0$ とするのは乾燥大気と同じであり、$q$ の全球平均 $q_0^0$ には拘束を課さない。地面温度 $T_s,T_d$ や海洋温度 $T_o$ の扱いは[地面](../tendency/ground.md)と [slab ocean](../tendency/slab-ocean.md) のとおりである。

## 物理過程

物理過程 $\mathcal{P}$ は[物理過程を評価する時刻](../tendency/physics-time-level.md)のとおり、RAW フィルター適用済みの前時刻 $\overline X^{n-1}$ の場で評価する。湿潤大気では次の過程を含む。

- [長波放射](../tendency/longwave-radiation.md)、[短波放射](../tendency/shortwave-radiation.md)、[オゾン](../tendency/ozone.md)
- [地面](../tendency/ground.md)または [slab ocean](../tendency/slab-ocean.md) との顕熱交換、および[蒸発](../tendency/evaporation.md)
- [Held–Suarez 強制](../tendency/Held-Suarez.md)の形の地表摩擦と、ケースによっては[上層の Rayleigh 摩擦](../tendency/upper-rayleigh-friction.md)
- [乾燥対流調節](../tendency/dry-convective-adjustment.md)（[湿潤対流調節](../tendency/moist-convective-adjustment.md#乾燥対流調節の変更)の節で述べる変更を含む）
- [湿潤対流調節](../tendency/moist-convective-adjustment.md)
- [大規模凝結](../tendency/large-scale-condensation.md)

放射、顕熱、蒸発、摩擦、乾燥対流調節は互いに独立に $\overline X^{n-1}$ から評価し、傾向を足し合わせる。対流調節と凝結は同じ不安定・過飽和を重複して除かないように、次の順序で逐次に評価する。

1. $\overline X^{n-1}$ の格子値 $\overline T_k,\overline q_k^+,\overline{\bm u}_k,\overline p_s$ から乾燥対流調節の傾向 $C^{\mathrm{dry}}_{T,k},C^{\mathrm{dry}}_{q,k}$ を作る。
2. 暫定場 $T^{(1)}_k=\overline T_k+2\Delta tC^{\mathrm{dry}}_{T,k}$、$q^{(1)}_k=\overline q_k^++2\Delta tC^{\mathrm{dry}}_{q,k}$ から湿潤対流調節の傾向 $C^{\mathrm{conv}}_{T,k},C^{\mathrm{conv}}_{q,k}$ と対流性降水 $P_{\mathrm{conv}}$ を作る。
3. 暫定場 $T^{(2)}_k=T^{(1)}_k+2\Delta tC^{\mathrm{conv}}_{T,k}$、$q^{(2)}_k=q^{(1)}_k+2\Delta tC^{\mathrm{conv}}_{q,k}$ から大規模凝結の傾向 $C^{\mathrm{ls}}_{T,k},C^{\mathrm{ls}}_{q,k}$ と降水 $P_{\mathrm{ls}}$ を作る。

暫定場に使う $2\Delta t$ は LeapFrog で場を進める時間幅で、初期化の 2 ステップでは $F$ の第 1 引数を $\Delta t$ とみなして同じ規則を使う。暫定場には放射、顕熱、蒸発の傾向は加えない。暫定場を作るのは格子空間だけで済むので、スペクトルへの変換回数は増えない。温度と比湿の傾向はすべて格子で足し合わせ、他の項とまとめて 1 回だけスペクトルへ変換する。

$\mathcal{P}$ に含まれる各過程は、大気・地表面の間で交換する水とエネルギーについて次を満たす。蒸発で地表面が失う潜熱 $LE$ は同じ $E$ として大気最下層の比湿に入り、凝結・対流で大気から除かれた水蒸気はすべて降水 $P_{\mathrm{conv}}+P_{\mathrm{ls}}$ として地表へ落ち、その潜熱は同じ気柱の温度に加わる。したがって水と湿潤エンタルピー $c_pT+Lq$ の収支は、放射と顕熱を除いて気柱ごとに閉じる。

## 初期化

初期化の手順は[乾燥大気](./dry.md#初期化)と同じで、$X^0$ に $q^0$ を含める。$q^0$ はケースごとに与える。$\Phi_s$ は積分中ずっと固定する。
