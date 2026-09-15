# 乾燥対流調節

地表加熱や上層の放射冷却によって生じた乾燥静的に不安定な成層を緩和で中立へ戻す。すなわち各気柱で温位が鉛直に一定となる参照温度 $T_{\mathrm{ref},k}$ を作り、温度を時定数 $\tau$ でそこへ近づける：

$$
\left(\frac{\partial T_k}{\partial t}\right)_{\mathrm{conv}}=-\frac{T_k-T_{\mathrm{ref},k}}{\tau}.
$$

参照温度は気柱のエンタルピーを保存するように決める。瞬時に中立へ戻す hard adjustment とは違い、有限の $\tau$ で緩和する。水蒸気は扱わないので、潜熱も降水も考えない。運動量の鉛直混合も行わず、変更するのは温度 $T_k$ だけである。

鉛直層の添字は [乾燥大気](./dry.md) と同じく上端から地表へ向かって $k=1,2,\dots,N$ とする。

## 温位と安定度

各格子点で界面気圧 $p_{k+1/2}=A_{k+1/2}+B_{k+1/2}p_s$、層の厚さ $\Delta p_k$、$\alpha_k$ を [乾燥大気](./dry.md) と同様に計算する。full level の気圧は [Held–Suarez 強制](./Held-Suarez.md) と同じく

$$
p_k=p_{k+1/2}\exp(-\alpha_k)
$$

とする。$p_0=10^5\,\mathrm{Pa}$、$\kappa=2/7$ として Exner 関数と温位を

$$
\Pi_k=\left(\frac{p_k}{p_0}\right)^\kappa,\qquad
\theta_k=\frac{T_k}{\Pi_k}
$$

と置く。$k$ が大きいほど下の層なので、隣り合う 2 層 $k,k+1$ について

$$
\theta_k<\theta_{k+1}
$$

が成り立つとき、その 2 層は乾燥静的に不安定であるとする。$\theta_k=\theta_{k+1}$ は中立として扱い、調節しない。したがって安定な気柱では

$$
\theta_1\ge\theta_2\ge\dots\ge\theta_N
$$

が成り立つ。

## 参照プロファイル

### エンタルピー保存

乾燥空気の単位面積当たりのエンタルピーは $(c_p/g)\sum_kT_k\Delta p_k$ である。連続する層 $k_t\le k\le k_b$ をひとつの混合層とみなし、その中で温位を一定値 $\theta_c$ にそろえる。このとき混合の前後でエンタルピーが変わらない条件

$$
\sum_{k=k_t}^{k_b}\Pi_k\theta_c\Delta p_k=\sum_{k=k_t}^{k_b}T_k\Delta p_k
$$

から

$$
\theta_c=\frac{\sum_{k=k_t}^{k_b}T_k\Delta p_k}{\sum_{k=k_t}^{k_b}\Pi_k\Delta p_k}
$$

を得る。つまり $\theta_c$ は重み $w_k=\Pi_k\Delta p_k$ による $\theta_k$ の加重平均である。

### 混合層の決定

混合層は、隣り合うブロックの温位を比較して不安定なものを併合していく方法（pool adjacent violators）で決める。各ブロックは上端 $k_t$、下端 $k_b$ と

$$
W=\sum_{k=k_t}^{k_b}\Pi_k\Delta p_k,\qquad
H=\sum_{k=k_t}^{k_b}T_k\Delta p_k,\qquad
\theta_c=\frac{H}{W}
$$

を持つ。地表から上へ次の手順で処理する。

1. 空のスタックを用意する。
2. $k=N,N-1,\dots,1$ の順に、第 $k$ 層だけからなるブロック（$k_t=k_b=k$、$W=\Pi_k\Delta p_k$、$H=T_k\Delta p_k$）をスタックに積む。
3. 積むたびに、スタックに 2 個以上のブロックがあり、一番上のブロックの温位がそのすぐ下のブロックの温位より小さい間、すなわち $\theta_c^{\mathrm{upper}}<\theta_c^{\mathrm{lower}}$ である間は両者を併合する。併合したブロックは $k_t$ を上側から、$k_b$ を下側から取り、$W,H$ はそれぞれの和とする。
4. 全層を処理し終えたとき、スタックに残ったブロックが混合層である。

この手順で得られるブロックの温位は上から下へ単調非増加になり、安定条件を満たす。層を上から処理しても結果は変わらない。計算量は 1 気柱あたり $O(N)$ である。

参照温度は、第 $k$ 層が属するブロックの温位 $\theta_c$ を使って

$$
T_{\mathrm{ref},k}=\Pi_k\theta_c
$$

とする。1 層だけのブロックでは $T_{\mathrm{ref},k}=T_k$ となるので、不安定な層を含まない気柱は変化しない。構成により、各ブロックの中で

$$
\sum_{k=k_t}^{k_b}(T_k-T_{\mathrm{ref},k})\Delta p_k=0
$$

が成り立つ。したがって対流加熱の気柱積分は 0 である。

## 時間積分への組み込み

### 1 つ前の時刻で評価する傾向項

対流調節の項は、[乾燥大気](./dry.md) の温度の傾向項 $\mathcal{R}(X^n)$ に加える。ただし他の項とは異なり、現在の時刻 $X^n$ ではなく、RAW フィルター適用済みの 1 つ前の時刻 $\overline X^{n-1}$ で評価する。すなわち各格子点で $\overline T_k^{n-1}$ と $(\ln p_s)^{\overline{n-1}}$ から $\Pi_k$、$\Delta p_k$ と参照温度 $T_{\mathrm{ref},k}[\overline X^{n-1}]$ を作り、

$$
C_k=-\frac{\overline T_k^{n-1}-T_{\mathrm{ref},k}[\overline X^{n-1}]}{\tau}
$$

とする。緩和の時定数は

$$
\tau=4\,\mathrm{hours}
$$

とする。そして温度の傾向項を

$$
\mathcal{R}_{T_k}(X^n)\;\to\;\mathcal{R}_{T_k}(X^n)+C_k
$$

と置き換える。$\zeta,\delta,\ln p_s$ の傾向項と、重力波の線形演算子 $G$ は変更しない。$C_k$ は $\mathcal{R}$ に含まれるので、重力波の陰的解法、超粘性、RAW フィルターは [乾燥大気](./dry.md) の手順のまま適用する。地面温度 $T_s,T_d$ の傾向項には加えない。放射と Held–Suarez 強制の項は、従来どおり $X^n$ で評価する。

$X^n$ で評価しない理由は、減衰項を LeapFrog の中心差分で扱うと計算モードが増幅するためである。

### 格子とスペクトルの変換

$\overline T_k^{n-1}$（$k=1,\dots,N$）と $(\ln p_s)^{\overline{n-1}}$ は、傾向項を評価する各ステップでスペクトルから格子へ変換する。得られた $C_k$ は放射の加熱率と同じく温度の傾向の格子値に足し、他の項とまとめて 1 回だけスペクトルへ変換する。そのため、スペクトルへの変換回数は増えない。不安定な層を含まない格子点では $C_k=0$ となる。

### 初期化

最初の 2 ステップでは $F$ の第 2 引数を 1 つ前の時刻として用いる。すなわち

$$
\begin{aligned}
X^{1/2}&=F_\mathrm{noRAW}\left(\frac{\Delta t}{4},X^0,X^0\right),\\
X^1&=F\left(\frac{\Delta t}{2},X^0,X^{1/2}\right)
\end{aligned}
$$

のどちらでも $C_k$ は $X^0$ で評価する。
