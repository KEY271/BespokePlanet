# 湿潤対流調節

条件付き不安定な気柱を、Frierson (2007) の簡略化 Betts–Miller 型の対流調節で緩和する。最下層から持ち上げた空気塊の湿潤断熱線を参照温度とし、参照温度に対する一定の相対湿度を参照比湿として、温度と比湿を時定数 $\tau_c$ でそこへ近づける。空気塊に正の浮力がある（CAPE が正の）気柱だけを対象とし、参照プロファイルは気柱の湿潤エンタルピー $c_pT+Lq$ を保存するように平行移動し、除かれた水蒸気を対流性降水とする。雲は扱わず、降水は直ちに地表へ落ちる。定数、飽和比湿、相当温位、持ち上げ凝結高度は[飽和比湿](./saturation-specific-humidity.md)のとおりである。

鉛直層の添字は[乾燥大気](../dynamics/dry.md)と同じく上端から地表へ向かって $k=1,2,\dots,N$ とし、full level の気圧は $p_k=p_{k+1/2}\exp(-\alpha_k)$ とする。

## 参照プロファイル

### 空気塊の持ち上げ

最下層 $k=N$ の温度と比湿 $T_N,q_N$ から空気塊を持ち上げる。$k=N-1,N-2,\dots,1$ の順に、まず乾燥断熱で持ち上げた温度

$$
T^{\mathrm{dry}}_k=T_N\left(\frac{p_k}{p_N}\right)^\kappa
$$

を計算する。$q_N\le q_s(T^{\mathrm{dry}}_k,p_k)$ なら空気塊はまだ飽和していないので

$$
T_{\mathrm{ref},k}=T^{\mathrm{dry}}_k,\qquad q_{p,k}=q_N
$$

とする。ここで $q_{p,k}$ は空気塊の比湿である。初めて $q_N>q_s(T^{\mathrm{dry}}_k,p_k)$ となった層 $k_L$ で、[飽和比湿](./saturation-specific-humidity.md#持ち上げ凝結高度)の手順により持ち上げ凝結高度 $T_L,p_L$ と空気塊の相当温位

$$
\theta_e^L=\theta_e(T_L,q_N,p_L)
$$

を求める。$p_k<p_L$ となる $k\le k_L$ の層では空気塊は飽和しており、飽和した空気塊の相当温位が保存する条件

$$
\theta_e\left(T_{\mathrm{ref},k},\,q_s(T_{\mathrm{ref},k},p_k),\,p_k\right)=\theta_e^L,\qquad
q_{p,k}=q_s(T_{\mathrm{ref},k},p_k)
$$

から $T_{\mathrm{ref},k}$ を決める。乾燥断熱線と湿潤断熱線の切り替えは層単位で行い、持ち上げ凝結高度を層の途中に補間することはしない。$\theta_e^L$ を地表の $\theta_e(T_N,q_N,p_N)$ ではなく持ち上げ凝結高度で評価するのは、簡略化した $\theta_e$ の式が未飽和の乾燥断熱上昇の間は一定でないためである。

飽和した層の $T_{\mathrm{ref},k}$ は、

$$
f(T)=\ln\theta_e(T,q_s(T,p_k),p_k)-\ln\theta_e^L
$$

の根として求める。$f$ は $T$ について単調増加で、根は区間

$$
T^{\mathrm{dry}}_k\le T_{\mathrm{ref},k}\le T_{\mathrm{ref},k+1}
$$

にある。下限は、乾燥断熱線に沿って $q_s/T$ が温度とともに増えることから $f(T^{\mathrm{dry}}_k)\le0$ となるためである。上限は、ひとつ下の層の参照温度をそのまま気圧 $p_k$ に置くと温位と飽和比湿の両方が大きくなり $f(T_{\mathrm{ref},k+1})\ge0$ となるためで、$k+1$ が未飽和の層なら $T_{\mathrm{ref},k+1}=T^{\mathrm{dry}}_{k+1}$（$k+1=N$ なら $T_N$）が同じ乾燥断熱線上にあることから従う。区間の上端 $T_{\mathrm{ref},k+1}$ を初期値として Newton 法

$$
T\leftarrow T-\frac{f(T)}{\partial\ln\theta_e/\partial T}
$$

を反復する。各反復で $f$ の符号により区間を狭め、Newton の更新が区間の外に出たときは区間の中点に置き換える。$|\Delta T|<10^{-3}\,\mathrm{K}$ または区間の幅が $10^{-3}\,\mathrm{K}$ 未満になったら収束とし、最大 20 回反復する。20 回で収束しなければ区間の中点を採用する。$\ln\theta_e$ は飽和比湿が大きい下層では $T$ について下に凸で、その場合 Newton 法は区間の上端から単調に収束するが、飽和比湿が小さい上層では凸性が保証されないので区間による保護を残す。

参照比湿は参照温度に対する相対湿度を一定値 $\mathrm{RH}_c$ として

$$
q_{\mathrm{ref},k}=\mathrm{RH}_cq_s(T_{\mathrm{ref},k},p_k),\qquad\mathrm{RH}_c=0.7
$$

とする。最下層でも同じ式を使い、$T_{\mathrm{ref},N}=T_N$、$q_{p,N}=q_N$ である。$q_N\le0$ のとき空気塊はどの層でも飽和せず、参照温度は全層で乾燥断熱線になる。

### 対流の深さ

空気塊の浮力は仮想温度で判定する。空気塊と環境の仮想温度を

$$
T_{v,\mathrm{ref},k}=(1+\delta_vq_{p,k})T_{\mathrm{ref},k},\qquad
T_{v,k}=(1+\delta_vq_k^+)T_k
$$

とし、$T_{v,\mathrm{ref},k}>T_{v,k}$ の層を浮力が正の層とする。[乾燥対流調節](#仮想温位による安定度)が仮想温位で安定度を判定するのと合わせてあり、$q=0$ のときは温度による判定に一致する。

$k=N-1,N-2,\dots,1$ の順に調べ、初めて浮力が正になった層を自由対流高度 $k_f$ とする。浮力が正の層がひとつもないときは CAPE が 0 以下であり、この気柱では湿潤対流調節を行わない。$k_f$ から上へ、浮力が正である間 $k$ を進め、初めて浮力が正でなくなる層の 1 つ下、すなわち $k_f$ から連続する浮力が正の層のうち最も上の層を浮力がなくなる高さ（対流頂）$k_t$ とする。

$$
k_t=\min\{k\le k_f:\;T_{v,\mathrm{ref},j}>T_{v,j}\ \text{for all}\ k\le j\le k_f\}
$$

$k_t$ より上に再び浮力が正になる層があっても、そこは対流層に含めない。一方 $k_f$ と $N$ の間に浮力が負の層があっても、$k_t\le k\le N$ の全層を対流層とする。この決め方は Frierson (2007) の level of zero buoyancy と同じであり、CAPE をこの連続した層で

$$
\mathrm{CAPE}=\sum_{k=k_t}^{k_f}R\left(T_{v,\mathrm{ref},k}-T_{v,k}\right)\ln\frac{p_{k+1/2}}{p_{k-1/2}}
$$

と定義すれば、湿潤対流調節を行う条件は $\mathrm{CAPE}>0$ と同値である。乾燥静的に不安定な成層は[乾燥対流調節](./dry-convective-adjustment.md)が扱うので、湿潤対流調節はそれを前提に条件付き不安定だけを扱う。

## エンタルピーの保存と降水

対流層 $k_t\le k\le N$ で、参照プロファイルへ緩和したときの温度と比湿の変化を

$$
\Delta T_k=T_{\mathrm{ref},k}-T_k,\qquad
\Delta q_k=q_{\mathrm{ref},k}-q_k
$$

と置く。気柱の質量を重みにした和を

$$
\Sigma_p=\sum_{k=k_t}^N\Delta p_k,\qquad
\Sigma_T=\sum_{k=k_t}^N\Delta T_k\Delta p_k,\qquad
\Sigma_q=\sum_{k=k_t}^N\Delta q_k\Delta p_k
$$

とする。比湿の変化だけから決まる降水と、温度の変化だけから決まる降水をそれぞれ

$$
P_q=-\frac{\Sigma_q}{g\tau_c},\qquad
P_T=\frac{c_p}{L}\frac{\Sigma_T}{g\tau_c}
$$

とする。$P_q>0$ は参照プロファイルが気柱を乾かすこと、$P_T>0$ は参照プロファイルが気柱を暖めることを意味する。一般には $P_q\ne P_T$ なので、次の 3 つに場合分けして参照プロファイルを補正し、両者を一致させる。

- $P_T\le0$：参照温度が平均として環境より冷たい。凝結の潜熱を放出しながら気柱を冷やすことはできないので、$P_q$ の符号にかかわらずこの気柱では湿潤対流調節を行わず、傾向をすべて 0 とする。
- $P_T>0$ かつ $P_q>0$：深い対流。
- $P_T>0$ かつ $P_q\le0$：浅い対流。

この場合分けは Frierson (2007) および Isca の実装と同じである。

### 深い対流

$P_T>0$ かつ $P_q>0$ のときは、参照温度を一定値 $\Delta T_{\mathrm{s}}$ だけ平行移動し、参照比湿は変えない。

$$
T_{\mathrm{ref},k}\leftarrow T_{\mathrm{ref},k}+\Delta T_{\mathrm{s}},\qquad
\Delta T_{\mathrm{s}}=\frac{1}{\Sigma_p}\left(-\frac{L}{c_p}\Sigma_q-\Sigma_T\right)
$$

これにより移動後の $\Sigma_T$ は $-(L/c_p)\Sigma_q$ に等しくなり、湿潤エンタルピーの気柱積分

$$
\sum_{k=k_t}^N\left(c_p\Delta T_k+L\Delta q_k\right)\Delta p_k=0
$$

が成り立つ。対流性降水は

$$
P_{\mathrm{conv}}=P_q>0
$$

である。

### 浅い対流

$P_T>0$ かつ $P_q\le0$、すなわち $\Sigma_q\ge0$ のときは、気柱全体としては参照比湿のほうが湿っているので、降水を伴わない浅い対流として扱う。対流の深さを、比湿の変化の気柱積分がちょうど 0 になるところまで減らし、その上でエンタルピーの気柱積分が変わらないように参照温度を平行移動する。まず地表から上へ向かう比湿の変化の累積

$$
s_k=\sum_{j=k}^N\Delta q_j\Delta p_j,\qquad s_{N+1}=0
$$

を $k=N,N-1,\dots,k_t$ の順に調べ、初めて $s_k\ge0$ になった層を $k^*$ とする。$s_{k_t}=\Sigma_q\ge0$ なので $k^*$ は必ず見つかる。$k^*=N$ のとき、すなわち最下層が参照比湿より乾いているときは対流層がなくなるので、この気柱では調節を行わない。

$k^*<N$ のときは新しい対流層を $k^*\le k\le N$ とし、$s_{k^*+1}<0\le s_{k^*}$ から $\Delta q_{k^*}\Delta p_{k^*}\ge-s_{k^*+1}>0$ であることを使って、符号が変わる層 $k^*$ の変化量だけを

$$
\Delta q_{k^*}\leftarrow c\,\Delta q_{k^*},\qquad
\Delta T_{k^*}\leftarrow c\,\Delta T_{k^*},\qquad
c=\frac{-s_{k^*+1}}{\Delta q_{k^*}\Delta p_{k^*}}\in(0,1]
$$

と縮める。これは対流頂を層 $k^*$ の途中に置くことに相当し、Isca の shallower 版と同じ扱いである。縮めた後は $\sum_{k=k^*}^N\Delta q_k\Delta p_k=0$ となる。参照比湿 $q_{\mathrm{ref},k}=q_k+\Delta q_k$ は、$k^*$ では $q_{k^*}$ と $\mathrm{RH}_cq_s(T_{\mathrm{ref},k^*},p_{k^*})$ の間、それより下の層では $\mathrm{RH}_cq_s(T_{\mathrm{ref},k},p_k)$ のままなので、参照比湿を一様に持ち上げる方法と違って参照温度に対して過飽和な参照比湿を作らず、降水を伴わないはずの浅い対流が[大規模凝結](./large-scale-condensation.md)の降水に変わることがない。

次に $k^*\le k\le N$ で

$$
\Sigma_p'=\sum_{k=k^*}^N\Delta p_k,\qquad
\Sigma_T'=\sum_{k=k^*}^N\Delta T_k\Delta p_k
$$

を計算し直し、温度の変化を一定値だけ平行移動する。

$$
\Delta T_k\leftarrow\Delta T_k+\Delta T_{\mathrm{s}},\qquad
\Delta T_{\mathrm{s}}=-\frac{\Sigma_T'}{\Sigma_p'}\qquad(k^*\le k\le N)
$$

移動後は比湿とエンタルピーの気柱積分の変化がどちらも 0 となり、比湿とエンタルピーは気柱内で鉛直に再配分されるだけである。対流性降水は

$$
P_{\mathrm{conv}}=0
$$

とする。

## 傾向項

補正後の変化量 $\Delta T_k,\Delta q_k$ を使い、対流層の各層で

$$
C^{\mathrm{conv}}_{T,k}=\frac{\Delta T_k}{\tau_c},\qquad
C^{\mathrm{conv}}_{q,k}=\frac{\Delta q_k}{\tau_c},\qquad
\tau_c=2\,\mathrm{hours}
$$

とする。対流層の外、および湿潤対流調節を行わない気柱では $C^{\mathrm{conv}}_{T,k}=C^{\mathrm{conv}}_{q,k}=0$ である。構成により

$$
\sum_k\left(c_pC^{\mathrm{conv}}_{T,k}+LC^{\mathrm{conv}}_{q,k}\right)\frac{\Delta p_k}{g}=0,\qquad
P_{\mathrm{conv}}=-\sum_kC^{\mathrm{conv}}_{q,k}\frac{\Delta p_k}{g}
$$

が成り立つ。運動量は混合しない。

## 時間積分への組み込み

湿潤対流調節は[湿潤大気](../dynamics/moist.md#物理過程)の順序のとおり、$\overline X^{n-1}$ に乾燥対流調節の傾向を $2\Delta t$ かけて加えた暫定場 $T^{(1)}_k,q^{(1)}_k$ から評価する。暫定場の $q^{(1)}_k$ が負の点では 0 とみなす。$p_k$ と $\Delta p_k$ は $(\ln p_s)^{\overline{n-1}}$ から作る。得られた $C^{\mathrm{conv}}_{T,k},C^{\mathrm{conv}}_{q,k}$ は温度と比湿の傾向の格子値に足し、他の項とまとめて 1 回だけスペクトルへ変換する。重力波の陰的演算子は変更しない。最初の 2 ステップでは、どちらの段階でも $X^0$ から作った暫定場で評価する。

$X^n$ で評価しない理由は、[乾燥対流調節](./dry-convective-adjustment.md)と同じく、緩和項を LeapFrog の中心差分で扱うと計算モードが増幅するためである。

## 乾燥対流調節の変更

湿潤大気では[乾燥対流調節](./dry-convective-adjustment.md)を次の 2 点だけ変更する。$q=0$ のときは元の手順に戻る。

### 仮想温位による安定度

安定度の判定には温位 $\theta_k$ の代わりに仮想温位

$$
\theta_{v,k}=\frac{T_{v,k}}{\Pi_k}=(1+\delta_vq_k^+)\theta_k
$$

を使い、隣り合う 2 層 $k,k+1$ について $\theta_{v,k}<\theta_{v,k+1}$ のとき乾燥静的に不安定とする。混合層を決める pool adjacent violators の各ブロックには $W,H$ に加えて

$$
Q=\sum_{k=k_t}^{k_b}q_k^+\Delta p_k,\qquad
\Delta P=\sum_{k=k_t}^{k_b}\Delta p_k
$$

を持たせ、ブロックの比湿と仮想温位を

$$
q_c=\frac{Q}{\Delta P},\qquad
\theta_{v,c}=(1+\delta_vq_c)\frac{H}{W}
$$

とする。併合の判定 $\theta^{\mathrm{upper}}_{v,c}<\theta^{\mathrm{lower}}_{v,c}$ にはこの $\theta_{v,c}$ を使い、併合したブロックの $Q,\Delta P$ はそれぞれの和とする。併合により $\theta_{v,c}$ は 2 つのブロックの値の間にあるとは限らないが、違反がある限り併合を続けるので、手順は有限回で終わり、残ったブロックの $\theta_{v,c}$ は上から下へ単調非増加になる。

### 比湿の混合

参照温度 $T_{\mathrm{ref},k}=\Pi_k\theta_c$、$\theta_c=H/W$ は元のままとし、これに加えて参照比湿をブロック内で一様に

$$
q_{\mathrm{ref},k}=q_c
$$

とする。比湿の傾向

$$
C^{\mathrm{dry}}_{q,k}=-\frac{q_k^+-q_{\mathrm{ref},k}}{\tau},\qquad\tau=4\,\mathrm{hours}
$$

を温度の傾向 $C^{\mathrm{dry}}_{T,k}$ と同じ $\tau$ で加える。構成により各ブロックで $\sum_k(q_k^+-q_{\mathrm{ref},k})\Delta p_k=0$ なので、乾燥対流調節は気柱の水蒸気量を変えず、降水も潜熱の放出も伴わない。混合によって上部の層が飽和を超えたときは、その過飽和を同じステップの[大規模凝結](./large-scale-condensation.md)が除く。
