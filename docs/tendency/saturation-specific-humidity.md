# 飽和比湿

水蒸気に関わる定数と、[蒸発](./evaporation.md)、[大規模凝結](./large-scale-condensation.md)、[湿潤対流調節](./moist-convective-adjustment.md)で共通に使う飽和比湿とその温度微分、相当温位の定義をまとめる。鉛直層の添字と界面気圧は[長波放射](./longwave-radiation.md)と同じとし、full level の気圧は[乾燥対流調節](./dry-convective-adjustment.md)と同じく

$$
p_k=p_{k+1/2}\exp(-\alpha_k)
$$

とする。

## 定数

| 記号 | 値 | 意味 |
| --- | --- | --- |
| $R$ | $287\,\mathrm{J\,kg^{-1}\,K^{-1}}$ | 乾燥空気の気体定数 |
| $R_v$ | $461.5\,\mathrm{J\,kg^{-1}\,K^{-1}}$ | 水蒸気の気体定数 |
| $\varepsilon$ | $R/R_v\approx0.622$ | 気体定数の比 |
| $\delta_v$ | $1/\varepsilon-1\approx0.608$ | 仮想温度の係数 |
| $c_p$ | $1004\,\mathrm{J\,kg^{-1}\,K^{-1}}$ | 乾燥空気の定圧比熱 |
| $L$ | $2.5\times10^6\,\mathrm{J\,kg^{-1}}$ | 凝結の潜熱 |
| $e_0$ | $610.78\,\mathrm{Pa}$ | $T_0$ における飽和水蒸気圧 |
| $T_0$ | $273.16\,\mathrm{K}$ | 基準温度 |

$L$ は温度によらない定数とし、氷相は扱わないので昇華の潜熱は使わない。$c_p$ と $\kappa=2/7$ は水蒸気を含んでも乾燥空気の値のままとする。

## 飽和水蒸気圧

飽和水蒸気圧は、$L$ を一定とした Clausius–Clapeyron の式を積分した

$$
e_s(T)=e_0\exp\left[\frac{L}{R_v}\left(\frac{1}{T_0}-\frac{1}{T}\right)\right]
$$

で与える。$L/R_v\approx5417\,\mathrm{K}$ である。

## 飽和比湿

気圧 $p$、温度 $T$ における飽和比湿は

$$
q_s(T,p)=\frac{\varepsilon e_s(T)}{p-(1-\varepsilon)e_s(T)}
$$

とする。水蒸気が希薄な範囲では $q_s\approx\varepsilon e_s/p$ であるが、分母の補正は計算量に影響しないので残す。$p\le e_s(T)$ となる高温では分母が意味を持たないので $q_s=1$ とし、そこでは次の温度微分も $\partial q_s/\partial T=0$ とする。$q\le1$ なのでこの層は常に未飽和として扱われ、凝結は起こらない。

飽和比湿の温度微分は

$$
\frac{\partial q_s}{\partial T}
=q_s\frac{L}{R_vT^2}\frac{p}{p-(1-\varepsilon)e_s}
$$

である。これは凝結量を Newton 法で求めるときに使う。$e_s$ は $T$ について下に凸で、$q_s$ は $e_s$ の単調増加かつ下に凸な関数なので、$q_s$ も $T$ について下に凸である。

第 $k$ 層の飽和比湿は $q_{s,k}=q_s(T_k,p_k)$、地表面の飽和比湿は地表面温度 $T_s$（slab ocean では $T_o$）と地表面気圧から $q_s(T_s,p_s)$ とする。相対湿度は $q_k/q_{s,k}$ で定義する。

## 相当温位と持ち上げ凝結高度

湿潤対流調節で湿潤断熱線を求めるために、相当温位を

$$
\theta_e(T,q,p)=T\left(\frac{p_0}{p}\right)^\kappa\exp\left(\frac{Lq}{c_pT}\right),\qquad p_0=10^5\,\mathrm{Pa}
$$

と定義する。飽和した空気塊では $q=q_s(T,p)$ を代入する。この式は、飽和した空気塊が擬断熱的に上昇するときの熱力学第一法則 $c_p\,d\ln T-R\,d\ln p=-L\,dq_s/T$ を、$L\,dq_s/T\approx d(Lq_s/T)$ と近似して積分したものであり、$(Lq_s/T)(dT/T)$ の項を落としている。この近似の範囲で、飽和した空気塊の $\theta_e(T,q_s(T,p),p)$ は擬断熱上昇の間一定とみなす。一方、未飽和の空気塊が乾燥断熱的に上昇するときは $q$ が一定でも $T$ が下がるので、この式の値は一定ではなく、地表から持ち上げ凝結高度までの上昇で 1 K 程度変わる。そのため湿潤断熱線の起点にする $\theta_e$ は地表の値ではなく、次に述べる持ち上げ凝結高度で評価する。

飽和した空気塊について $\ln\theta_e$ の温度微分は

$$
\frac{\partial\ln\theta_e}{\partial T}
=\frac{1}{T}+\frac{L}{c_pT}\left(\frac{\partial q_s}{\partial T}-\frac{q_s}{T}\right)
$$

であり、$L/(R_vT)\gg1$ なので常に正である。したがって $p$ を固定したとき $\theta_e(T,q_s(T,p),p)$ は $T$ について単調増加であり、与えられた $\theta_e$ に対する飽和温度は一意に定まる。

### 持ち上げ凝結高度

温度 $T_N$、比湿 $q_N>0$ の空気塊を気圧 $p_N$ から乾燥断熱的に持ち上げたとき、ちょうど飽和する気圧 $p_L$ と温度 $T_L$ を持ち上げ凝結高度とする。すなわち

$$
q_s(T_L,p_L)=q_N,\qquad p_L=p_N\left(\frac{T_L}{T_N}\right)^{1/\kappa}
$$

を満たす $T_L$ を求める。$q_N\ge q_s(T_N,p_N)$ なら空気塊は最初から飽和しているので $T_L=T_N$、$p_L=p_N$ とする。そうでなければ

$$
g(T)=\ln q_s\!\left(T,\,p_N\left(\frac{T}{T_N}\right)^{1/\kappa}\right)-\ln q_N
$$

を初期値 $T_N$ の Newton 法で解く。乾燥断熱線に沿った微分は

$$
\frac{dg}{dT}=\frac{p}{p-(1-\varepsilon)e_s}\,\frac{1}{T}\left(\frac{L}{R_vT}-\frac{1}{\kappa}\right)
$$

で、$L/(R_vT)\approx18\gg1/\kappa=3.5$ なので正である。$g$ は $T$ について単調増加かつこの温度範囲では上に凸なので、$g(T_N)>0$ の初期値から始めた Newton 法は 1 回目で根を下側に越えたあと単調に収束する。$|\Delta T|<10^{-3}\,\mathrm{K}$ で収束とし、最大 20 回反復する。持ち上げた空気塊の相当温位は

$$
\theta_e^L=\theta_e(T_L,q_N,p_L)
$$

とし、これを持ち上げ凝結高度より上で飽和した空気塊が保存する量とする。
