# 長波放射

鉛直層の添字は上端から地表へ向かって $k=1,2,\dots,N$ とし、界面気圧を $p_{k+1/2}$、層の厚さを

$$
\Delta p_k=p_{k+1/2}-p_{k-1/2}
$$

とする。大気上端は $p_{1/2}=p_T$、地表は $p_{N+1/2}=p_s$ である。基準気圧は $p_0=10^5\,\mathrm{Pa}$ とする。

## 長波の光学的厚さ

長波放射は波長によらない灰色大気として、上向きと下向きの 2-stream で扱う。各大気層は長波を吸収し、Kirchhoff の法則に従って同じ放射率で放射する。第 $k$ 層の長波光学的厚さは、よく混合された吸収体による項、水蒸気による項、[オゾン](./ozone.md)による項の和

$$
\Delta\tau^{\mathrm{LW}}_k
=\Delta\tau^{\mathrm{CO_2}}_k+\Delta\tau^{\mathrm{H_2O}}_k+\Delta\tau^{\mathrm{O_3,LW}}_k
$$

とする。前の 2 項は Byrne & O'Gorman (2013) が Frierson et al. (2006) の灰色放射を水蒸気に依存するように拡張したもの（Isca の `two_stream_gray_rad` における `byrne` スキーム）に従い、光学的厚さの気圧微分を

$$
\frac{d\tau}{dp}=\frac{a\mu+bq}{p_0}
$$

と置く。すなわち

$$
\Delta\tau^{\mathrm{CO_2}}_k=a\mu\frac{\Delta p_k}{p_0},\qquad
\Delta\tau^{\mathrm{H_2O}}_k=bq_k\frac{\Delta p_k}{p_0}
$$

で、係数は

$$
a=0.1627,\qquad b=1997.9,\qquad \mu=1,\qquad p_0=10^5\,\mathrm{Pa}
$$

とする。$a$ は $\mathrm{CO_2}$ のように大気に一様に混合した吸収体を表し、単位質量当たりの吸収が一定なので光学的厚さは層の質量 $\Delta p_k/g$ に比例する。$b$ は水蒸気の単位質量当たりの吸収を表し、比湿 $q_k$（無次元、$\mathrm{kg\,kg^{-1}}$）に掛けることで、層に含まれる水蒸気の質量 $q_k\Delta p_k/g$ に比例する光学的厚さを与える。$a,b,\mu$ は Isca の `byrne` スキームの既定値（`bog_a`, `bog_b`, `bog_mu`）をそのまま使う。$\mu$ は一様吸収体の量を現在の地球に対する比で変えるための係数で、ここでは $1$ に固定する。分母は局所の $p_s$ ではなく定数 $p_0$ とし、$p_s$ が変わっても単位質量当たりの吸収が変わらないようにする。地球相当の設定では $a\mu$ の寄与は気柱全体で $0.16$ にすぎず、長波の光学的厚さの大部分は水蒸気の項が担う。

層の透過率と放射率は

$$
t_k=\exp\left(-\Delta\tau^{\mathrm{LW}}_k\right),\qquad
\varepsilon_k=1-t_k
$$

である。オゾンを含む全光学的厚さに対して Kirchhoff の法則を適用する。放射が鉛直に進むものとして拡散係数による光路長の補正は行わない。$a,b$ はこの形の 2-stream のフラックスに対して定めた値なので、拡散係数を別に掛けてはならない。

### 比湿の与え方

$q_k$ には、ケースによって次のいずれかを使う。

- 比湿を予報する[湿潤ケース](../cases/moist.md)では、[物理過程を評価する時刻](./physics-time-level.md)のとおり前時刻の比湿の正の部分 $\overline q_k^{+,n-1}=\max(\overline q_k^{n-1},0)$ を使う。負の比湿の扱いは[湿潤大気](../dynamics/moist.md#負の比湿)のとおりである。これにより、温度が上がって比湿が増えると長波の光学的厚さが増え、大気上端から出る長波が減る簡易的な水蒸気フィードバックが働く。飽和比湿は Clausius–Clapeyron の式に従い $1\,\mathrm{K}$ 当たり $L/(R_vT^2)\approx6.5\,\%$（$T=288\,\mathrm{K}$）増えるので、相対湿度がほぼ保たれる限り $\Delta\tau^{\mathrm{H_2O}}$ も同じ割合で増える。
- 比湿を予報しない[放射ケース](../cases/radiation.md)と [slab ocean ケース](../cases/slab-ocean.md)では、次の基準水蒸気分布 $\overline q^{\mathrm{ref}}_k$ を使う。

### 基準水蒸気分布

乾燥ケースでは、水蒸気の平均的な放射効果を表す固定の分布

$$
q^{\mathrm{ref}}(p)=q_0\left(\frac{p}{p_s}\right)^3,\qquad q_0=0.010
$$

を与える。第 $k$ 層に使う値は層の質量重み付き平均

$$
\overline q^{\mathrm{ref}}_k
=\frac{1}{\Delta p_k}\int_{p_{k-1/2}}^{p_{k+1/2}}q^{\mathrm{ref}}(p)\,dp
=\frac{q_0}{4p_s^3\Delta p_k}\left(p_{k+1/2}^4-p_{k-1/2}^4\right)
$$

とする。これを $\Delta\tau^{\mathrm{H_2O}}_k$ に代入すると

$$
\Delta\tau^{\mathrm{H_2O}}_k
=\frac{bq_0}{4}\frac{p_s}{p_0}\left[\left(\frac{p_{k+1/2}}{p_s}\right)^4-\left(\frac{p_{k-1/2}}{p_s}\right)^4\right]
$$

となり、Frierson et al. (2006) が水蒸気を模して置いた $(p/p_s)^4$ に比例する光学的厚さと同じ形になる。$q_0$ は、この分布の可降水量

$$
W^{\mathrm{ref}}=\int_0^{p_s}\frac{q^{\mathrm{ref}}}{g}dp=\frac{q_0p_s}{4g}\approx25.5\,\mathrm{kg\,m^{-2}}
$$

が地球の全球平均の可降水量（約 $25\,\mathrm{kg\,m^{-2}}$）に一致するように選ぶ。基準分布には緯度依存性を持たせない。Frierson et al. は赤道 $\tau_0=6$、極 $\tau_0=1.5$ の緯度依存性を与えているが、この緯度構造は湿潤ケースでは予報した比湿から自然に現れるので、乾燥ケースには入れない。

$p_s=p_0$ のとき、$N=12$ の各層の光学的厚さは次のとおりである。$\tau_{k+1/2}$ は上端から層下端までの累積で、オゾンは含まない。

| $k$ | 層 [hPa] | $\overline q^{\mathrm{ref}}_k$ | $\Delta\tau^{\mathrm{CO_2}}_k$ | $\Delta\tau^{\mathrm{H_2O}}_k$ | $\tau_{k+1/2}$ | $\varepsilon_k$ |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 1--3 | 0 | 0.000 | 0.000 | 0.000 | 0.000 |
| 2 | 3--10 | 0 | 0.001 | 0.000 | 0.001 | 0.001 |
| 3 | 10--50 | 0.000 00 | 0.007 | 0.000 | 0.008 | 0.007 |
| 4 | 50--100 | 0.000 00 | 0.008 | 0.000 | 0.017 | 0.009 |
| 5 | 100--180 | 0.000 03 | 0.013 | 0.005 | 0.034 | 0.018 |
| 6 | 180--280 | 0.000 13 | 0.016 | 0.025 | 0.076 | 0.041 |
| 7 | 280--400 | 0.000 41 | 0.020 | 0.097 | 0.193 | 0.110 |
| 8 | 400--520 | 0.000 99 | 0.020 | 0.237 | 0.450 | 0.227 |
| 9 | 520--650 | 0.002 03 | 0.021 | 0.526 | 0.997 | 0.422 |
| 10 | 650--770 | 0.003 60 | 0.020 | 0.864 | 1.881 | 0.587 |
| 11 | 770--880 | 0.005 64 | 0.018 | 1.240 | 3.138 | 0.716 |
| 12 | 880--1000 | 0.008 34 | 0.020 | 1.999 | 5.157 | 0.867 |

地表までの総光学的厚さは

$$
\tau_s=a\mu+\frac{bq_0}{4}\approx0.16+5.0=5.16
$$

である。

### パラメータの根拠

以前は水蒸気の項を持たず、オゾン以外の総光学的厚さを $\tau_s=1$ として気圧に比例して配分していた。この設定では上端から測った光学的厚さが $1$ に達するのは地表だけであり、大気上端から出る長波の大部分を地表と下層がそのまま放つので、温室効果が地球より大幅に弱かった。実際、slab ocean ケースの 5 年目の平衡状態は全球平均で $\braket{T_o}\approx275\,\mathrm{K}$、OLR $\approx239\,\mathrm{W\,m^{-2}}$（放射温度 $T_e\approx255\,\mathrm{K}$）で、温室効果は約 $20\,\mathrm{K}$ しかなかった。地球では $T_s\approx288\,\mathrm{K}$、温室効果は約 $33\,\mathrm{K}$ である。

新しい設定では、地球の全球平均に近い可降水量に対して $\tau_s\approx5.2$ となり、上端から測った光学的厚さが $1$ になる高度は $p/p_s\approx0.65$（約 $650\,\mathrm{hPa}$、高度約 $3.5\,\mathrm{km}$）で、雲を含む地球の実効放射高度（約 $5\,\mathrm{km}$）よりやや低いが同程度である。参考として、Frierson et al. (2006) の全球平均 $\tau_0$ は $4.5$、上の係数を地球の全球平均可降水量 $25\,\mathrm{kg\,m^{-2}}$ に当てると $\tau_s\approx5.1$ であり、地球相当の灰色大気の総光学的厚さは $4.5$--$6$ 程度と見込まれる。$a\mu$ が小さいので、水蒸気の少ない上部対流圏より上は長波に対してほぼ透明であり、放射冷却は水蒸気のある層に集中する。

放射対流平衡では、対流圏の気温減率が放射平衡より緩いので、地表温度は純放射平衡の $T_s^4=T_e^4(1+\tau_s/2)$（$\tau_s=5.16$ で $350\,\mathrm{K}$）よりはるかに低くなる。湿潤ケースで熱帯の最下層が $q_N\approx0.02$ に達すると最下層だけで $\Delta\tau_N\approx4.8$、気柱で $\tau_s\approx10$ 程度になりうる。灰色大気には大気の窓がないので、水蒸気による光学的厚さの増加は実際の地球より強く効く。したがって水蒸気フィードバックは地球より強めに出る傾向があり、暖かい気候では暴走しやすい。湿潤ケースでは日平均の $\braket{T_o}$ と OLR で平衡温度を監視し、5 年平均の $\braket{T_o}$ が地球の $288\,\mathrm{K}$ から大きく外れる場合は、$q_0$ ではなく $b$ または $\mu$ を見直す。$\alpha_s=0.3$ は雲を含む地球の惑星アルベドの代わりであり、雲の放射効果は長波にも入れない。

## 長波フラックス

上向き長波フラックスを $F^\uparrow_{k+1/2}$、下向き長波フラックスを $F^\downarrow_{k+1/2}$ とする。大気層の温度を $T_k$、Stefan–Boltzmann 定数を $\sigma=5.670374419\times10^{-8}\,\mathrm{W\,m^{-2}\,K^{-4}}$ とすれば、各層を通るフラックスは

$$
F^\uparrow_{k-1/2}
=t_kF^\uparrow_{k+1/2}+(1-t_k)\sigma T_k^4,
$$

$$
F^\downarrow_{k+1/2}
=t_kF^\downarrow_{k-1/2}+(1-t_k)\sigma T_k^4
$$

と計算する。上端から入射する長波はないものとし、上側の境界条件は

$$
F^\downarrow_{1/2}=0
$$

とする。地表は長波について黒体とし、浅い地面層の温度を $T_s$（slab ocean では $T_o$）として下側の境界条件を

$$
F^\uparrow_{N+1/2}=\sigma T_s^4
$$

とする。このため、下向き長波は地表ですべて吸収される。

正味長波フラックスを上向き正で

$$
F^{\mathrm{LW}}_{k+1/2}
=F^\uparrow_{k+1/2}-F^\downarrow_{k+1/2}
$$

と定義する。放射による第 $k$ 層の温度変化は

$$
\left(\frac{\partial T_k}{\partial t}\right)_{\mathrm{LW}}
=\frac{g}{c_p\Delta p_k}
\left(F^{\mathrm{LW}}_{k+1/2}-F^{\mathrm{LW}}_{k-1/2}\right)
$$

で与える。ただし $c_p=1004.0\,\mathrm{J\,kg^{-1}\,K^{-1}}$ である。

## 数値的な注意

- 光学的厚さが $1$ を超える層は、層全体が温度 $T_k$ の一様な放射体として扱われる。$N=12$ では最下層（$880$--$1000\,\mathrm{hPa}$）が基準分布で $\Delta\tau_N\approx2.1$、湿潤な熱帯では $5$ 程度になり、放射率がほぼ $1$ になる。このとき地表に届く下向き長波は $\sigma T_N^4$ にほぼ等しく、最下層内の温度構造は表現されない。層の放射と吸収は同じ $t_k$ で計算しているので、粗い鉛直解像度でもエネルギー収支は閉じたままである。
- 層の冷却率は $(1-t_k)/\Delta p_k\le(a\mu+bq_k)/p_0$ で抑えられるので、光学的厚さを増やしても $\Delta p_k$ の小さい上層で 1 ステップ当たりの冷却が発散することはない。乾燥した上層では $\Delta\tau_k\propto\Delta p_k$ となり、冷却率は $\Delta p_k$ によらず $(g/c_p)(a\mu/p_0)\sigma T_k^4$ 程度（$T_k=255\,\mathrm{K}$ で約 $0.3\,\mathrm{K\,day^{-1}}$）に収まる。
- $q_k$ に $\overline q^{+,n-1}_k$ を使うので、光学的厚さと放射フラックスはどちらも同じ前時刻の場から作られ、[物理過程を評価する時刻](./physics-time-level.md)の議論はそのまま成り立つ。
- 実装では、放射ルーチンは比湿の列を省略可能な引数として受け取る。湿潤ケースは前時刻の比湿の格子値を渡し、乾燥ケースは引数を省くことで基準水蒸気分布の層平均を内部で作る。基準分布の層平均を明示的に渡した場合の結果は、省いた場合とビット単位で一致する。

## 参考文献

- Frierson, D. M. W., I. M. Held, and P. Zurita-Gotor, 2006: A gray-radiation aquaplanet moist GCM. Part I: Static stability and eddy scale. *J. Atmos. Sci.*, **63**, 2548--2566.
- Byrne, M. P., and P. A. O'Gorman, 2013: Land--ocean warming contrast over a wide range of climates: Convective quasi-equilibrium theory and idealized simulations. *J. Climate*, **26**, 4000--4016.
- Vallis, G. K., et al., 2018: Isca, v1.0: a framework for the global modelling of the atmospheres of Earth and other planets at varying levels of complexity. *Geosci. Model Dev.*, **11**, 843--859.
