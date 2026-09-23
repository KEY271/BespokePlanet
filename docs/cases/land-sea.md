# 陸海ケース

[湿潤ケース](./moist.md)の aquaplanet に、[地形](../dynamics/topography.md)で生成した大陸と山脈を加え、地表を[陸と海の混合](../tendency/land-sea-surface.md)で扱う計算ケースの定義である。力学の傾向項は[湿潤大気](../dynamics/moist.md)のとおりに解く。$T=31$ と $T=63$ の両方で同じ惑星を計算する。

## 湿潤ケースとの差分

次の点を湿潤ケースから変更する。

- 地表ジオポテンシャル $\Phi_s$ に、[地形](../dynamics/topography.md#既定の惑星)の既定の惑星から生成した値を使う。$\Phi_s$ は積分中ずっと固定する。
- 同じ生成から得た陸面率 $f_L$ を[陸海タイル](../tendency/land-sea-surface.md)に渡し、各温度で計算したフラックスを面積平均する。陸は[二層の地面](../tendency/ground.md)と[バケツ](../tendency/bucket.md)、海は深さ 30 m の [slab ocean](../tendency/slab-ocean.md) と[海氷](../tendency/sea-ice.md)である。湿潤ケースと同じ [Q flux](../tendency/q-flux.md) を海だけに与え、海面積平均を引いて全球の海への積分を 0 にする。$T_L,T_d,T_o,A,V,W,S$ を予報する。
- [降雪と積雪](../tendency/snow.md)を有効にする。$263.15\,\mathrm K$ より低い層の大規模凝結を雪とし、$278.15\,\mathrm K$ より暖かい層で融かす。地面に届いた雪は陸の積雪 $S$ にため、積雪率 $f=S/(S+S_0)$ で陸のアルベドを $0.4f$ 上げ、陸の顕熱・蒸発を $1-f$ 倍にし、陸温度が $275\,\mathrm K$ を超えた分の熱で融かす。海に落ちた雪は混合層の熱で融かし、海氷上には積もらせない。
- 初期値は[地形上の初期状態](../dynamics/topography.md#地形上の初期状態)のとおり、放射ケースの基本場を気圧の関数として読み、地表面気圧を $\Phi^\ast(\varphi,p_s)=\Phi_s$ から決める。
- 出力に陸面率・地表高度の静的な場と、陸と海で分けた全球平均を加える。

[雲](../tendency/cloud.md)、[短波放射](../tendency/shortwave-radiation.md)、[長波放射](../tendency/longwave-radiation.md)、[オゾン](../tendency/ozone.md)、[乾燥対流調節](../tendency/dry-convective-adjustment.md)、[湿潤対流調節](../tendency/moist-convective-adjustment.md)、[大規模凝結](../tendency/large-scale-condensation.md)、地表摩擦、[上層の Rayleigh 摩擦](../tendency/upper-rayleigh-friction.md)、[暦と軌道](../calendar.md)（地軸の傾き $23.4^\circ$）、自転角速度、計算期間、タイムステップ、超粘性の時定数は湿潤ケースと同じとする。すべての物理過程は[物理過程を評価する時刻](../tendency/physics-time-level.md)のとおり $\overline X^{n-1}$ の場で評価する。

## 地表のパラメータ

| 量 | 陸 | 海 |
| --- | --- | --- |
| 熱容量 | $C_s=2\times10^6\,\mathrm{J\,m^{-2}\,K^{-1}}$（浅い層）、$C_d=2\times10^7$（深い層） | $C_o=1.2558\times10^8\,\mathrm{J\,m^{-2}\,K^{-1}}$ |
| 地中との熱交換 | $K_{sd}=2\,\mathrm{W\,m^{-2}\,K^{-1}}$ | なし |
| 短波アルベド（雲を含まない） | $\alpha_L=0.2$ | $\alpha_o=0.06$ |
| 水 | $W_{\max}=150\,\mathrm{kg\,m^{-2}}$、$W_0=75\,\mathrm{kg\,m^{-2}}$、$\beta_L=W/W_{\max}$ | 無限、$\beta_o=1$ |
| 雪 | 積雪 $S$（水当量、初期値 0）、$S_0=50\,\mathrm{kg\,m^{-2}}$、アルベド $\alpha_L+0.4f$、顕熱・蒸発 $\times(1-f)$、融雪 $T_L>275\,\mathrm K$ | ためない。混合層の熱で融かす |

熱容量は陸・海それぞれの面積当たりとし、海岸でも混合しない。蒸発・顕熱・上向き長波は各温度で計算して面積平均する。氷のアルベドは 0.60、昇華・凝華は 0。雲の反射は別に計算するので、アルベドは地表そのものの値である。$W_{\max}$ が陸の乾燥の時間スケールを決める主要なパラメータである。

## 解像度

湿潤ケースと同じ $T=31$ に加えて、$T=63$ でも計算する。両方の解像度で同じ地形の定義を使うので、$z_s$ と $f_L$ の切断の影響（負の行き過ぎ、RMS 差）は $T=63$ で小さい。鉛直層数 $N=12$、タイムステップ $\Delta t=1200\,\mathrm{s}$ はどちらの解像度でも変えない。$T=63$ の放射ケースが同じ $\Delta t$ で走っているので、地形による重力波は陰解法で扱える範囲にある。移流の CFL は最大風速の診断で監視する。

$T=31$ と $T=63$ はそれぞれコマンドライン引数 `land` と `land-t63` で走らせ、出力先は `moist_land_sea_t31` と `moist_land_sea_t63` である。

## 初期値

[地形上の初期状態](../dynamics/topography.md#地形上の初期状態)のとおり、格子の $\Phi_s$ から各格子点の $p_s$ を決め、full level の気圧で JW の温度と釣り合った東西風を評価する。比湿の初期値は湿潤ケースと同じ相対湿度 $\mathrm{RH}_0=0.7$、$p_{\mathrm{top}}=200\,\mathrm{hPa}$ とし、各格子点の $p_s$ から作った full level の気圧を使う。陸の $T_L,T_d$ は最下層の温度 $T_N$ に等しくし、海洋混合層の温度 $T_o$ は[湿潤ケース](./moist.md#初期値)と同じ解析的な分布 $T_p+(T_e-T_p)\cos^2\varphi$ とする。海氷は $A=V=0$ から始める。陸温度には凍結制限を加えない。陸を含む格子のバケツは $W=W_0=75\,\mathrm{kg\,m^{-2}}$、純海洋格子では $W=0$ とする。積雪は全格子で $S=0$ から始める。海上では $\Phi_s$ の切断のリップルを除いて湿潤ケースの初期値と同じである。

## 時間

湿潤ケースと同じく $5$ 年間、すなわち $5\times360$ 太陽日シミュレーションする。タイムステップは $\Delta t=1200\,\mathrm{s}$ とする。

## 出力

湿潤ケースの出力に次を加える。ファイル名の規則、日平均・月平均・毎年の瞬時値の定義は湿潤ケースと同じとする。

積分の開始時に一度だけ、時間に依らない格子の場を書く。

- 陸面率 $f_L(\lambda,\varphi)$（`land_fraction.bin`）
- 海面積当たりの [Q flux](../tendency/q-flux.md#出力) $Q(\lambda,\varphi)$（`ocean_q_flux.bin`）
- 格子の地表高度 $z_s(\lambda,\varphi)$（切断後、m、`surface_height.bin`）

1 日ごとの日平均量には、[陸と海の混合](../tendency/land-sea-surface.md#診断)の定義で次を加える。全球平均の表面温度は陸・開水面・氷表面の面積平均である。

- 陸温度と海水温の平均 $\braket{T_L}_L,\braket{T_o}_O$
- 陸上・海上の降水 $\braket{P}_L,\braket{P}_O$
- 陸上・海上の蒸発 $\braket{E}_L,\braket{E}_O$
- 陸面貯水量 $\braket{W}_L$、蒸発効率 $\braket{\beta_L}_L$、流出 $\braket{R}_L$、乾燥域率
- バケツの水収支残差の陸面平均と最大絶対値
- 地中温度の陸面積平均 $\braket{T_d}_L$
- 降雪・大気中の融雪の全球平均、陸の降雪 $\braket{S_{\rm ls}}_L$、積雪 $\braket{S}_L$、積雪率 $\braket{f}_L$、融雪 $\braket{M_s}_L$、積雪の収支残差の最大絶対値（[降雪と積雪の出力](../tendency/snow.md#10-出力と診断)）

30 日ごとの月平均量には次を加える。

- 各格子点の地中温度の月平均 $T_d(\lambda,\varphi)$
- 各格子点の地上気温 $T_a=T_N(p_s/p_N)^\kappa$ の月平均（`monthly_surface_air_temperature_m{month:04d}.bin`、K、[地面](../tendency/ground.md)の顕熱が使う外挿温度）
- 各格子点の陸面貯水量 $W(\lambda,\varphi)$、蒸発効率 $\beta_L(\lambda,\varphi)$、流出 $R(\lambda,\varphi)$ の月平均
- 各格子点の積雪 $S$、積雪率 $f$、降雪 $S_{\rm ls}$、融雪 $M_s$ の月平均

毎年の 4/1 の瞬時値には、放射ケースと同じく格子の $T_d$ と、再開用の格子の $W$、$S$ と積雪率 $f$ を加える。$T_s,T_d,W,S$ はいずれも格子量なのでスペクトルの出力はない。

メタデータには、地形の形状一覧とパラメータ、全球の陸面率 $\braket{f_L}$、格子の $z_s$ の最大・最小と解析形との RMS 差、陸と海のパラメータの表、バケツの $W_{\max}$ と $W_0$、雪のパラメータ、切断波数を書く。
