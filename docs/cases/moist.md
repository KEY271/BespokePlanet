# 湿潤ケース

[Slab ocean ケース](./slab-ocean.md)に水蒸気を加え、地軸の傾きを地球の値に戻した計算ケースの定義である。力学の傾向項は[湿潤大気](../dynamics/moist.md)のとおりに解き、全球を深さ 30 m の [slab ocean](../tendency/slab-ocean.md) で覆った aquaplanet とする。

## Slab ocean ケースとの差分

次の点を slab ocean ケースから変更する。

- 予報変数に比湿 $q$ を加え、[湿潤大気](../dynamics/moist.md)の方程式を解く。仮想温度、比湿の移流、比湿への超粘性と RAW フィルターを有効にする。
- 海面からの[蒸発](../tendency/evaporation.md)を有効にする。海面は常に飽和しているとみなし $\beta=1$ とする。潜熱フラックス $LE$ は海洋の熱収支から差し引く。
- [湿潤対流調節](../tendency/moist-convective-adjustment.md)と[大規模凝結](../tendency/large-scale-condensation.md)を有効にする。[乾燥対流調節](../tendency/dry-convective-adjustment.md)は[湿潤対流調節](../tendency/moist-convective-adjustment.md#乾燥対流調節の変更)の節の変更を加えた上で引き続き用いる。
- 地表摩擦は放射ケースから引き継ぐ [Held–Suarez 強制](../tendency/Held-Suarez.md#rayleigh-摩擦)の Rayleigh 摩擦（$\sigma_b=0.7$、時定数 $1\,\mathrm{day}$）をそのまま使い、拡散型の境界層スキームは持たない。
- [上層の Rayleigh 摩擦](../tendency/upper-rayleigh-friction.md)は slab ocean ケースから引き継ぐ。上端 2 層（$1$--$3\,\mathrm{hPa}$ と $3$--$10\,\mathrm{hPa}$）の風を時定数 $1\,\mathrm{day}$ で減衰させ、放射と湿潤過程による上層風の加速を抑える。全層一様の発散への超粘性（$\tau_\delta=1\,\mathrm{hour}$）も併用する。
- [長波放射](../tendency/longwave-radiation.md)の水蒸気による光学的厚さ $\Delta\tau^{\mathrm{H_2O}}_k=bq_k\Delta p_k/p_0$ に、slab ocean ケースの固定した[基準水蒸気分布](../tendency/longwave-radiation.md#基準水蒸気分布)ではなく、予報した比湿の正の部分 $\overline q_k^{+,n-1}$ を使う。これにより、比湿の変動が長波放射に反映される簡易的な水蒸気フィードバックを持つ。
- [雲](../tendency/cloud.md)の診断を有効にし、対流調節と大規模凝結の後の相対湿度と対流性降水から気柱の実効雲量を診断して、[短波放射](../tendency/shortwave-radiation.md)で下向き短波を一回だけ反射する（雲量 1 で $\alpha_c=0.43$）。雲による短波の吸収は扱わない。海面のアルベドは雲を繰り込んだ $0.3$ から開いた海面の $\alpha_o=0.06$ に変える。長波には雲を入れない。
- 地軸の傾きを slab ocean ケースの $0$ から[放射ケース](./radiation.md)と同じ $\varepsilon=23.4^\circ$（[暦と軌道](../calendar.md)）に戻す。したがって日変化と季節変化の両方を持ち、積分開始の 4 月 1 日は北半球の春分である。海洋の熱容量は 30 m の水柱のものなので、季節変化に対する海面温度の応答は地球の海洋より速い。

[オゾン](../tendency/ozone.md)、長波の係数 $a,b,\mu$、海面との顕熱交換、地表摩擦、[暦と軌道](../calendar.md)、自転角速度、計算期間、タイムステップおよび出力間隔は slab ocean ケースと同じとする。短波は水蒸気に依存させず、雲の放射効果は短波にだけ入れる。長波は灰色大気で大気の窓を持たないので、水蒸気フィードバックは地球より強めに出る傾向がある（[長波放射](../tendency/longwave-radiation.md#パラメータの根拠)）。Frierson 型の湿潤 GCM が持つ拡散型の境界層スキームは持たず、地表のバルクフラックス、対流調節、Rayleigh 摩擦で代替しているので、Frierson や Isca のモデルの再現ではない。すべての物理過程は[物理過程を評価する時刻](../tendency/physics-time-level.md)のとおり $\overline X^{n-1}$ の場で評価し、蒸発 → 対流調節・大規模凝結 → 雲量の診断 → 放射の順序（[湿潤大気](../dynamics/moist.md#物理過程)）で逐次に評価する。

## 解像度

他のケースは切断波数 $T=63$ で計算しているが、湿潤ケースは $T=31$ で計算する。比湿の追加で移流とスペクトル変換が増え、対流調節と大規模凝結の分だけ 1 ステップあたりの計算量も増えるためである。鉛直層数 $N=12$、タイムステップ $\Delta t=1200\,\mathrm{s}$、超粘性の時定数は $T=31$ でも変えない。超粘性の減衰率は切断波数 $n_\mathrm{max}$ で規格化してあるので、$T=31$ では最短波長の減衰時定数が同じまま、実長さスケールに対する減衰が強くなる。

## 初期値

$\zeta,\delta,T,\ln p_s$ と $\Phi_s=0$、海洋温度 $T_o$ の初期値は slab ocean ケースと同じとする。比湿の初期値は、相対湿度 $\mathrm{RH}_0=0.7$ の対流圏として

$$
q^0_k=
\begin{cases}
\mathrm{RH}_0\,q_s(T^0_k,p_k),&p_k\ge p_{\mathrm{top}},\\
0,&p_k<p_{\mathrm{top}},
\end{cases}
\qquad p_{\mathrm{top}}=200\,\mathrm{hPa}
$$

とする。$T^0_k$ は Jablonowski–Williamson（JW）の温度、$p_k$ は full level の気圧、$q_s$ は[飽和比湿](../tendency/saturation-specific-humidity.md)である。$p_{\mathrm{top}}$ は JW の対流圏界面 $\eta_t=0.2$ に対応する。JW の成層圏は上端に向かって暖かく（$2\,\mathrm{hPa}$ で約 $262\,\mathrm{K}$）、そこでは $p\le e_s$ となって飽和比湿が $1$ に張り付き意味を持たないので、成層圏は乾燥した状態から始める。$q^0$ は格子で計算してからスペクトルに変換するので、切断誤差により格子で見た初期比湿は $\mathrm{RH}_0q_s$ からわずかにずれ、局所的に負になりうる。

この初期値は JW の温度場に対して定めるので、最下層の温度が約 $307\,\mathrm{K}$ に達する JW の赤道では $q^0_N\approx0.026$、可降水量約 $81\,\mathrm{kg\,m^{-2}}$、水蒸気による長波光学的厚さ約 $16$ となり、極（最下層約 $224\,\mathrm{K}$）では可降水量約 $2\,\mathrm{kg\,m^{-2}}$ になる。全球平均では可降水量約 $46\,\mathrm{kg\,m^{-2}}$、光学的厚さ約 $9$ で、slab ocean ケースが使う[基準水蒸気分布](../tendency/longwave-radiation.md#基準水蒸気分布)（$25.5\,\mathrm{kg\,m^{-2}}$、$5.0$）より湿っている。JW は力学のテスト用の基本場であって放射対流平衡ではないので、この差は積分の初期に降水と放射によって調整される。$q^0=0$ から始める場合のように大気が一時的に長波に透明になることはないが、初期の数十日は平衡状態の統計に含めない。

$q^0\neq0$ なので仮想温度は温度より $\delta_vq^0T$ だけ高くなる。その大きさは赤道の最下層で約 $5\,\mathrm{K}$、緯度 $45^\circ$ で約 $0.6\,\mathrm{K}$、極ではほぼ $0$ である。[放射ケース](./radiation.md)で $\Phi_s=0$ に釣り合わせた東西風はこの分だけ釣り合いからずれるが、ずれは JW の温度場の南北コントラスト（約 $80\,\mathrm{K}$）に比べて小さく、重力波として数日で調整されるので、5 年間の積分では無視できる。

## 時間

slab ocean ケースと同じく $5$ 年間、すなわち $5\times360$ 太陽日シミュレーションする。タイムステップは $\Delta t=1200\,\mathrm{s}$ とする。

## 出力

slab ocean ケースの出力に、水蒸気に関する量を追加する。降水と蒸発は $\mathrm{kg\,m^{-2}\,s^{-1}}$ で集計し、ファイルには $\mathrm{mm\,day^{-1}}$（$86400$ 倍）で書く。可降水量は

$$
W=\sum_{k=1}^N\frac{q_k^+\Delta p_k}{g}
$$

で、単位は $\mathrm{kg\,m^{-2}}$ である。[負の比湿](../dynamics/moist.md#負の比湿)を許しているので、水収支の確認には符号付きの気柱水蒸気量と負の部分の量

$$
W_\pm=\sum_{k=1}^N\frac{q_k\Delta p_k}{g},\qquad
W_-=-\sum_{k=1}^N\frac{\min(q_k,0)\Delta p_k}{g},\qquad
W=W_\pm+W_-
$$

も別に集計する。物理過程の傾向は符号付きの $q$ に加わるので、$E-P$ と対応するのは $W$ ではなく $W_\pm$ である。

1 日ごとの日平均量には、次を追加する。

- 全球平均の降水 $\braket{P}=\braket{P_{\mathrm{conv}}+P_{\mathrm{ls}}}$、およびその内訳 $\braket{P_{\mathrm{conv}}}$、$\braket{P_{\mathrm{ls}}}$
- 全球平均の蒸発 $\braket{E}$
- 全球平均の潜熱フラックス $\braket{LE}$
- 全球平均の可降水量 $\braket{W}$、符号付きの気柱水蒸気量 $\braket{W_\pm}$、負の部分 $\braket{W_-}$
- 全球平均の全雲量 $\braket{C}$

これらは放射の全球平均と同じくオンラインで集計する。定常状態では $\braket{P}$ と $\braket{E}$ が一致するので、その差が水収支の指標になる。また日平均の $\braket{W_\pm}$ の日々の差分と $\braket{E}-\braket{P}$ の差は、移流の離散化と RAW フィルターによる水蒸気の非保存の大きさを表す。

1 日ごとの統計量には、日平均のほかに、その日の最大風速とその場所を追加する。

- その日の最大風速 $|\bm{u}|_{\mathrm{max}}=\max\sqrt{u^2+v^2}$
- 最大をとった格子点の経度 $\lambda$ と緯度 $\varphi$（度）
- 最大をとった鉛直層の番号 $k$ と、その full level の $\eta_k$

最大は、その日に含まれる 72 ステップすべてについて、全格子点・全鉛直層にわたって取る。$u,v$ は非線形項の評価のためどのステップでも格子上にあるので、追加のスペクトル変換はいらない。同じ最大値が複数の格子点で出た場合は、走査順で最初に見つかった格子点を記録する。格子は[八面体ガウス格子](../dynamics/octahedral-gaussian-grid.md)で緯度ごとに経度の数が違うので、添字ではなく経度・緯度の値そのものを書く。この最大風速は、上層の Rayleigh 摩擦が有効な湿潤ケースでも上層風と移流 CFL の安定性を監視するために使う。

30 日ごとの月平均量には、次を追加する。

- 各格子点の降水の月平均 $P(\lambda,\varphi)$
- 各格子点の蒸発の月平均 $E(\lambda,\varphi)$
- 各格子点の可降水量の月平均 $W(\lambda,\varphi)$
- 各格子点の全雲量の月平均 $C(\lambda,\varphi)$
- 経度方向に平均した比湿の月平均 $[q](\varphi,\eta)$
- 偏差の積の月平均 $[v'q'](\varphi,\eta)$

毎年の 4/1 の瞬時値に以下を追加する。

- 現時刻 $X^n$ のスペクトル比湿 $q$（格子の $q$ も合わせて出力する）
- 現時刻 $X^n$ の海洋温度 $T_o$。$T_o$ は格子量なので格子の値だけを書く（放射ケースの $T_s$ と同じファイル名）。
- 格子の雲量 $C(\lambda,\varphi)$。予報変数ではないので格子の値だけを書き、その時刻で終わるステップの傾向評価で診断した直近の値を使う。初期の瞬時値では 0 である。
- 時刻、すなわちモデル時刻（積分開始からの経過秒）とステップ数。日変化と軌道上の位置、月・年の境界はここから定まる。

ファイル名の規則は slab ocean ケースと同じとする。出力先は `moist_slab_ocean` である。
