# 湿潤ケース

[Slab ocean ケース](./slab-ocean.md)に水蒸気を加えた計算ケースの定義である。力学の傾向項は[湿潤大気](../dynamics/moist.md)のとおりに解き、全球を深さ 30 m の [slab ocean](../tendency/slab-ocean.md) で覆った aquaplanet とする。

## Slab ocean ケースとの差分

次の点を slab ocean ケースから変更する。

- 予報変数に比湿 $q$ を加え、[湿潤大気](../dynamics/moist.md)の方程式を解く。仮想温度、比湿の移流、比湿への超粘性と RAW フィルターを有効にする。
- 海面からの[蒸発](../tendency/evaporation.md)を有効にする。海面は常に飽和しているとみなし $\beta=1$ とする。潜熱フラックス $LE$ は海洋の熱収支から差し引く。
- [湿潤対流調節](../tendency/moist-convective-adjustment.md)と[大規模凝結](../tendency/large-scale-condensation.md)を有効にする。[乾燥対流調節](../tendency/dry-convective-adjustment.md)は[湿潤対流調節](../tendency/moist-convective-adjustment.md#乾燥対流調節の変更)の節の変更を加えた上で引き続き用いる。
- 地表摩擦は放射ケースから引き継ぐ [Held–Suarez 強制](../tendency/Held-Suarez.md#rayleigh-摩擦)の Rayleigh 摩擦（$\sigma_b=0.7$、時定数 $1\,\mathrm{day}$）をそのまま使い、拡散型の境界層スキームは持たない。
- [上層の Rayleigh 摩擦](../tendency/upper-rayleigh-friction.md)は使わない。上端 2 層の風を陽に減衰させることはせず、上層の発散の抑制は全層一様の超粘性（$\tau_\delta=1\,\mathrm{hour}$）だけに任せる。

[長波放射](../tendency/longwave-radiation.md)、[短波放射](../tendency/shortwave-radiation.md)、[オゾン](../tendency/ozone.md)、海面との顕熱交換、地表摩擦、[暦と軌道](../calendar.md)、自転角速度、地軸の傾き $0$、計算期間、タイムステップおよび出力間隔は slab ocean ケースと同じとする。放射は水蒸気の量に依存させない。すなわち長波の光学的厚さ $\tau_s=1$ は水蒸気の平均的な効果を含んだ固定値とみなし、水蒸気の変動による放射のフィードバックと雲の放射効果は扱わない。したがってこのケースは現実的な気候感度を持つモデルではなく、潜熱と大規模力学の相互作用を放射フィードバックから切り離して調べる、湿潤力学の理想化モデルである。Frierson 型の湿潤 GCM が持つ拡散型の境界層スキームも持たず、地表のバルクフラックス、対流調節、Rayleigh 摩擦で代替しているので、Frierson や Isca のモデルの再現ではない。すべての物理過程は[物理過程を評価する時刻](../tendency/physics-time-level.md)のとおり $\overline X^{n-1}$ の場で評価し、対流調節と大規模凝結は[湿潤大気](../dynamics/moist.md#物理過程)の順序で逐次に評価する。

## 解像度

他のケースは切断波数 $T=63$ で計算しているが、湿潤ケースはまず $T=31$ で計算する。比湿の追加で移流とスペクトル変換が増え、対流調節と大規模凝結の分だけ 1 ステップあたりの計算量も増えるので、まずは $T=31$ の 5 年積分で水収支と気候が破綻していないことを確かめ、その後で $T=63$ に上げる。鉛直層数 $N=12$、タイムステップ $\Delta t=1200\,\mathrm{s}$、超粘性の時定数は $T=31$ でも変えない。超粘性の減衰率は切断波数 $n_\mathrm{max}$ で規格化してあるので、$T=31$ では最短波長の減衰時定数が同じまま、実長さスケールに対する減衰が強くなる。

## 初期値

$\zeta,\delta,T,\ln p_s$ と $\Phi_s=0$、海洋温度 $T_o$ の初期値は slab ocean ケースと同じとする。比湿の初期値は全層で

$$
q^0=0
$$

とする。$q^0=0$ では仮想温度が温度に一致するので、[放射ケース](./radiation.md)で $\Phi_s=0$ に釣り合わせた東西風は初期に釣り合ったままである。水蒸気は海面からの蒸発によって供給され、数十日程度で大気に行き渡る。5 年間の積分ではこの spin-up は無視できる。

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

これらは放射の全球平均と同じくオンラインで集計する。定常状態では $\braket{P}$ と $\braket{E}$ が一致するので、その差を水収支の確認に使う。また日平均の $\braket{W_\pm}$ の日々の差分と $\braket{E}-\braket{P}$ の差は、移流の離散化と RAW フィルターによる水蒸気の非保存の大きさを表すので、これが $\braket{P}$ に比べて十分小さいことを確認する。$\braket{W_-}$ が $\braket{W}$ に比べて無視できないほど大きくなれば、正値性を保つ移流や水分補正の導入を検討する。

1 日ごとの統計量には、日平均のほかに、その日の最大風速とその場所を追加する。

- その日の最大風速 $|\bm{u}|_{\mathrm{max}}=\max\sqrt{u^2+v^2}$
- 最大をとった格子点の経度 $\lambda$ と緯度 $\varphi$（度）
- 最大をとった鉛直層の番号 $k$ と、その full level の $\eta_k$

最大は、その日に含まれる 72 ステップすべてについて、全格子点・全鉛直層にわたって取る。$u,v$ は非線形項の評価のためどのステップでも格子上に用意されているので、その値を走査して最大値とその添字を更新すればよく、追加のスペクトル変換はいらない。同じ最大値が複数の格子点で出た場合は、走査順で最初に見つかった格子点を記録する。格子は[八面体ガウス格子](../dynamics/octahedral-gaussian-grid.md)で緯度ごとに経度の数が違うので、添字ではなく経度・緯度の値そのものを書く。この最大風速は、湿潤ケースで上層の Rayleigh 摩擦を外したことによる上層風の暴走を見張るためにも使う。

30 日ごとの月平均量には、次を追加する。

- 各格子点の降水の月平均 $P(\lambda,\varphi)$
- 各格子点の蒸発の月平均 $E(\lambda,\varphi)$
- 各格子点の可降水量の月平均 $W(\lambda,\varphi)$
- 経度方向に平均した比湿の月平均 $[q](\varphi,\eta)$
- 偏差の積の月平均 $[v'q'](\varphi,\eta)$

毎年の 4/1 の瞬時値に以下を追加する。

- 現時刻 $X^n$ のスペクトル比湿 $q$（格子の $q$ も合わせて出力する）
- 現時刻 $X^n$ のスペクトル海洋温度 $T_o$。
- 時刻、すなわちモデル時刻（積分開始からの経過秒）とステップ数。日変化と軌道上の位置、月・年の境界はここから定まる。

ファイル名の規則は slab ocean ケースと同じとする。

## 実装後の確認項目

湿潤ケースを走らせる前に、次を単体試験で確認する。

- [大規模凝結](../tendency/large-scale-condensation.md)の後で各層の $|q_k-q_s(T_k,p_k)|\le10^{-4}q_s$ となり、層ごとの湿潤エンタルピー $c_pT_k+Lq_k$ が変わらないこと。
- 地面と最下層の温位が等しい乾燥断熱中立な気柱で、顕熱フラックス $H_{sa}$（$H_{oa}$）がちょうど 0 になること。
- 空気塊に正の浮力を持つ層がない気柱、および $P_T\le0$ の気柱で、[湿潤対流調節](../tendency/moist-convective-adjustment.md)の傾向がすべて 0 になること。
- 深い対流で $\sum_k(c_pC^{\mathrm{conv}}_{T,k}+LC^{\mathrm{conv}}_{q,k})\Delta p_k=0$ と $P_{\mathrm{conv}}=-\sum_kC^{\mathrm{conv}}_{q,k}\Delta p_k/g>0$、浅い対流で $\sum_kC^{\mathrm{conv}}_{q,k}\Delta p_k=0$ と $\sum_kC^{\mathrm{conv}}_{T,k}\Delta p_k=0$ が丸め誤差の範囲で成り立つこと。
- 数ステップの全球積分で $\braket{W_\pm}$ の変化が $\braket{E}-\braket{P}$ に、移流の離散化誤差の範囲で一致すること。
