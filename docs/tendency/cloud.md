# 雲

雲を予報変数としては持たず、各タイムステップで格子点ごとに気柱一つの実効的な雲量（気柱の面積のうち雲が占める割合）$C\in[0,1]$ を、相対湿度と降水から診断する。層ごとの雲量や雲の高さは区別しない。凝結した水は[大規模凝結](./large-scale-condensation.md)と[湿潤対流調節](./moist-convective-adjustment.md)のとおり直ちに降水として地表へ落ちるので、雲水・雲氷の量は追跡せず、雲は放射のためだけの診断量である。灰色スキームでは、診断した雲量を[短波放射](./shortwave-radiation.md)で下向き短波を反射するためにだけ使い、[長波放射](./longwave-radiation.md)には影響させない。[帯域別放射](./band-radiation.md)では、雲量を大規模雲と対流雲に分けてそれぞれ高さを持たせ、長波と短波の両方に効かせる（[帯域別放射の雲](#帯域別放射の雲)）。鉛直層の添字は[長波放射](./longwave-radiation.md)、full level の気圧 $p_k$ と飽和比湿 $q_s$ は[飽和比湿](./saturation-specific-humidity.md)と同じとする。

## 診断に使う場

雲量は、[湿潤大気](../dynamics/moist.md#物理過程)の順序のとおり乾燥対流調節・湿潤対流調節・大規模凝結の傾向を $2\Delta t$ かけて加えた暫定場

$$
T^{(3)}_k=T^{(2)}_k+2\Delta tC^{\mathrm{ls}}_{T,k},\qquad
q^{(3)}_k=\max\left(q^{(2)}_k+2\Delta tC^{\mathrm{ls}}_{q,k},\,0\right)
$$

から診断する。$T^{(2)}_k,q^{(2)}_k$ は[大規模凝結](./large-scale-condensation.md#評価する場)が評価に使った暫定場で、$p_k$ は $(\ln p_s)^{\overline{n-1}}$ から作る。大規模凝結の後なので、各層は許容誤差の範囲で飽和を超えていない（$q^{(3)}_k\le(1+10^{-4})q_s(T^{(3)}_k,p_k)$）。降水には同じステップの湿潤対流調節の対流性降水 $P_{\mathrm{conv}}$ を使う。

前時刻の場 $\overline X^{n-1}$ ではなく暫定場を使うのは、対流調節と大規模凝結が同じステップで除いた不安定と過飽和を反映した相対湿度で雲を診断し、大規模凝結で飽和に達した気柱をそのまま曇天として扱うためである。

## 相対湿度による雲量

各層の相対湿度を

$$
\mathrm{RH}_k=\min\left(1,\;\frac{q^{(3)}_k}{q_s\!\left(T^{(3)}_k,p_k\right)}\right)
$$

とし、気柱の最大値

$$
\mathrm{RH}_{\max}=\max_{1\le k\le N}\mathrm{RH}_k
$$

から、Slingo (1987) の層状雲の式の形で

$$
C^{\mathrm{RH}}=\left[\max\left(0,\;\frac{\mathrm{RH}_{\max}-\mathrm{RH}_c}{1-\mathrm{RH}_c}\right)\right]^2,\qquad
\mathrm{RH}_c=0.8
$$

とする。$\mathrm{RH}_{\max}\le\mathrm{RH}_c$ では雲はなく、どこかの層が飽和していれば $C^{\mathrm{RH}}=1$ になる。二乗にするのは $\mathrm{RH}_c$ の近くで雲量が滑らかに立ち上がるようにするためである。気柱の最大値を使うのは、雲がもっとも湿った層に生じ、その層が気柱を覆うと考えるためで、層ごとに雲量を作って重ねる代わりに最大重なりを仮定したことに相当する。

$\mathrm{RH}_c$ は湿潤対流調節の参照相対湿度 $0.7$ より大きくとる。こうすると、対流調節で参照プロファイルに緩和された気柱がそれだけで雲に覆われることはなく、対流域の雲は後述の対流性の雲量が担う。大規模凝結が起きた気柱は飽和した層を持つので $C^{\mathrm{RH}}=1$ であり、大規模な降水域は曇天になる。

$p_k\le e_s(T^{(3)}_k)$ となって $q_s=1$ に固定される暖かい成層圏では $\mathrm{RH}_k\approx0$ で、最大値に寄与しない。負の比湿の層は $q^{(3)}_k=0$ とみなす。高度による $\mathrm{RH}_c$ の変化や、下層の逆転層による層積雲の補正（Isca の SimCloud が持つもの）は入れない。

## 対流性降水による雲量

対流性の雲量は Slingo (1987) に従い、対流性降水 $P_{\mathrm{conv}}$（$\mathrm{mm\,day^{-1}}$ に換算した値）の対数から

$$
C^{\mathrm{conv}}=
\begin{cases}
\min\left\{C^{\mathrm{conv}}_{\max},\;\max\left[0,\;c_0+c_1\ln\dfrac{P_{\mathrm{conv}}}{P_0}\right]\right\},&P_{\mathrm{conv}}>0,\\[8pt]
0,&P_{\mathrm{conv}}=0
\end{cases}
$$

とする。係数は

$$
c_0=0.245,\qquad c_1=0.125,\qquad P_0=1\,\mathrm{mm\,day^{-1}},\qquad C^{\mathrm{conv}}_{\max}=0.8
$$

である。$P_{\mathrm{conv}}<0.14\,\mathrm{mm\,day^{-1}}$ では 0、$1\,\mathrm{mm\,day^{-1}}$ で $0.245$、$10\,\mathrm{mm\,day^{-1}}$ で $0.53$、$85\,\mathrm{mm\,day^{-1}}$ 以上で上限の $0.8$ になる。[浅い対流](./moist-convective-adjustment.md#浅い対流)は $P_{\mathrm{conv}}=0$ なので雲を作らない。

## 実効雲量

気柱の実効雲量は両者の大きい方

$$
C=\max\left(C^{\mathrm{RH}},\;C^{\mathrm{conv}}\right)
$$

とする。$0\le C\le1$ である。[短波放射](./shortwave-radiation.md)はこの $C$ を使い、地表に達する前の下向き短波のうち $C\alpha_c$（$\alpha_c=0.43$）を一回だけ反射する。雲の高さは持たないので、反射する位置による加熱の違いは表れない。全層で $\mathrm{RH}_k\le\mathrm{RH}_c$ かつ $P_{\mathrm{conv}}=0$ の気柱では $C=0$ となり、短波のフラックスは雲を持たない式とビット単位で一致する。反射短波は雲の反射と地表の反射の和 $[C\alpha_c+(1-C\alpha_c)\alpha_s]F^{\downarrow\mathrm{SW}}$ で、$C$ について単調増加なので短波入射の $\alpha_c+(1-\alpha_c)\alpha_s$ 倍を超えない。

雲は比湿と対流性降水を必要とするので、比湿を予報するケースでだけ有効にできる。比湿を予報しない[放射ケース](../cases/radiation.md)と [slab ocean ケース](../cases/slab-ocean.md)では雲を診断せず、$C=0$ として短波を計算する。

## 帯域別放射の雲

[帯域別放射](./band-radiation.md#6-雲の種類と副気柱)では、上の $C^{\mathrm{RH}}$ と $C^{\mathrm{conv}}$ から、気柱を晴天・大規模雲・対流雲の 3 つの副気柱に分ける。

$$
a_{\mathrm{cv}}=C^{\mathrm{conv}},\qquad
a_{\mathrm{ls}}=\max\left(C^{\mathrm{RH}}-C^{\mathrm{conv}},\,0\right),\qquad
a_{\mathrm{clr}}=1-a_{\mathrm{ls}}-a_{\mathrm{cv}}
$$

面積率の和 $a_{\mathrm{ls}}+a_{\mathrm{cv}}$ は実効雲量 $C$ に等しい。二つの雲が水平に重なる部分は、雲頂の高い対流雲の副気柱に含める。

- 対流雲は、同じステップの[湿潤対流調節](./moist-convective-adjustment.md)の対流頂 $k_t$ の層に置く。$C^{\mathrm{conv}}>0$ は降水のある深い対流のときだけなので、そのとき $k_t$ は定まっている。
- 大規模雲は、$\mathrm{RH}_k$ が気柱の最大値 $\mathrm{RH}_{\max}$ をとる層のうち最も上の層に置く。

放射の扱いは[帯域別放射](./band-radiation.md#63-雲の放射特性)のとおりで、どちらの雲も長波では黒体、短波では雲頂にある吸収のない δ-Eddington の層（光学的厚さ $\tau=3.9$、非対称因子 $g=0.85$）とする。$\tau$ は ISCCP の全球平均の雲の光学的厚さ、$g$ は水滴の散乱の値で、大気上端の収支に合わせて調整した値ではない。雲量の出力 $C$ は変わらない。帯域別放射を使うケースでは、面積率 $a_{\mathrm{ls}}$、$a_{\mathrm{cv}}$ と雲の層の気圧も出力する（同ページの 10.5 節）。

## 評価する時刻と順序

雲は[物理過程を評価する時刻](./physics-time-level.md)の例外として、$\overline X^{n-1}$ そのものではなく、そこから同じステップの対流調節と大規模凝結で進めた暫定場と、同じステップの対流性降水から診断する。したがって物理過程の評価順は

1. [蒸発](./evaporation.md)
2. [乾燥対流調節](./dry-convective-adjustment.md) → [湿潤対流調節](./moist-convective-adjustment.md) → [大規模凝結](./large-scale-condensation.md)
3. 雲量の診断
4. [陸面のバケツ](./bucket.md)の更新（陸面を持つケースのみ）
5. [長波放射](./longwave-radiation.md)・[短波放射](./shortwave-radiation.md)と地表の熱収支

とする。放射は雲量以外の入力（温度、比湿、地表温度、風）を従来どおり $\overline X^{n-1}$ からとり、雲量だけを 3 の結果から受け取る。雲量は暫定場から作るが放射の傾向は $\overline X^{n-1}$ に加わるので、LeapFrog の計算モードを増幅させる減衰項にはならない。暫定場に使う時間幅 $2\Delta t$ は、最初の 2 ステップでは[大規模凝結](./large-scale-condensation.md)と同じく $F$ の第 1 引数を $\Delta t$ とみなした値、すなわち各段階が実際に進める幅に読み替える。

雲量は時間方向に平滑化せず、毎ステップ診断し直す。$\mathrm{RH}_c$ 付近で雲量が滑らかに立ち上がるようにしてあるので、雲量がステップごとに 0 と 1 の間で振動することはない。

## 診断

雲量 $C(\lambda,\varphi)$ を格子で持ち、[湿潤ケース](../cases/moist.md#出力)の出力には日平均の全球平均 $\braket{C}$、月平均の $C(\lambda,\varphi)$、毎年の瞬時値の $C(\lambda,\varphi)$ が含まれる。瞬時値の $C$ は予報変数ではないので、瞬時値の時刻で終わるステップの傾向評価（その 1 ステップ前の場から作った暫定場）で診断した直近の値を書き、初期の瞬時値では未診断なので 0 を書く。惑星アルベド $\braket{F^{\mathrm{SW}}_{\mathrm{refl}}}/\braket{F^{\downarrow\mathrm{SW}}_{1/2}}$ は短波入射と反射短波の出力から求める。

## パラメータの根拠

- $\mathrm{RH}_c=0.8$ は Slingo (1987) の低層・高層の層状雲の値で、湿潤対流調節の参照相対湿度 $0.7$ より大きくとる。
- $c_0,c_1,C^{\mathrm{conv}}_{\max}$ は Slingo (1987) の対流性の雲の値をそのまま使う。
- 灰色スキームの雲の短波の反射率 $\alpha_c=0.43$ は[短波放射](./shortwave-radiation.md#雲による反射)で与える。地球の全球平均の雲量 $0.6$–$0.7$ に対して、雲の反射だけで惑星アルベドの $0.26$–$0.30$ を与える値であり、これまで地表アルベド $\alpha_s=0.3$ に繰り込んでいた雲の反射を雲に移すことに相当する。晴天の反射（Rayleigh 散乱と地表の反射）も含めた値なので、それらを別に計算する帯域別放射では使わない（[帯域別放射](./band-radiation.md#81-実装後の積分)）。

$\braket{C}$ と惑星アルベドを地球の値に近づけるための主要な調整パラメータは $\mathrm{RH}_c$ である。

## 参考文献

- Slingo, J. M., 1987: The development and verification of a cloud prediction scheme for the ECMWF model. *Quart. J. Roy. Meteor. Soc.*, **113**, 899--927.
- Liu, Q., et al., 2021: SimCloud version 1.0: a simple diagnostic cloud scheme for idealized climate models. *Geosci. Model Dev.*, **14**, 2801--2826.
