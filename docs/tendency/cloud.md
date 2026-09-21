# 雲

雲を予報変数としては持たず、各タイムステップで格子点ごとに気柱一つの実効的な雲量（気柱の面積のうち雲が占める割合）$C\in[0,1]$ を、相対湿度と降水から診断する。層ごとの雲量や雲の高さは区別しない。凝結した水は[大規模凝結](./large-scale-condensation.md)と[湿潤対流調節](./moist-convective-adjustment.md)のとおり直ちに降水として地表へ落ちるので、雲水・雲氷の量は追跡せず、雲は放射のためだけの診断量である。診断した雲量は[短波放射](./shortwave-radiation.md)で下向き短波を反射するためにだけ使い、[長波放射](./longwave-radiation.md)には当面影響させない。鉛直層の添字は[長波放射](./longwave-radiation.md)、full level の気圧 $p_k$ と飽和比湿 $q_s$ は[飽和比湿](./saturation-specific-humidity.md)と同じとする。

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

とする。$0\le C\le1$ である。[短波放射](./shortwave-radiation.md)はこの $C$ を使い、地表に達する前の下向き短波のうち $C\alpha_c$（$\alpha_c=0.43$）を一回だけ反射する。雲の高さは持たないので、反射する位置による加熱の違いは表れない。

## 評価する時刻と順序

雲は[物理過程を評価する時刻](./physics-time-level.md)の例外として、$\overline X^{n-1}$ そのものではなく、そこから同じステップの対流調節と大規模凝結で進めた暫定場と、同じステップの対流性降水から診断する。したがって物理過程の評価順は

1. [蒸発](./evaporation.md)
2. [乾燥対流調節](./dry-convective-adjustment.md) → [湿潤対流調節](./moist-convective-adjustment.md) → [大規模凝結](./large-scale-condensation.md)
3. 雲量の診断
4. [長波放射](./longwave-radiation.md)・[短波放射](./shortwave-radiation.md)と地表の熱収支

とする。放射は雲量以外の入力（温度、比湿、地表温度、風）を従来どおり $\overline X^{n-1}$ からとり、雲量だけを 3 の結果から受け取る。雲量は暫定場から作るが放射の傾向は $\overline X^{n-1}$ に加わるので、LeapFrog の計算モードを増幅させる減衰項にはならない。暫定場に使う時間幅は最初の 2 ステップでは $F$ の第 1 引数に読み替える。

雲量は時間方向に平滑化せず、毎ステップ診断し直す。$\mathrm{RH}_c$ 付近で雲量が滑らかに立ち上がるようにしてあるので、雲量がステップごとに 0 と 1 の間で振動することはない。

## 診断

雲量 $C(\lambda,\varphi)$ を格子で持ち、[湿潤ケース](../cases/moist.md#出力)の出力に日平均の全球平均 $\braket{C}$、月平均の $C(\lambda,\varphi)$、毎年の瞬時値の $C(\lambda,\varphi)$ を加える。瞬時値の $C$ は予報変数ではないので、瞬時値の時刻で終わるステップの傾向評価（その 1 ステップ前の場から作った暫定場）で診断した直近の値を書き、初期の瞬時値では未診断なので 0 を書く。惑星アルベド $\braket{F^{\mathrm{SW}}_{\mathrm{refl}}}/\braket{F^{\downarrow\mathrm{SW}}_{1/2}}$ は既存の短波入射と反射短波の出力から求める。

## パラメータの根拠

- $\mathrm{RH}_c=0.8$ は Slingo (1987) の低層・高層の層状雲の値で、湿潤対流調節の参照相対湿度 $0.7$ より大きくとる。
- $c_0,c_1,C^{\mathrm{conv}}_{\max}$ は Slingo (1987) の対流性の雲の値をそのまま使う。
- 雲の短波の反射率 $\alpha_c=0.43$ は[短波放射](./shortwave-radiation.md#雲による反射)で与える。地球の全球平均の雲量 $0.6$–$0.7$ に対して、雲の反射だけで惑星アルベドの $0.26$–$0.30$ を与える値であり、これまで地表アルベド $\alpha_s=0.3$ に繰り込んでいた雲の反射を雲に移すことに相当する。

5 年積分の後に、$\braket{C}$ と惑星アルベドを見て $\mathrm{RH}_c$ を見直す。$\braket{C}$ が地球の $0.6$–$0.7$ から大きく外れる場合は、まず $\mathrm{RH}_c$ を変える。

## 実装

- 雲の診断は、格子列ごとに $T^{(3)}_k,q^{(3)}_k,p_k,P_{\mathrm{conv}}$ から $C$ を返す純粋な列関数として一つのモジュールにまとめる。呼び出しは、対流調節と大規模凝結を逐次に評価している列ループの末尾で行う。暫定場はその列ループの中でしか持たないので、暫定場を格子全体で保持する配列は増やさない。
- 物理過程の作業領域に雲量 $C$ の格子を置き、放射の列ルーチンに $C$ を省略可能な引数として渡す。傾向の評価順を蒸発 → 対流調節・大規模凝結 → 放射に変える。放射が蒸発から受け取る潜熱フラックスの扱いは変えない。
- 雲の設定型（有効・無効、$\mathrm{RH}_c$、$c_0,c_1,P_0,C^{\mathrm{conv}}_{\max}$）を他の物理過程と同じくデフォルト初期化子に一か所で書く。雲アルベド $\alpha_c$ は放射の設定に置く。雲は比湿と対流性降水を必要とするので、湿潤過程が有効なケースでだけ有効にできる。
- 比湿を予報しない[放射ケース](../cases/radiation.md)と [slab ocean ケース](../cases/slab-ocean.md)では雲を診断せず、$C=0$ として短波を計算する。

## 実装後の確認項目

- 全層で $\mathrm{RH}_k\le\mathrm{RH}_c$ かつ $P_{\mathrm{conv}}=0$ の気柱で $C=0$ となり、短波のフラックスが雲を渡さない場合とビット単位で一致すること。
- ある層が飽和している気柱で $C=1$ となること。
- $P_{\mathrm{conv}}$ を増やすと $C^{\mathrm{conv}}$ が単調に増え、$0.14\,\mathrm{mm\,day^{-1}}$ 未満で 0、$85\,\mathrm{mm\,day^{-1}}$ 以上で $0.8$ になること。
- 数ステップの全球積分で $0\le\braket{C}\le1$ で、反射短波が短波入射の $\alpha_c+(1-\alpha_c)\alpha_s$ 倍を超えないこと。

## 参考文献

- Slingo, J. M., 1987: The development and verification of a cloud prediction scheme for the ECMWF model. *Quart. J. Roy. Meteor. Soc.*, **113**, 899--927.
- Liu, Q., et al., 2021: SimCloud version 1.0: a simple diagnostic cloud scheme for idealized climate models. *Geosci. Model Dev.*, **14**, 2801--2826.
