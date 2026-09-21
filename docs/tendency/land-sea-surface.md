# 陸と海の混合

[地面](./ground.md)の二層モデルと [slab ocean](./slab-ocean.md) を、格子点ごとに与えた陸面率 $f_L(\lambda,\varphi)$ で一つの地表収支に混ぜる。$f_L$ は[地形](../dynamics/topography.md)の生成で決まる時間変化しない格子量で、$f_L=1$ が陸、$f_L=0$ が海、$0<f_L<1$ が海岸である。放射・顕熱・蒸発のバルク式は変えず、地表の熱容量、アルベド、湿り具合、地中との熱交換だけを $f_L$ で線形に混ぜる。

## 混ぜ方

各格子点の地表は温度 $T_s$ を一つだけ持ち、その単位面積当たりの熱容量を

$$
C=f_LC_s+(1-f_L)C_o
$$

とする。$C_s=2\times10^6\,\mathrm{J\,m^{-2}\,K^{-1}}$ は[地面](./ground.md)の浅い層、$C_o=\rho_oc_oh_o=1.2558\times10^8\,\mathrm{J\,m^{-2}\,K^{-1}}$ は深さ $30\,\mathrm{m}$ の [slab ocean](./slab-ocean.md) の値である。同様に短波アルベドと[蒸発](./evaporation.md)の湿り具合を

$$
\alpha=f_L\alpha_L+(1-f_L)\alpha_o,\qquad
\beta=f_L\beta_L+(1-f_L)\beta_o
$$

とし、地中との熱交換係数を

$$
K=f_LK_{sd}
$$

とする。$K_{sd}=2\,\mathrm{W\,m^{-2}\,K^{-1}}$ は地面の値である。既定値は

$$
\alpha_L=0.2,\qquad\alpha_o=0.06,\qquad\beta_L=0.5,\qquad\beta_o=1
$$

とする。アルベドは雲を含まない地表そのものの値で、雲の反射は[雲](./cloud.md)で診断した雲量から[短波放射](./shortwave-radiation.md)が別に計算する。$\alpha_o=0.06$ は開いた海面の代表値、$\alpha_L=0.2$ は雪氷のない陸面の代表値である。比湿を予報せず雲を診断しない[放射ケース](../cases/radiation.md)と [slab ocean ケース](../cases/slab-ocean.md)では、雲の反射を地表に繰り込んだ従来の $\alpha_L=\alpha_o=0.3$ を使う。$\beta_L$ は土壌水分を予報しない定数で、陸を海と区別する主要なパラメータである。雪氷、海氷、河川や流出は扱わない。

## 地表収支

地表と地中の温度は

$$
C\frac{\partial T_s}{\partial t}
=F^{\mathrm{SW}}_{\mathrm{abs}}+F^\downarrow_{N+1/2}-\sigma T_s^4-H_{sd}-H_{sa}-LE,\qquad
H_{sd}=K(T_s-T_d),
$$

$$
C_d\frac{\partial T_d}{\partial t}=H_{sd}
$$

に従う。$F^{\mathrm{SW}}_{\mathrm{abs}}=(1-\alpha)F^{\downarrow\mathrm{SW}}_{N+1/2}$ で、$F^{\downarrow\mathrm{SW}}_{N+1/2}$ は[短波放射](./shortwave-radiation.md)の雲による反射を差し引いた後に地表に達する短波（同文書の $F^{\downarrow\mathrm{SW}}_{\mathrm{sfc}}$）と読む。反射短波は雲の反射と地表の反射 $\alpha F^{\downarrow\mathrm{SW}}_{\mathrm{sfc}}$ の和である。顕熱フラックス $H_{sa}$ は[地面](./ground.md)の式のまま、蒸発 $E$ は[蒸発](./evaporation.md)の式で $\beta$ に上の混ぜた値を使う。大気最下層に加える顕熱・潜熱・長波と、地表から差し引く値には同じ前時刻のフラックスを使うので、大気・地表・地中を合わせたエネルギー収支は $f_L$ の値に関わらず閉じる。

深い層の熱容量 $C_d=2\times10^7\,\mathrm{J\,m^{-2}\,K^{-1}}$ は $f_L$ で割らない。海（$f_L=0$）では $K=0$ となって $T_d$ は初期値に留まり、収支に現れない。海岸（$0<f_L<1$）では単位面積当たり $C_d$ の深い層が $f_LK_{sd}$ で浅い層とつながることになり、陸の部分だけで見た深い層の熱容量は $C_d/f_L$ に相当する。これは物理的には粗いが、エネルギーは保存し、$f_L\to0$ で特異にならない。

### 極限での一致

$f_L\equiv1$ では $C=C_s$、$\alpha=\alpha_L$、$\beta=\beta_L$、$K=K_{sd}$ となって[地面](./ground.md)の二層モデルに、$f_L\equiv0$ では $C=C_o$、$\alpha=\alpha_o$、$\beta=1$、$K=0$ となって [slab ocean](./slab-ocean.md) にそれぞれ戻る。$1\cdot C_s+0\cdot C_o$ と $0\cdot C_s+1\cdot C_o$ は浮動小数点でも厳密に $C_s$、$C_o$ に等しいので、[放射ケース](../cases/radiation.md)を $f_L\equiv1$、[slab ocean ケース](../cases/slab-ocean.md)と[湿潤ケース](../cases/moist.md)を $f_L\equiv0$ として走らせても、それぞれのケースのアルベドの値を使う限り結果はビット単位で変わらない。既存ケースはこの混ぜ方を通して計算し、この一致を確認項目にする。

## 一枚の地表として扱う理由

海岸の格子点で陸と海の温度を別々に持つ二枚のタイルにすれば、それぞれの熱容量で温度を進めて、フラックスを $f_L$ で面積平均できる。しかし $T_s$ は他の予報変数と同じくスペクトルで持ち、RAW フィルターを掛けているので、タイルを増やすと状態変数と出力が二重になる。まず $T_s$ を一つにして混ぜた熱容量で進め、海岸の幅を[地形](../dynamics/topography.md#解像度との関係)の条件で滑らかにとることで対応する。

## 数値上の注意

$C_s$ と $C_o$ は約 $60$ 倍違うので、$\partial T_s/\partial t$ は海岸を横切って急に変わる。$T_s$ の傾向は格子で計算してからスペクトルに射影するので、海岸が切断波数に対して急すぎると $T_s$ に Gibbs のリップルが出る。$T_s$ には超粘性を掛けないので、リップルは RAW フィルター以外では減衰しない。そこで

- 陸面率の海岸の幅は $w\ge180^\circ/T$ とし、陸面率そのものは切断せず解析形のまま格子で使う。
- 月平均の $T_s(\lambda,\varphi)$ を見て、海岸の海側に陸の日変化が漏れていないこと、$T_s$ が海岸で単調に変わることを確認する。

陸では $C_s$ が小さいので日変化の振幅が大きくなる。正味フラックスの振幅を $300\,\mathrm{W\,m^{-2}}$ とすると $T_s$ の日変化は $300/(C_s\cdot2\pi/86400)\approx2\,\mathrm{K}$ の桁で、長波と顕熱の負のフィードバックでさらに抑えられる。タイムステップ $\Delta t=1200\,\mathrm{s}$ に対して $C_s/(4\sigma T_s^3)\approx4\times10^5\,\mathrm{s}$ なので、前時刻で評価する陽的な長波冷却は安定である。

## 評価する時刻

$C,\alpha,\beta,K$ は時間変化しないので、フラックスの評価時刻は変わらない。すべて[物理過程を評価する時刻](./physics-time-level.md)のとおり $\overline X^{n-1}$ の場で評価し、$T_s,T_d$ は[地面](./ground.md)と同じ LeapFrog で進め、超粘性と重力波の処理は適用せず、RAW フィルターは適用する。

## 診断

陸と海で分けた全球平均を診断に加える。面積平均を $\braket{\cdot}$ として

$$
\braket{X}_L=\frac{\braket{f_LX}}{\braket{f_L}},\qquad
\braket{X}_O=\frac{\braket{(1-f_L)X}}{\braket{1-f_L}}
$$

と定義し、$T_s$、降水 $P$、蒸発 $E$ について $\braket{\cdot}_L$ と $\braket{\cdot}_O$ をオンラインで集計する。$\braket{f_L}=0$ または $1$ の場合は分母が $0$ になる側を出力しない。
