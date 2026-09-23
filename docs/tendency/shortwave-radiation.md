# 短波放射

状態: **灰色スキーム**。放射を用いるケースは[帯域別放射](./band-radiation.md)を使う。このページの灰色スキームは、設定 `radiation_config%scheme = radiation_scheme_gray` のときの放射で、比較と回帰テストのために残している。[大気上端の入射](#大気上端の入射)の太陽の位置と天頂角の計算、地表アルベドは帯域別放射でも共通である。

鉛直層の添字と界面気圧は[長波放射](./longwave-radiation.md)と同じとする。

短波は散乱を扱わず、全短波のうち割合 $f_{\mathrm{UV}}=0.02$ だけをオゾンが吸収する UV 成分として扱う。残りの $1-f_{\mathrm{UV}}$ は大気を透過する。オゾンの光学的厚さは[オゾン](./ozone.md)で与える。雲は[雲](./cloud.md)で診断した気柱の実効雲量 $C$ に基づいて、地表に達する前の下向き短波を一回だけ反射する。雲による短波の吸収と、上向き短波の大気中での再吸収・再反射は扱わない。太陽定数は

$$
S_0=1361\,\mathrm{W\,m^{-2}}
$$

とし、軌道離心率は 0 とする。このため太陽と惑星の距離による $S_0$ の変化は考えない。一方、地軸の傾き、惑星の自転、緯度、経度、昼夜および季節による入射角の変化は陽に計算する。

## 大気上端の入射

経度を $\lambda$、緯度を $\varphi$、シミュレーション開始からの経過時間を $t$ とする。Octahedral Gaussian Grid の第 $j$ 緯度リングに $M_j$ 点があるとき、各格子点では

$$
\sin\varphi_j=\mu_j,\qquad
\lambda_{i,j}=\frac{2\pi(i-1)}{M_j}
$$

を使う。[暦と軌道](../calendar.md)の太陽の赤緯 $\delta_\odot$ と時角 $H$ から、太陽天頂角 $\zeta$ の余弦を

$$
\mu_0(\lambda,\varphi,t)=\cos\zeta
=\sin\varphi\sin\delta_\odot
+\cos\varphi\cos\delta_\odot\cos H
$$

と計算し、大気上端に入る下向き短波フラックスを

$$
F^{\downarrow\mathrm{SW}}_{1/2}(\lambda,\varphi,t)
=S_0\max(0,\mu_0)
$$

とする。$\mu_0\leq0$ の夜側では、すべての短波フラックスと短波加熱を 0 とする。短波フラックスは時間平均・日平均・東西平均を事前に取らず、放射傾向を評価する各タイムステップで、その時刻と各格子点の緯度・経度を用いて再計算する。

## オゾンによる吸収

大気上端に入る短波を、UV 成分と非 UV 成分に

$$
F^{\downarrow\mathrm{UV}}_{1/2}
=f_{\mathrm{UV}}F^{\downarrow\mathrm{SW}}_{1/2},\qquad
F^{\downarrow\mathrm{nonUV}}_{1/2}
=(1-f_{\mathrm{UV}})F^{\downarrow\mathrm{SW}}_{1/2}
$$

と分ける。太陽天頂角は大気上端の入射量にだけ反映し、オゾン中の光路は常に鉛直とみなす。したがって各層の UV 透過率は

$$
t^{\mathrm{UV}}_k
=\exp\left(-\Delta\tau^{\mathrm{O_3,SW}}_k\right)
$$

であり、下向き UV フラックスを上端から地表へ

$$
F^{\downarrow\mathrm{UV}}_{k+1/2}
=t^{\mathrm{UV}}_kF^{\downarrow\mathrm{UV}}_{k-1/2}
$$

と順に計算する。非 UV 成分は全層で変化しない。第 $k$ 層が吸収する UV フラックスは

$$
A^{\mathrm{UV}}_k
=F^{\downarrow\mathrm{UV}}_{k-1/2}
-F^{\downarrow\mathrm{UV}}_{k+1/2}
$$

であり、対応する温度変化は

$$
\left(\frac{\partial T_k}{\partial t}\right)_{\mathrm{SW}}
=\frac{g}{c_p\Delta p_k}A^{\mathrm{UV}}_k
$$

とする。雲は大気を加熱しないので、短波による大気の加熱はこの UV 吸収だけである。大気に与える放射温度変化は、長波と短波の和

$$
\left(\frac{\partial T_k}{\partial t}\right)_{\mathrm{rad}}
=\left(\frac{\partial T_k}{\partial t}\right)_{\mathrm{LW}}
+\left(\frac{\partial T_k}{\partial t}\right)_{\mathrm{SW}}
$$

である。オゾン層を通過した全短波は

$$
F^{\downarrow\mathrm{SW}}_{N+1/2}
=F^{\downarrow\mathrm{nonUV}}_{1/2}
+F^{\downarrow\mathrm{UV}}_{N+1/2}
$$

である。

## 雲による反射

雲の短波アルベドを

$$
\alpha_c=0.43
$$

とし、[雲](./cloud.md)の実効雲量 $C$ を使って、オゾン層を通過した下向き短波のうち

$$
F^{\mathrm{SW}}_{\mathrm{cld}}=C\alpha_cF^{\downarrow\mathrm{SW}}_{N+1/2}
$$

を雲が反射し、残りの

$$
F^{\downarrow\mathrm{SW}}_{\mathrm{sfc}}=(1-C\alpha_c)F^{\downarrow\mathrm{SW}}_{N+1/2}
$$

が地表に達するものとする。雲は高さを持たないので、オゾンを含む全層の下、地表の直上で反射するとみなす。オゾンは $1$--$100\,\mathrm{hPa}$ にあり雲は対流圏にあるので、オゾン吸収を先に、雲の反射を後に置くこの順序は物理的な並びと一致する。雲による反射は雲量に線形で、雲量 $C=1$ の気柱の反射率がちょうど $\alpha_c$ になる。雲の光学的厚さや太陽天頂角による $\alpha_c$ の変化、層状雲と対流性の雲の区別は入れない。雲を診断しないケースでは $C=0$ である。

$\alpha_c=0.43$ は、地球の全球平均の雲量 $0.6$–$0.7$ に対して雲の反射だけで惑星アルベドの $0.26$–$0.30$ を与える値で、これまで地表アルベドに繰り込んでいた雲の反射を雲に移すことに相当する。

## 地表

地表が吸収する短波フラックスは

$$
F^{\mathrm{SW}}_{\mathrm{abs}}
=(1-\alpha_s)F^{\downarrow\mathrm{SW}}_{\mathrm{sfc}}
=(1-\alpha_s)(1-C\alpha_c)F^{\downarrow\mathrm{SW}}_{N+1/2}
$$

とする。地表で反射する短波は $\alpha_sF^{\downarrow\mathrm{SW}}_{\mathrm{sfc}}$ である。

### 地表アルベド

地表アルベドは雲を含まない地表そのものの値とし、[陸と海の混合](./land-sea-surface.md)のとおり陸面率 $f_L$ で

$$
\alpha_s=f_L\alpha_L+(1-f_L)\alpha_o,\qquad
\alpha_L=0.2,\qquad\alpha_o=0.06
$$

と混ぜる。$\alpha_o=0.06$ は開いた海面の代表値、$\alpha_L=0.2$ は雪氷のない陸面（植生と裸地）の代表値である。雪氷、海氷、天頂角による海面アルベドの増加は扱わない。全球を海で覆った[湿潤ケース](../cases/moist.md)では $\alpha_s=\alpha_o=0.06$ である。

比湿を予報せず雲を診断しない[放射ケース](../cases/radiation.md)と [slab ocean ケース](../cases/slab-ocean.md)では、雲の反射を地表に繰り込んだ

$$
\alpha_s=0.3
$$

を従来どおり使う。これは雲を含む地球の惑星アルベドの値であり、$C=0$ と組み合わせて雲のある大気の全球平均の短波収支を代用する。

## 反射短波と収支

雲と地表で反射した短波はどちらも大気中で再吸収・再反射されずに大気上端から出るものとし、反射短波を

$$
F^{\mathrm{SW}}_{\mathrm{refl}}
=F^{\mathrm{SW}}_{\mathrm{cld}}+\alpha_sF^{\downarrow\mathrm{SW}}_{\mathrm{sfc}}
=\left[C\alpha_c+(1-C\alpha_c)\alpha_s\right]F^{\downarrow\mathrm{SW}}_{N+1/2}
$$

とする。したがって $F^{\mathrm{SW}}_{\mathrm{refl}}$ は雲と地表での反射量の和と、大気上端から出る反射量の両方を表す。診断の反射短波として出力するのはこの値である。各カラムで

$$
F^{\downarrow\mathrm{SW}}_{1/2}
=\sum_{k=1}^{N}A^{\mathrm{UV}}_k
+F^{\mathrm{SW}}_{\mathrm{abs}}
+F^{\mathrm{SW}}_{\mathrm{refl}}
$$

が成り立つ。惑星アルベドは $F^{\mathrm{SW}}_{\mathrm{refl}}/F^{\downarrow\mathrm{SW}}_{1/2}$ で、オゾン吸収を無視すれば $C\alpha_c+(1-C\alpha_c)\alpha_s$ である。海上（$\alpha_s=0.06$）で $C=1$ なら $0.464$、$C=0.65$ なら $0.323$、$C=0$ なら $0.06$ になる。

円軌道かつ空間的に一定な $S_0$ という仮定の下では、任意の時刻における解析的な全球平均入射量は

$$
\left\langle F^{\downarrow\mathrm{SW}}_{1/2}\right\rangle=\frac{S_0}{4}
$$

となる。オゾン吸収と雲の反射はこの大気上端入射量を変えず、地表へ到達する短波だけを減らす。鉛直透過率は太陽天頂角に依存しないため、オゾンによる UV 吸収の解析的な全球平均は

$$
\left\langle\sum_{k=1}^{N}A^{\mathrm{UV}}_k\right\rangle
=f_{\mathrm{UV}}\frac{S_0}{4}
\left(1-e^{-\tau_{\mathrm{O_3,SW}}}\right)
=5.286599\,\mathrm{W\,m^{-2}}
$$

となり、雲量に依存しない。Gaussian 求積による大気上端入射量の全球平均は離散化誤差の範囲で $S_0/4$ に一致し、上のカラム短波エネルギー収支は雲のある気柱を含む各気柱で丸め誤差の範囲で閉じる。ただし各格子点に与えるフラックスそのものはこの全球平均値ではなく、上記の瞬時値である。

## 評価する時刻

太陽の位置は傾向を加える時刻 $t^n$ から、オゾンの光学的厚さに使う界面気圧は[物理過程を評価する時刻](./physics-time-level.md)のとおり $(\ln p_s)^{\overline{n-1}}$ から作る。雲量 $C$ は[雲](./cloud.md)のとおり、同じステップの対流調節と大規模凝結で進めた暫定場と対流性降水から診断した値を使う。地表アルベドは時間変化しない。

## 数値的な注意

- 雲の反射は $F^{\downarrow\mathrm{SW}}_{N+1/2}$ に係数 $(1-C\alpha_c)$ を掛けるだけで、オゾン吸収の計算と大気の加熱は雲がないときと変わらない。同じ入射に対して $C$ を増やすと、反射短波は単調に増え、地表の吸収は単調に減り、大気の UV 加熱は変わらない。$C=0$ のとき、地表に達する短波、地表の吸収、反射短波は雲を入れる前の式にビット単位で一致する。したがって雲を診断しない放射ケースと slab ocean ケースの結果は変わらない。
- 雲と地表で 2 段に反射するが、上向きの短波は扱わないので、雲の底と地表の間の多重反射は含まない。地表アルベドが小さいのでこの寄与は $C\alpha_c\alpha_s(1-C\alpha_c)\lesssim1.5\,\%$ である。
