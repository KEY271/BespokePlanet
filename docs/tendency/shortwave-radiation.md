# 短波放射

鉛直層の添字と界面気圧は[長波放射](./longwave-radiation.md)と同じとする。

短波は散乱を扱わず、全短波のうち割合 $f_{\mathrm{UV}}=0.02$ だけをオゾンが吸収する UV 成分として扱う。残りの $1-f_{\mathrm{UV}}$ は大気を透過する。オゾンの光学的厚さは[オゾン](./ozone.md)で与える。太陽定数は

$$
S_0=1361\,\mathrm{W\,m^{-2}}
$$

とし、軌道離心率は 0 とする。このため太陽と惑星の距離による $S_0$ の変化は考えない。一方、地軸の傾き、惑星の自転、緯度、経度、昼夜および季節による入射角の変化は陽に計算する。

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

とする。$\mu_0\leq0$ の夜側では、すべての短波フラックスと短波加熱を 0 とする。

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

とする。大気に与える放射温度変化は、長波と短波の和

$$
\left(\frac{\partial T_k}{\partial t}\right)_{\mathrm{rad}}
=\left(\frac{\partial T_k}{\partial t}\right)_{\mathrm{LW}}
+\left(\frac{\partial T_k}{\partial t}\right)_{\mathrm{SW}}
$$

である。

地表に達する全短波は

$$
F^{\downarrow\mathrm{SW}}_{N+1/2}
=F^{\downarrow\mathrm{nonUV}}_{1/2}
+F^{\downarrow\mathrm{UV}}_{N+1/2}
$$

である。

地表の短波アルベドを $\alpha_s=0.3$ とし、地表が吸収する短波フラックスは

$$
F^{\mathrm{SW}}_{\mathrm{abs}}
=(1-\alpha_s)F^{\downarrow\mathrm{SW}}_{N+1/2}
$$

とする。反射短波は

$$
F^{\mathrm{SW}}_{\mathrm{refl}}
=\alpha_sF^{\downarrow\mathrm{SW}}_{N+1/2}
$$

とし、上向きの反射短波が大気中で再吸収される効果は扱わない。したがって $F^{\mathrm{SW}}_{\mathrm{refl}}$ は地表での反射量と大気上端から出る反射量の両方を表す。短波フラックスは時間平均・日平均・東西平均を事前に取らず、放射傾向を評価する各タイムステップで、その時刻と各格子点の緯度・経度を用いて再計算する。

円軌道かつ空間的に一定な $S_0$ という仮定の下では、任意の時刻における解析的な全球平均入射量は

$$
\left\langle F^{\downarrow\mathrm{SW}}_{1/2}\right\rangle=\frac{S_0}{4}
$$

となる。オゾン吸収はこの大気上端入射量を変えず、地表へ到達する短波だけを減らす。各カラムでは

$$
F^{\downarrow\mathrm{SW}}_{1/2}
=\sum_{k=1}^{N}A^{\mathrm{UV}}_k
+F^{\mathrm{SW}}_{\mathrm{abs}}
+F^{\mathrm{SW}}_{\mathrm{refl}}
$$

が成り立つ。鉛直透過率は太陽天頂角に依存しないため、オゾンによる UV 吸収の解析的な全球平均は

$$
\left\langle\sum_{k=1}^{N}A^{\mathrm{UV}}_k\right\rangle
=f_{\mathrm{UV}}\frac{S_0}{4}
\left(1-e^{-\tau_{\mathrm{O_3,SW}}}\right)
=5.286599\,\mathrm{W\,m^{-2}}
$$

となる。数値実装では、Gaussian 求積による大気上端入射量の全球平均が離散化誤差の範囲で $S_0/4$ に一致すること、および上のカラム短波エネルギー収支をテストで確認する。ただし各格子点に与えるフラックスそのものはこの全球平均値ではなく、上記の瞬時値である。
