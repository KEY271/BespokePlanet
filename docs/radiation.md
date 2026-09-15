# 日変化・季節変化を含む放射

鉛直層の添字は上端から地表へ向かって $k=1,2,\dots,N$ とし、界面気圧を $p_{k+1/2}$、層の厚さを

$$
\Delta p_k=p_{k+1/2}-p_{k-1/2}
$$

とする。大気上端は $p_{1/2}=p_T$、地表は $p_{N+1/2}=p_s$ である。

## 長波フラックス

長波放射は上向きと下向きの 2-stream で扱う。各大気層は長波を吸収し、Kirchhoff の法則に従って同じ放射率で放射する。地表から大気上端までの長波の総光学的厚さを $\tau_s=1$ とし、第 $k$ 層の光学的厚さを

$$
\Delta\tau_k=\tau_s\frac{\Delta p_k}{p_s-p_T}
$$

とする。したがって、層の長波透過率 $t_k$ は

$$
t_k=\exp(-\Delta\tau_k)
$$

である。ここでは放射が鉛直に進むものとして、拡散係数による光路長の補正は行わない。

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

とする。地表は長波について黒体とし、浅い地面層の温度を $T_s$ として下側の境界条件を

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
\left(\frac{\partial T_k}{\partial t}\right)_{\mathrm{rad}}
=\frac{g}{c_p\Delta p_k}
\left(F^{\mathrm{LW}}_{k+1/2}-F^{\mathrm{LW}}_{k-1/2}\right)
$$

で与える。ただし $c_p=1004.0\,\mathrm{J\,kg^{-1}\,K^{-1}}$ である。

## 短波放射

短波について大気は完全に透明で、吸収も散乱も行わない。太陽定数は

$$
S_0=1361\,\mathrm{W\,m^{-2}}
$$

とし、軌道離心率は 0 とする。このため太陽と惑星の距離による $S_0$ の変化は考えない。一方、地軸の傾き、惑星の自転、緯度、経度、昼夜および季節による入射角の変化は陽に計算する。

経度を $\lambda$、緯度を $\varphi$、シミュレーション開始からの経過時間を $t$ とする。Octahedral Gaussian Grid の第 $j$ 緯度リングに $M_j$ 点があるとき、各格子点では

$$
\sin\varphi_j=\mu_j,\qquad
\lambda_{i,j}=\frac{2\pi(i-1)}{M_j}
$$

を使う。このケース専用の暦では、1 太陽日を

$$
P_{\mathrm{day}}=86400\,\mathrm{s}
$$

とし、1 年を 360 太陽日、1 ヶ月を 30 太陽日とする。したがって 1 年は 12 ヶ月で、公転周期は

$$
P_{\mathrm{orb}}=360P_{\mathrm{day}}
$$

である。この 360 日暦はこの放射ケースにだけ適用し、他の計算ケースの時間設定は変更しない。円軌道上の公転黄経を

$$
L(t)=2\pi\frac{t}{P_{\mathrm{orb}}}\pmod{2\pi}
$$

と置く。地軸の傾きは $\varepsilon=23.4^\circ$ とする。赤道座標における太陽の赤経 $\alpha_\odot$ と赤緯 $\delta_\odot$ は

$$
\alpha_\odot(t)=\operatorname{atan2}
\left(\cos\varepsilon\sin L,\cos L\right),
$$

$$
\sin\delta_\odot(t)=\sin\varepsilon\sin L,\qquad
\cos\delta_\odot(t)=\sqrt{1-\sin^2\delta_\odot}
$$

から求める。$t=0$ は 1 年目の 4 月 1 日 00:00:00 とする。この時刻を北半球の春分とし、経度 $\lambda=0$ がちょうど南中するように自転位相を定める。自転角を

$$
\theta(t)=\Omega t,\qquad
\Omega=2\pi\left(\frac{1}{P_{\mathrm{day}}}+\frac{1}{360P_{\mathrm{day}}}\right)
$$

とする。自転角速度に平均公転角速度を加えることで、平均太陽に対して 1 年にちょうど 360 回自転する。

この $\Omega$ は放射ケースの惑星渦度 $f=2\Omega\sin\varphi$ と初期値の釣り合いにも用いる。他の計算ケースでは、それぞれ従来の自転角速度を変更せずに用いる。各格子点の太陽時角は

$$
H(\lambda,t)=\theta(t)+\lambda-\alpha_\odot(t)
$$

である。したがって太陽天頂角 $\zeta$ の余弦を

$$
\mu_0(\lambda,\varphi,t)=\cos\zeta
=\sin\varphi\sin\delta_\odot
+\cos\varphi\cos\delta_\odot\cos H
$$

と計算し、大気上端および地表に入る下向き短波フラックスを

$$
F^{\downarrow\mathrm{SW}}(\lambda,\varphi,t)
=S_0\max(0,\mu_0)
$$

とする。$\mu_0\leq0$ の夜側では短波フラックスは 0 である。地表の短波アルベドを $\alpha_s=0.3$ とし、地表が吸収する短波フラックスは

$$
F^{\mathrm{SW}}_{\mathrm{abs}}
=(1-\alpha_s)F^{\downarrow\mathrm{SW}}
$$

とする。短波フラックスは時間平均・日平均・東西平均を事前に取らず、放射傾向を評価する各タイムステップで、その時刻と各格子点の緯度・経度を用いて再計算する。

円軌道かつ空間的に一定な $S_0$ という仮定の下では、任意の時刻における解析的な全球平均入射量は

$$
\left\langle F^{\downarrow\mathrm{SW}}\right\rangle=\frac{S_0}{4}
$$

となる。数値実装では、Gaussian 求積による全球平均が離散化誤差の範囲でこの値に一致することをテストで確認する。ただし各格子点に与えるフラックスそのものはこの全球平均値ではなく、上記の瞬時値である。

## 地面温度

地面は、温度 $T_s$、単位面積当たり熱容量 $C_s$ の浅い上層と、温度 $T_d$、単位面積当たり熱容量 $C_d$ の深い下層に分ける。両層間の下向き熱フラックスを

$$
H_{sd}=K_{sd}(T_s-T_d)
$$

とする。ここで $K_{sd}=2\,\mathrm{W\,m^{-2}\,K^{-1}}$ は層間の熱交換係数である。

地面から大気への上向き熱フラックスを

$$
H_{sa}=\frac{p_s}{RT_N}c_pC_H\sqrt{u_N^2+v_N^2+U_g^2}(T_s-T_N)
$$

とする。ただし $C_H=10^{-3}$ で、$U_g=1\,\mathrm{m/s}$ は乱流を表す。そこで

$$
\left(\frac{\partial T_N}{\partial t}\right)_{sa}=\frac{g}{c_p\Delta p_N}H_{sa}
$$

を最下層の温度変化に加える。この寄与は重力波の陰的演算子には含めない。

浅い層と深い層の温度は

$$
C_s\frac{\partial T_s}{\partial t}
=F^{\mathrm{SW}}_{\mathrm{abs}}+F^\downarrow_{N+1/2}-\sigma T_s^4-H_{sd}-H_{sa},
$$

$$
C_d\frac{\partial T_d}{\partial t}=H_{sd}
$$

に従って時間発展させる。熱容量は

$$
C_s=2\times10^6\,\mathrm{J\,m^{-2}\,K^{-1}},\quad
C_d=2\times10^7\,\mathrm{J\,m^{-2}\,K^{-1}}
$$

重力波・超粘性の処理については $T_s,T_d$ には適用しない。ただし RAW フィルターについては適用する。

## Rayleigh 摩擦

Rayleigh 摩擦について、[Held-Suarez 強制](./Held-Suarez.md)と同じものを加える。

## 初期値

[Jablonowski–Williamson の初期状態](Jablonowski-Williamson.md)に従って速度と温度を初期化する。地表面気圧の初期値は $p_s=10^5\,\mathrm{Pa}$ とする。ただし地形については

$$
\Phi_s=0
$$

とし、平坦地形にする。また、$T_s,T_d$ は $T_N$ と等しくする。

## 時間・出力

$5$ 年間、すなわち $5\times360$ 太陽日シミュレーションする。タイムステップは $\Delta t=1200\,\mathrm{s}$ と少し長めに取る。このため 1 日は 72 ステップ、1 ヶ月は 2160 ステップ、1 年は 25920 ステップにそれぞれ厳密に一致する。

1日ごとに次の日平均量を出力する。

- 大気温度の全球平均 $\braket{T}$
- 地表温度の全球平均 $\braket{T_s}$
- 地面下層温度の全球平均 $\braket{T_d}$
- 大気の全球平均運動エネルギー
- 地表面気圧の全球平均 $\braket{p_s}$
- 全球の短波入射
- 全球の反射された短波
- 全球の出ていった長波

これらの量を求めるために別でスペクトル法を用いるとパフォーマンスが悪化するので計算結果を流用してオンラインで集計する。

日平均は各太陽日の左端を含み右端を含まない全タイムステップ（72 サンプル）を等重みで平均し、時刻の列にはその日の開始時刻を書く。

30 日ごとに次の月平均量を出力する。最初の月は 1 年目の 4 月 1 日から 4 月 30 日までである。

- 各格子点の地表温度の月平均 $T_s(\lambda,\varphi)$
- 各格子点の地表面気圧の月平均 $p_s(\lambda,\varphi)$
- 経度方向に平均した温度・流速の月平均 $[T](\varphi,\eta),[u](\varphi,\eta),[v](\varphi,\eta)$
- 上の平均からの偏差をプライムをつけて表したときの $[u'v'](\varphi,\eta),[v'T'](\varphi,\eta)$

これらもオンラインで効率的に計算する。

月平均は各 30 日区間の左端を含み右端を含まない全タイムステップ（2160 サンプル）を等重みで平均する。`m0001` は 1 年目の 4 月、`m0002` は 5 月に対応し、以後も30日ごとに連番とする。

毎年の 4/1 に次の瞬時値を出力する。

- スペクトル $\zeta,\delta,T,\ln p_s$
- 格子 $\zeta,\delta,T,u,v,\ln p_s,T_s,T_d$

初期値を含め、最終的に6つの瞬時値が出るようにする。

瞬時値のファイル名には `yearly_` を付け、開始時刻を `y0001`、以後の各 4 月 1 日を `y0002` から `y0006` として保存する。
