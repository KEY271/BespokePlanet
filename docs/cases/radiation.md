# 放射ケース

日変化・季節変化を含む放射を有効にした計算ケースの定義である。力学の傾向項は [乾燥大気](../dynamics/dry.md)のとおりに解く。

## 物理過程

力学に加えて、[長波放射](../tendency/longwave-radiation.md)、[短波放射](../tendency/shortwave-radiation.md)、[オゾン](../tendency/ozone.md)、[地面](../tendency/ground.md)、[上層の Rayleigh 摩擦](../tendency/upper-rayleigh-friction.md)を有効にする。比湿は予報しないので、長波の水蒸気による光学的厚さには固定した[基準水蒸気分布](../tendency/longwave-radiation.md#基準水蒸気分布) $\overline q^{\mathrm{ref}}_k$（$q_0=0.010$、可降水量約 $25\,\mathrm{kg\,m^{-2}}$）を使う。同じ理由で[雲](../tendency/cloud.md)は診断せず、地表アルベド $\alpha_s=0.3$ に雲の反射を繰り込む。[Held–Suarez 強制](../tendency/Held-Suarez.md)の Newton 緩和は使わない。これらはすべて[物理過程を評価する時刻](../tendency/physics-time-level.md)のとおり $\overline X^{n-1}$ の場で評価する。[暦と軌道](../calendar.md)の 360 日暦と、そこから定まる自転角速度 $\Omega$ は、この放射ケースにだけ適用する。他のケースの時間設定と自転角速度は変更しない。

## 乾燥対流調節

放射によって対流不安定になりがちなので、[乾燥対流調節](../tendency/dry-convective-adjustment.md)を行う。

## 超粘性

放射ケースでも、相対渦度、発散、温度には [乾燥大気](../dynamics/dry.md#超粘性)と同じ超粘性を全層で適用する。上層だけラプラシアンの次数を変えるスポンジ型超粘性は用いず、上層発散の抑制は[上層の Rayleigh 摩擦](../tendency/upper-rayleigh-friction.md)で行う。

## 初期値

[Jablonowski–Williamson の初期状態](../dynamics/Jablonowski-Williamson.md)（以下 JW）をもとに初期化する。温度と地表面気圧 $p_s=p_0=10^5\,\mathrm{Pa}$ は JW と同じとし、$\Omega$ には放射ケースの値を用いる。ただし地形については

$$
\Phi_s=0
$$

とし、平坦地形にする。JW の東西風は $\Phi_s\neq0$ を前提に釣り合っているので、そのままでは $\Phi_s=0$ と整合しない。そこで温度はそのままにして、東西風の方を $\Phi_s=0$ に釣り合うよう調整する。

### 釣り合った東西風

JW の $\eta_v,A(\varphi),B(\varphi)$ を用い、

$$
c(\eta)=\cos^{3/2}\eta_v,\qquad
c_s=\cos^{3/2}\eta_{v,s},\qquad
\eta_{v,s}=\frac{\pi}{2}(1-\eta_0)
$$

と置く。JW の基本場の東西風と、その地表 $\eta=1$ での値を

$$
u_{\mathrm{JW}}(\varphi,\eta)=u_0c(\eta)\sin^2(2\varphi),\qquad
u_{\mathrm{JW},s}(\varphi)=u_0c_s\sin^2(2\varphi)
$$

とする。JW のジオポテンシャルの緯度依存部分を

$$
G(\varphi,\eta)=u_0c(\eta)\left[u_0c(\eta)A(\varphi)+a\Omega B(\varphi)\right]
$$

と書くと、JW の基本場は $\Phi=\overline\Phi(\eta)+G(\varphi,\eta)$ であり、$\Phi_s=G(\varphi,1)$ である。温度を変えずに $\Phi_s=0$ とすると、静水圧の関係から全層のジオポテンシャルが一様に $G(\varphi,1)$ だけずれ、

$$
\Phi(\varphi,\eta)=\overline\Phi(\eta)+G(\varphi,\eta)-G(\varphi,1)
$$

となる。$p_s$ が一様なので $\eta$ 面は等圧面であり、$v=0$ の東西対称な流れが定常であるための条件は、曲率項を含めた傾度風バランス

$$
\left(2\Omega\sin\varphi+\frac{u\tan\varphi}{a}\right)u
=-\frac{1}{a}\frac{\partial\Phi}{\partial\varphi}
$$

である。JW の東西風は $G(\varphi,\eta)$ に対してこの関係を厳密に満たすので、

$$
F(u)=\left(2\Omega\sin\varphi+\frac{u\tan\varphi}{a}\right)u
$$

と置けば、求める東西風 $u_{\mathrm{bal}}$ は

$$
F(u_{\mathrm{bal}})=F(u_{\mathrm{JW}})-F(u_{\mathrm{JW},s})
$$

を満たす。$a/\tan\varphi$ を掛けると $u$ の 2 次方程式

$$
\left(u_{\mathrm{bal}}+a\Omega\cos\varphi\right)^2
=\left(a\Omega\cos\varphi\right)^2+D,\\
\begin{aligned}
D(\varphi,\eta)&=\left(u_{\mathrm{JW}}-u_{\mathrm{JW},s}\right)
\left(u_{\mathrm{JW}}+u_{\mathrm{JW},s}+2a\Omega\cos\varphi\right)\\
&=u_0^2\left[c(\eta)-c_s\right]\sin^2(2\varphi)
\left\{\left[c(\eta)+c_s\right]\sin^2(2\varphi)+\frac{2a\Omega}{u_0}\cos\varphi\right\}
\end{aligned}
$$

になる。$\Phi_s=0$ のときに風速 0 へ連続につながる根を選び、極付近での桁落ちを避けるため有理化して

$$
u_{\mathrm{bal}}(\varphi,\eta)
=\frac{D(\varphi,\eta)}
{a\Omega\cos\varphi+\sqrt{\left(a\Omega\cos\varphi\right)^2+D(\varphi,\eta)}}
$$

とする。$|\eta_v|$ は $0\le\eta\le1$ で $\eta=1$ のとき最大なので $c(\eta)\ge c_s$ であり、$D\ge0$ となって根号の中は常に正である。また $\sin^2(2\varphi)$ の因子により、$u_{\mathrm{bal}}$ は赤道と極で 0 になる。

$u_{\mathrm{bal}}$ は地表 $\eta=1$ で 0 となり、上空ほど JW の東西風に近づく。$N=12$ の full level では、最大値はおおむね緯度 $45^\circ$ にあり、最上層 $\eta=0.002$ で約 $23.2\,\mathrm{m\,s^{-1}}$（JW では約 $31.1$）、$\eta=0.23$ で約 $27.2\,\mathrm{m\,s^{-1}}$（約 $35.0$）、最下層 $\eta=0.94$ で約 $3.0\,\mathrm{m\,s^{-1}}$（約 $11.3$）となる。

東西風には、局所的な擾乱を加えず

$$
u=u_\mathrm{bal},v=0
$$

を初期値とし、$\zeta,\delta$ をスペクトル法で求める。温度は JW の式をそのまま用いる。また、$T_s,T_d$ は $T_N$ と等しくする。

## 時間

[暦と軌道](../calendar.md)の暦に従い、$5$ 年間、すなわち $5\times360$ 太陽日シミュレーションする。タイムステップは $\Delta t=1200\,\mathrm{s}$ とする。このため 1 日は 72 ステップ、1 ヶ月は 2160 ステップ、1 年は 25920 ステップにそれぞれ厳密に一致する。

## 出力

1日ごとに次の日平均量を出力する。

- 大気温度の全球平均 $\braket{T}$
- 地表温度の全球平均 $\braket{T_s}$
- 地面下層温度の全球平均 $\braket{T_d}$
- 大気の全球平均運動エネルギー
- 地表面気圧の全球平均 $\braket{p_s}$
- 大気上端における全球の短波入射
- 大気上端から出る全球の反射短波
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
