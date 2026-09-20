# 暦と軌道

経度を $\lambda$、シミュレーション開始からの経過時間を $t$ とする。

暦では、1 太陽日を

$$
P_{\mathrm{day}}=86400\,\mathrm{s}
$$

とし、1 年を 360 太陽日、1 ヶ月を 30 太陽日とする。したがって 1 年は 12 ヶ月で、公転周期は

$$
P_{\mathrm{orb}}=360P_{\mathrm{day}}
$$

である。円軌道上の公転黄経を

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

この $\Omega$ は惑星渦度 $f=2\Omega\sin\varphi$ と初期値の釣り合いにも用いる。各格子点の太陽時角は

$$
H(\lambda,t)=\theta(t)+\lambda-\alpha_\odot(t)
$$

である。
