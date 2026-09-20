# Held-Suarez 強制

外力と熱を加えて Hadley 循環などが再現されるか確認するための傾向を与える。

## 熱

平衡温度を

$$
T_{\mathrm{eq},k}=\max\left[200\,\mathrm{K},\left\{315\,\mathrm{K}-60\,\mathrm{K}\sin^2\varphi-10\,\mathrm{K}\ln\frac{p_{k+1/2}\exp(-\alpha_k)}{p_0}\cos^2\varphi\right\}\left(\frac{p_{k+1/2}\exp(-\alpha_k)}{p_0}\right)^\kappa\right]
$$

と定める。$p_0=10^5\,\mathrm{hPa},\kappa=2/7$ である。そして温度の時間微分の右辺に

$$
-k_{T,k}(T_k-T_{\mathrm{eq},k})
$$

を加える。比例定数は

$$
k_{T,k}=(40\,\mathrm{day})^{-1}+[(4\,\mathrm{day})^{-1}-(40\,\mathrm{day})^{-1}]\max\left\{0,\frac{\sigma_k-\sigma_b}{1-\sigma_b}\right\}\cos^4\varphi
$$

とする。ただし

$$
\sigma_k=\frac{p_{k+1/2}\exp(-\alpha_k)}{p_s},\quad\sigma_b=0.7
$$

である。

## Rayleigh 摩擦

外力 $\bm{F}_k$ として、各格子点で

$$
\bm{F}_k=-k_{v,k}\bm{u}_k
$$

を加える。比例定数は

$$
k_{v,k}=(1\,\mathrm{day})^{-1}\max\left\{0,\frac{\sigma_k-\sigma_b}{1-\sigma_b}\right\}
$$

とする。

## 強制を評価する時刻

Newton 緩和と Rayleigh 摩擦は、RAW フィルター適用済みの 1 つ前の時刻 $\overline X^{n-1}$ の温度、風速、地表面気圧から評価する。最初の 2 ステップでは、どちらの段階でも初期場 $X^0$ を用いる。移流などの力学的な傾向は従来どおり現在時刻 $X^n$ で評価する。
