# Held-Suarez 強制

外力と熱を加えて Hadley 循環などが再現されるか確認する。

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

## 初期状態

初期状態は $u=v=0$ で静止した等温大気とする。ただし対称性を破るために小さな擾乱を加え

$$
T=264\,\mathrm{K}+0.1\,\mathrm{K}\cos\varphi\sin\lambda
$$

とする。圧力は $p_s=10^5\,\mathrm{Pa}$ とし、ジオポテンシャルは $\Phi_s=0$ の平坦な地形を与える。

## 時間

計算がうまく行っているか確認するため、とりあえず200日間で試す。出力は5日ごとに保存する。タイムステップは $\Delta t=900\,\mathrm{s}$ で他の計算と同じとする。
