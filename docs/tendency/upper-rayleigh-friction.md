# 上層の Rayleigh 摩擦

上端の2層（$k=1,2$）に Rayleigh 摩擦を加え、放射によって励起される上層の発散を減衰させる。

外力を

$$
\bm{F}_k=-k_{v,k}^{(\mathrm{top})}\overline{\bm{u}}_k^{n-1}
$$

とし、比例定数を

$$
k_{v,k}^{(\mathrm{top})}
=(1\,\mathrm{day})^{-1}I_{\mathrm{top},k},\\
I_{\mathrm{top},k}=
\begin{cases}
1, & k=1,2,\\
0, & k=3,\dots,12
\end{cases}
$$

とする。

上端 2 層はそれぞれ $1$-$3\,\mathrm{hPa}$ と $3$-$10\,\mathrm{hPa}$ に対応し、風速の減衰時定数は $1\,\mathrm{day}$ である。上層側の追加係数は各層内で一定なので、上端 2 層の相対渦度と発散には前時刻の値からそれぞれ $-\overline\zeta_k^{n-1}/(1\,\mathrm{day})$ と $-\overline\delta_k^{n-1}/(1\,\mathrm{day})$ が加わる。
