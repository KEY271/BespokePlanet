# 長波放射

鉛直層の添字は上端から地表へ向かって $k=1,2,\dots,N$ とし、界面気圧を $p_{k+1/2}$、層の厚さを

$$
\Delta p_k=p_{k+1/2}-p_{k-1/2}
$$

とする。大気上端は $p_{1/2}=p_T$、地表は $p_{N+1/2}=p_s$ である。

## 長波フラックス

長波放射は上向きと下向きの 2-stream で扱う。各大気層は長波を吸収し、Kirchhoff の法則に従って同じ放射率で放射する。オゾン以外の成分について、地表から大気上端までの長波の総光学的厚さを $\tau_s=1$ とし、第 $k$ 層の光学的厚さを

$$
\Delta\tau^{\mathrm{base}}_k=\tau_s\frac{\Delta p_k}{p_s-p_T}
$$

とする。さらに、[オゾン](./ozone.md)による長波光学的厚さ $\Delta\tau^{\mathrm{O_3,LW}}_k$ を加えて、層の全長波光学的厚さと透過率は

$$
\Delta\tau^{\mathrm{LW}}_k
=\Delta\tau^{\mathrm{base}}_k+\Delta\tau^{\mathrm{O_3,LW}}_k,\qquad
t_k=\exp(-\Delta\tau^{\mathrm{LW}}_k)
$$

である。オゾンを含む全光学的厚さに対して Kirchhoff の法則を適用し、層の放射率を $1-t_k$ とする。ここでは放射が鉛直に進むものとして、拡散係数による光路長の補正は行わない。

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
\left(\frac{\partial T_k}{\partial t}\right)_{\mathrm{LW}}
=\frac{g}{c_p\Delta p_k}
\left(F^{\mathrm{LW}}_{k+1/2}-F^{\mathrm{LW}}_{k-1/2}\right)
$$

で与える。ただし $c_p=1004.0\,\mathrm{J\,kg^{-1}\,K^{-1}}$ である。
