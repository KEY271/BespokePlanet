# Slab ocean

海洋を水平流・鉛直流を持たない一層の水柱として扱い、海面が受け取る[短波](./shortwave-radiation.md)と[長波](./longwave-radiation.md)、および大気との顕熱交換によって海面温度を時間発展させる。各格子点の水柱は互いに独立であり、海流、水平熱輸送、鉛直混合、海氷、蒸発および潜熱フラックスは扱わない。

## 海洋温度

slab ocean は深さ

$$
h_o=30\,\mathrm{m}
$$

の一様な水温 $T_o$ を持つ。海水密度を $\rho_o=1000\,\mathrm{kg\,m^{-3}}$、海水の比熱を $c_o=4186\,\mathrm{J\,kg^{-1}\,K^{-1}}$ とし、単位面積当たりの熱容量を

$$
C_o=\rho_oc_oh_o
=1.2558\times10^8\,\mathrm{J\,m^{-2}\,K^{-1}}
$$

とする。

海面から大気への上向き顕熱フラックスは[地面](./ground.md)と同じく

$$
H_{oa}=\frac{p_s}{RT_N}c_pC_H
\sqrt{u_N^2+v_N^2+U_g^2}(T_o-T_N)
$$

とする。ここで $C_H=10^{-3}$、$U_g=1\,\mathrm{m\,s^{-1}}$ である。大気最下層には

$$
\left(\frac{\partial T_N}{\partial t}\right)_{oa}
=\frac{g}{c_p\Delta p_N}H_{oa}
$$

を加える。この寄与は重力波の陰的演算子には含めない。

海洋温度は

$$
C_o\frac{\partial T_o}{\partial t}
=F^{\mathrm{SW}}_{\mathrm{abs}}
+F^\downarrow_{N+1/2}
-\sigma T_o^4
-H_{oa}
$$

に従う。海面が大気へ与える長波放射と顕熱には、海洋から差し引くものと同じフラックスを使う。したがって大気と slab ocean を合わせたエネルギー収支は閉じる。海底との熱交換や prescribed Q-flux は加えない。

右辺の放射フラックスと顕熱フラックスは、すべて[物理過程を評価する時刻](./physics-time-level.md)のとおり $\overline X^{n-1}$ の場で評価する。$T_o$ には重力波と超粘性を適用せず、RAW フィルターは適用する。
