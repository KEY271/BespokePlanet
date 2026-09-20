# 地面

地表が受け取る[短波](./shortwave-radiation.md)と[長波](./longwave-radiation.md)のフラックスを境界条件として、地面の温度を時間発展させる。

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

に従って時間発展させる。右辺の放射フラックス、$H_{sd}$、$H_{sa}$ はすべて[物理過程を評価する時刻](./physics-time-level.md)のとおり $\overline X^{n-1}$ の場で評価する。大気層に加える顕熱・放射フラックスと地面層から差し引く値には同じ前時刻のフラックスを使うので、大気・地表・地中を合わせたエネルギー収支は閉じたままである。熱容量は

$$
C_s=2\times10^6\,\mathrm{J\,m^{-2}\,K^{-1}},\quad
C_d=2\times10^7\,\mathrm{J\,m^{-2}\,K^{-1}}
$$

重力波・超粘性の処理については $T_s,T_d$ には適用しない。ただし RAW フィルターについては適用する。
