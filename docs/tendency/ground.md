# 地面

地表が受け取る[短波](./shortwave-radiation.md)と[長波](./longwave-radiation.md)のフラックスを境界条件として、地面の温度を時間発展させる。

## 地面温度

地面は、温度 $T_s$、単位面積当たり熱容量 $C_s$ の浅い上層と、温度 $T_d$、単位面積当たり熱容量 $C_d$ の深い下層に分ける。両層間の下向き熱フラックスを

$$
H_{sd}=K_{sd}(T_s-T_d)
$$

とする。ここで $K_{sd}=2\,\mathrm{W\,m^{-2}\,K^{-1}}$ は層間の熱交換係数である。

地面から大気への上向き熱フラックスを、最下層の空気を地表面気圧まで乾燥断熱的に外挿した温度と地面温度の差から

$$
H_{sa}=\rho_Nc_pC_H\sqrt{u_N^2+v_N^2+U_g^2}\left[T_s-T_N\left(\frac{p_s}{p_N}\right)^\kappa\right],\qquad
\rho_N=\frac{p_N}{RT_N}
$$

とする。ただし $C_H=10^{-3}$ で、$U_g=1\,\mathrm{m/s}$ は乱流を表す。$p_N=p_{N+1/2}\exp(-\alpha_N)$ は最下層 full level の気圧、$\kappa=2/7$ である。括弧内は温位の差 $(\theta_s-\theta_N)(p_s/p_0)^\kappa$ に等しく、地面と最下層の温位が等しい乾燥断熱中立な大気では $H_{sa}=0$ になる。最下層 full level は $p_s=1000\,\mathrm{hPa}$ のとき約 $939\,\mathrm{hPa}$ にあり、中立な大気でも $T_s-T_N\approx5\,\mathrm{K}$ となるので、温度差 $T_s-T_N$ をそのまま使うと中立成層でも上向きの顕熱フラックスが出てしまう。密度 $\rho_N$ も最下層の気圧と温度から作り、水蒸気を含むケースでも乾燥空気の式のままとする。そこで

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

$T_s,T_d$ は他の予報変数と同じ LeapFrog で

$$
T_s^{n+1}=\overline T_s^{n-1}+2\Delta t\left(\frac{\partial T_s}{\partial t}\right)\left[\overline X^{n-1}\right],\qquad
T_d^{n+1}=\overline T_d^{n-1}+2\Delta t\left(\frac{\partial T_d}{\partial t}\right)\left[\overline X^{n-1}\right]
$$

と進め、最初の 2 ステップも[乾燥大気](../dynamics/dry.md#初期化)の $F$ と同じ時間幅・同じ引数で扱う。重力波・超粘性の処理については $T_s,T_d$ には適用しない。ただし RAW フィルターについては適用する。
