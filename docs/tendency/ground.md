# 地面

地表が受け取る[短波](./shortwave-radiation.md)と[長波](./longwave-radiation.md)のフラックスを境界条件として、地面の温度を時間発展させる。

陸海タイル経路では以下の $T_s$ を陸温度 $T_L$ と読み、熱容量とフラックスを陸面積当たりで扱う。大気との交換と格子全体の蓄熱に陸面率を掛ける。海水と海氷は[別のタイル](./land-sea-surface.md)で進める。積雪があるケースでは顕熱 $H_{sa}$ を $1-f$ 倍にし、浅い層の更新後に $275\,\mathrm K$ を超えた分の熱を融雪に使う（[降雪と積雪](./snow.md#53-積雪と融雪の更新)）。

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

とする。ただし $C_H=10^{-3}$ で、$U_g=1\,\mathrm{m/s}$ は乱流を表す。$p_N=p_{N+1/2}\exp(-\alpha_N)$ は最下層 full level の気圧、$\kappa=2/7$ である。密度 $\rho_N$ も最下層の気圧と温度から作り、水蒸気を含むケースでも乾燥空気の式のままとする。括弧内は温位の差 $(\theta_s-\theta_N)(p_s/p_0)^\kappa$ に等しく、地面と最下層の温位が等しい乾燥断熱中立な大気では $H_{sa}=0$ になる。最下層 full level は $p_s=1000\,\mathrm{hPa}$ のとき約 $939\,\mathrm{hPa}$ にあり、中立な大気でも $T_s-T_N\approx5\,\mathrm{K}$ となるので、温度差 $T_s-T_N$ をそのまま使うと中立成層でも上向きの顕熱フラックスが出てしまう。この $H_{sa}$ による最下層の温度変化

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

に従って時間発展させる。蒸発を有効にしたケースでは、[蒸発](./evaporation.md)の潜熱フラックス $LE$ も $T_s$ の式の右辺から差し引く。右辺の放射フラックス、$H_{sd}$、$H_{sa}$ はすべて[物理過程を評価する時刻](./physics-time-level.md)のとおり $\overline X^{n-1}$ の場で評価する。大気層に加える顕熱・放射フラックスと地面層から差し引く値には同じ前時刻のフラックスを使うので、大気・地表・地中を合わせたエネルギー収支は閉じたままである。熱容量は

$$
C_s=2\times10^6\,\mathrm{J\,m^{-2}\,K^{-1}},\quad
C_d=2\times10^7\,\mathrm{J\,m^{-2}\,K^{-1}}
$$

$T_s,T_d$ は他の予報変数と同じ LeapFrog で

$$
T_s^{n+1}=\overline T_s^{n-1}+2\Delta t\left(\frac{\partial T_s}{\partial t}\right)\left[\overline X^{n-1}\right],\qquad
T_d^{n+1}=\overline T_d^{n-1}+2\Delta t\left(\frac{\partial T_d}{\partial t}\right)\left[\overline X^{n-1}\right]
$$

と進め、最初の 2 ステップも[乾燥大気](../dynamics/dry.md#初期化)の $F$ と同じ時間幅・同じ引数で扱う。$T_s,T_d$ は水平移流を持たないので、大気の予報変数と違ってスペクトルでは持たず、格子点ごとの値として保持する（[陸と海の混合](./land-sea-surface.md#t_s-を格子で持つ理由)）。スペクトル変換・重力波・超粘性の処理は適用せず、RAW フィルターは格子点ごとに適用する。
