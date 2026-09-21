# 物理過程を評価する時刻

力学以外の物理過程は、現在の時刻 $X^n$ ではなく、RAW フィルター適用済みの 1 つ前の時刻 $\overline X^{n-1}$ の場で統一して評価する。対象は長波・短波放射、乾燥対流調節、顕熱交換、浅い地面層と深い地面層の熱交換、Rayleigh 摩擦、および湿潤ケースの蒸発・湿潤対流調節・大規模凝結である。たとえば[長波放射](./longwave-radiation.md)では、$T_k$ と $T_s$ を $\overline T_k^{n-1}$ と $\overline T_s^{n-1}$ に読み替え、

$$
F^\uparrow_{k-1/2}
=t_kF^\uparrow_{k+1/2}+(1-t_k)\sigma\left(\overline T_k^{n-1}\right)^4,
$$

$$
F^\downarrow_{k+1/2}
=t_kF^\downarrow_{k-1/2}+(1-t_k)\sigma\left(\overline T_k^{n-1}\right)^4,
$$

$$
F^\uparrow_{N+1/2}=\sigma\left(\overline T_s^{n-1}\right)^4
$$

とする。$X^n$ で評価しない理由は、$\sigma T^4$、乾燥対流調節、顕熱交換、Rayleigh 摩擦のような減衰項を LeapFrog の中心差分で扱うと計算モードが増幅するためである。長波の冷却率は $g/(c_p\Delta p_k)$ に比例するので、$\Delta p_k$ の小さい上層ほど 1 ステップ当たりの増幅が速い。

層の光学的厚さと加熱率に使う $\Delta p_k$ も、$(\ln p_s)^{\overline{n-1}}$ から作る。[長波放射](./longwave-radiation.md)の水蒸気による光学的厚さに使う比湿も、同じ前時刻の $\overline q_k^{+,n-1}=\max(\overline q_k^{n-1},0)$ とする。したがって透過率 $t_k$ と $\sigma T^4$ はどちらも前時刻の場から作られる。短波のオゾン光学的厚さにも同じ前時刻の界面気圧を使う。顕熱フラックス $H_{sa}$ は $\overline p_s^{n-1},\overline T_N^{n-1},\overline T_s^{n-1},\overline{\bm u}_N^{n-1}$ から、地中フラックス $H_{sd}$ は $\overline T_s^{n-1},\overline T_d^{n-1}$ から作る。Rayleigh 摩擦の $\sigma_k$ と風速、および乾燥対流調節も同じ $\overline X^{n-1}$ を使う。太陽の位置だけは予報場ではないため、傾向を加える時刻 $t^n$ から計算する。短波の雲による反射に使う雲量は、[雲](./cloud.md)のとおり $\overline X^{n-1}$ から同じステップの対流調節と大規模凝結で進めた暫定場と、同じステップの対流性降水から診断する。放射のそれ以外の入力は $\overline X^{n-1}$ のままであり、放射の傾向も $\overline X^{n-1}$ に対して加える。

放射と地面フラックスの全項を単一の前時刻場から作るので、大気・浅い地面層・深い地面層の間で交換されるエネルギーは厳密に打ち消し合う。大気上端から出る長波 $F^\uparrow_{1/2}$ も同じ $\overline X^{n-1}$ の値であり、これをそのまま診断の OLR として出力する。

$\overline T_k^{n-1}$、前時刻の界面気圧と風速は、乾燥対流調節・放射・摩擦で共用する。

最初の 2 ステップでは、[乾燥対流調節](./dry-convective-adjustment.md#初期化)と同じく $F$ の第 2 引数を 1 つ前の時刻として用いる。すなわち $X^{1/2}$ と $X^1$ を作るどちらの段階でも、すべての物理過程を $X^0$ の場で評価する。
