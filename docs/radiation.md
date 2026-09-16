# 日変化・季節変化を含む放射

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

とする。さらに、後述するオゾン分布の規格化された層割合を $q_k$ とし、オゾンによる長波光学的厚さを

$$
\Delta\tau^{\mathrm{O_3,LW}}_k
=\tau_{\mathrm{O_3,LW}}q_k,\qquad
\tau_{\mathrm{O_3,LW}}=0.005
$$

とする。層の全長波光学的厚さと透過率は

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

### 物理過程を評価する時刻

放射ケースの力学以外の物理過程は、現在の時刻 $X^n$ ではなく、RAW フィルター適用済みの 1 つ前の時刻 $\overline X^{n-1}$ の場で統一して評価する。対象は長波・短波放射、乾燥対流調節、顕熱交換、浅い地面層と深い地面層の熱交換、および Rayleigh 摩擦である。たとえば長波では、上の $T_k$ と $T_s$ を $\overline T_k^{n-1}$ と $\overline T_s^{n-1}$ に読み替え、

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

層の光学的厚さと加熱率に使う $\Delta p_k$ も、$(\ln p_s)^{\overline{n-1}}$ から作る。短波のオゾン光学的厚さにも同じ前時刻の界面気圧を使う。顕熱フラックス $H_{sa}$ は $\overline p_s^{n-1},\overline T_N^{n-1},\overline T_s^{n-1},\overline{\bm u}_N^{n-1}$ から、地中フラックス $H_{sd}$ は $\overline T_s^{n-1},\overline T_d^{n-1}$ から作る。Rayleigh 摩擦の $\sigma_k$ と風速、および乾燥対流調節も同じ $\overline X^{n-1}$ を使う。太陽の位置だけは予報場ではないため、傾向を加える時刻 $t^n$ から計算する。

放射と地面フラックスの全項を単一の前時刻場から作るので、大気・浅い地面層・深い地面層の間で交換されるエネルギーは厳密に打ち消し合う。大気上端から出る長波 $F^\uparrow_{1/2}$ も同じ $\overline X^{n-1}$ の値であり、これをそのまま診断の OLR として出力する。

$\overline T_k^{n-1}$、前時刻の界面気圧と風速は、乾燥対流調節・放射・摩擦で共用する。

最初の 2 ステップでは、[乾燥対流調節](./dry-convective-adjustment.md#初期化)と同じく $F$ の第 2 引数を 1 つ前の時刻として用いる。すなわち $X^{1/2}$ と $X^1$ を作るどちらの段階でも、すべての物理過程を $X^0$ の場で評価する。

## 短波放射

短波は散乱を扱わず、全短波のうち割合 $f_{\mathrm{UV}}=0.02$ だけをオゾンが吸収する UV 成分として扱う。残りの $1-f_{\mathrm{UV}}$ は大気を透過する。オゾンは予報せず、気圧に固定した有効光学的厚さを与える。太陽定数は

$$
S_0=1361\,\mathrm{W\,m^{-2}}
$$

とし、軌道離心率は 0 とする。このため太陽と惑星の距離による $S_0$ の変化は考えない。一方、地軸の傾き、惑星の自転、緯度、経度、昼夜および季節による入射角の変化は陽に計算する。

経度を $\lambda$、緯度を $\varphi$、シミュレーション開始からの経過時間を $t$ とする。Octahedral Gaussian Grid の第 $j$ 緯度リングに $M_j$ 点があるとき、各格子点では

$$
\sin\varphi_j=\mu_j,\qquad
\lambda_{i,j}=\frac{2\pi(i-1)}{M_j}
$$

を使う。このケース専用の暦では、1 太陽日を

$$
P_{\mathrm{day}}=86400\,\mathrm{s}
$$

とし、1 年を 360 太陽日、1 ヶ月を 30 太陽日とする。したがって 1 年は 12 ヶ月で、公転周期は

$$
P_{\mathrm{orb}}=360P_{\mathrm{day}}
$$

である。この 360 日暦はこの放射ケースにだけ適用し、他の計算ケースの時間設定は変更しない。円軌道上の公転黄経を

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

この $\Omega$ は放射ケースの惑星渦度 $f=2\Omega\sin\varphi$ と初期値の釣り合いにも用いる。他の計算ケースでは、それぞれ従来の自転角速度を変更せずに用いる。各格子点の太陽時角は

$$
H(\lambda,t)=\theta(t)+\lambda-\alpha_\odot(t)
$$

である。したがって太陽天頂角 $\zeta$ の余弦を

$$
\mu_0(\lambda,\varphi,t)=\cos\zeta
=\sin\varphi\sin\delta_\odot
+\cos\varphi\cos\delta_\odot\cos H
$$

と計算し、大気上端に入る下向き短波フラックスを

$$
F^{\downarrow\mathrm{SW}}_{1/2}(\lambda,\varphi,t)
=S_0\max(0,\mu_0)
$$

とする。$\mu_0\leq0$ の夜側では、すべての短波フラックスと短波加熱を 0 とする。

オゾンの鉛直分布は、対数気圧に関する切断正規分布

$$
w_{\mathrm{O_3}}(p)=
\begin{cases}
\displaystyle
\frac{1}{p}\exp\left[-\frac{\left\{\ln(p/p_{\mathrm{O_3}})\right\}^2}
{2s_{\mathrm{O_3}}^2}\right],
&p_a\leq p\leq p_b,\\[8pt]
0,&\text{otherwise}
\end{cases}
$$

で与える。各定数は

$$
p_a=1\,\mathrm{hPa},\qquad
p_b=100\,\mathrm{hPa},\qquad
p_{\mathrm{O_3}}=10\,\mathrm{hPa},\qquad
s_{\mathrm{O_3}}=\ln 3
$$

とする。これによりオゾンの放射効果は $1$--$100\,\mathrm{hPa}$ に限られ、$10\,\mathrm{hPa}$ 付近で対数気圧当たりの光学的厚さが最大になる。第 $k$ 層に含まれる規格化されたオゾン割合を

$$
q_k=
\frac{\displaystyle\int_{p_{k-1/2}}^{p_{k+1/2}}w_{\mathrm{O_3}}(p)\,dp}
{\displaystyle\int_{p_a}^{p_b}w_{\mathrm{O_3}}(p)\,dp}
$$

とする。積分区間が $[p_a,p_b]$ の外にある部分では被積分関数を 0 とする。この規格化により $\sum_kq_k=1$ となる。UV に対するオゾン層全体の鉛直光学的厚さを

$$
\tau_{\mathrm{O_3,SW}}=1.5
$$

とし、第 $k$ 層の UV 光学的厚さを

$$
\Delta\tau^{\mathrm{O_3,SW}}_k
=\tau_{\mathrm{O_3,SW}}q_k
$$

とする。同じ $q_k$ を、前節の総光学的厚さ $\tau_{\mathrm{O_3,LW}}=0.005$ を持つオゾン長波放射にも用いる。

実装では、積分を誤差関数 `erf` によって解析的に計算した式を用いる。まず

$$
p_L=\max\{p_{k-1/2},p_a\},\quad p_U=\min\{p_{k+1/2},p_b\}
$$

として $p_U\le p_L$ なら $q_k=0$ とする。$p_U>p_L$ なら

$$
q_k=\frac{\mathrm{erf}[x(p_U)]-\mathrm{erf}[x(p_L)]}{\mathrm{erf}[x(p_b)]-\mathrm{erf}[x(p_a)]},\quad x(p)=\frac{\ln(p/p_\mathrm{O_3})}{\sqrt{2}s_\mathrm{O_3}}
$$

とする。

大気上端に入る短波を、UV 成分と非 UV 成分に

$$
F^{\downarrow\mathrm{UV}}_{1/2}
=f_{\mathrm{UV}}F^{\downarrow\mathrm{SW}}_{1/2},\qquad
F^{\downarrow\mathrm{nonUV}}_{1/2}
=(1-f_{\mathrm{UV}})F^{\downarrow\mathrm{SW}}_{1/2}
$$

と分ける。太陽天頂角は大気上端の入射量にだけ反映し、オゾン中の光路は常に鉛直とみなす。したがって各層の UV 透過率は

$$
t^{\mathrm{UV}}_k
=\exp\left(-\Delta\tau^{\mathrm{O_3,SW}}_k\right)
$$

であり、下向き UV フラックスを上端から地表へ

$$
F^{\downarrow\mathrm{UV}}_{k+1/2}
=t^{\mathrm{UV}}_kF^{\downarrow\mathrm{UV}}_{k-1/2}
$$

と順に計算する。非 UV 成分は全層で変化しない。第 $k$ 層が吸収する UV フラックスは

$$
A^{\mathrm{UV}}_k
=F^{\downarrow\mathrm{UV}}_{k-1/2}
-F^{\downarrow\mathrm{UV}}_{k+1/2}
$$

であり、対応する温度変化は

$$
\left(\frac{\partial T_k}{\partial t}\right)_{\mathrm{SW}}
=\frac{g}{c_p\Delta p_k}A^{\mathrm{UV}}_k
$$

とする。大気に与える放射温度変化は、長波と短波の和

$$
\left(\frac{\partial T_k}{\partial t}\right)_{\mathrm{rad}}
=\left(\frac{\partial T_k}{\partial t}\right)_{\mathrm{LW}}
+\left(\frac{\partial T_k}{\partial t}\right)_{\mathrm{SW}}
$$

である。

地表に達する全短波は

$$
F^{\downarrow\mathrm{SW}}_{N+1/2}
=F^{\downarrow\mathrm{nonUV}}_{1/2}
+F^{\downarrow\mathrm{UV}}_{N+1/2}
$$

である。

地表の短波アルベドを $\alpha_s=0.3$ とし、地表が吸収する短波フラックスは

$$
F^{\mathrm{SW}}_{\mathrm{abs}}
=(1-\alpha_s)F^{\downarrow\mathrm{SW}}_{N+1/2}
$$

とする。反射短波は

$$
F^{\mathrm{SW}}_{\mathrm{refl}}
=\alpha_sF^{\downarrow\mathrm{SW}}_{N+1/2}
$$

とし、上向きの反射短波が大気中で再吸収される効果は扱わない。したがって $F^{\mathrm{SW}}_{\mathrm{refl}}$ は地表での反射量と大気上端から出る反射量の両方を表す。短波フラックスは時間平均・日平均・東西平均を事前に取らず、放射傾向を評価する各タイムステップで、その時刻と各格子点の緯度・経度を用いて再計算する。

円軌道かつ空間的に一定な $S_0$ という仮定の下では、任意の時刻における解析的な全球平均入射量は

$$
\left\langle F^{\downarrow\mathrm{SW}}_{1/2}\right\rangle=\frac{S_0}{4}
$$

となる。オゾン吸収はこの大気上端入射量を変えず、地表へ到達する短波だけを減らす。各カラムでは

$$
F^{\downarrow\mathrm{SW}}_{1/2}
=\sum_{k=1}^{N}A^{\mathrm{UV}}_k
+F^{\mathrm{SW}}_{\mathrm{abs}}
+F^{\mathrm{SW}}_{\mathrm{refl}}
$$

が成り立つ。鉛直透過率は太陽天頂角に依存しないため、オゾンによる UV 吸収の解析的な全球平均は

$$
\left\langle\sum_{k=1}^{N}A^{\mathrm{UV}}_k\right\rangle
=f_{\mathrm{UV}}\frac{S_0}{4}
\left(1-e^{-\tau_{\mathrm{O_3,SW}}}\right)
=5.286599\,\mathrm{W\,m^{-2}}
$$

となる。数値実装では、Gaussian 求積による大気上端入射量の全球平均が離散化誤差の範囲で $S_0/4$ に一致すること、および上のカラム短波エネルギー収支をテストで確認する。ただし各格子点に与えるフラックスそのものはこの全球平均値ではなく、上記の瞬時値である。

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

に従って時間発展させる。右辺の放射フラックス、$H_{sd}$、$H_{sa}$ はすべて[物理過程を評価する時刻](#物理過程を評価する時刻)のとおり $\overline X^{n-1}$ の場で評価する。大気層に加える顕熱・放射フラックスと地面層から差し引く値には同じ前時刻のフラックスを使うので、大気・地表・地中を合わせたエネルギー収支は閉じたままである。熱容量は

$$
C_s=2\times10^6\,\mathrm{J\,m^{-2}\,K^{-1}},\quad
C_d=2\times10^7\,\mathrm{J\,m^{-2}\,K^{-1}}
$$

重力波・超粘性の処理については $T_s,T_d$ には適用しない。ただし RAW フィルターについては適用する。

## Rayleigh 摩擦

境界層には [Held–Suarez 強制](./Held-Suarez.md)と同じ Rayleigh 摩擦を加える。さらに、上端の2層（$k=1,2$）にも Rayleigh 摩擦を加え、放射によって励起される上層の発散を減衰させる。

放射ケースの外力を

$$
\bm{F}_k=-k_{v,k}^{(\mathrm{rad})}\overline{\bm{u}}_k^{n-1}
$$

とし、比例定数を

$$
k_{v,k}^{(\mathrm{rad})}
=(1\,\mathrm{day})^{-1}
\left[
\max\left\{0,\frac{\sigma_k-\sigma_b}{1-\sigma_b}\right\}
+I_{\mathrm{top},k}
\right],
$$

$$
I_{\mathrm{top},k}=
\begin{cases}
1, & k=1,2,\\
0, & k=3,\dots,12
\end{cases}
$$

とする。ここで $\sigma_k$ も前時刻の界面気圧と地表面気圧から求め、$\sigma_b=0.7$ は Held–Suarez 強制と同じ定義である。上端 2 層はそれぞれ $1$-$3\,\mathrm{hPa}$ と $3$-$10\,\mathrm{hPa}$ に対応し、風速の減衰時定数は $1\,\mathrm{day}$ である。上層側の追加係数は各層内で一定なので、上端 2 層の相対渦度と発散には前時刻の値からそれぞれ $-\overline\zeta_k^{n-1}/(1\,\mathrm{day})$ と $-\overline\delta_k^{n-1}/(1\,\mathrm{day})$ が加わる。

## 乾燥対流調節

放射によって対流不安定になりがちなので、[乾燥対流調節](./dry-convective-adjustment.md)を行う。

## 超粘性

放射ケースでも、相対渦度、発散、温度には [乾燥大気](./dry.md#超粘性)と同じ超粘性を全層で適用する。上層だけラプラシアンの次数を変えるスポンジ型超粘性は用いず、上層発散の抑制は前節の Rayleigh 摩擦で行う。

## 初期値

[Jablonowski–Williamson の初期状態](Jablonowski-Williamson.md)（以下 JW）をもとに初期化する。温度と地表面気圧 $p_s=p_0=10^5\,\mathrm{Pa}$ は JW と同じとし、$\Omega$ には前節の放射ケースの値を用いる。ただし地形については

$$
\Phi_s=0
$$

とし、平坦地形にする。JW の東西風は $\Phi_s\neq0$ を前提に釣り合っているので、そのままでは $\Phi_s=0$ と整合しない。そこで温度はそのままにして、東西風の方を $\Phi_s=0$ に釣り合うよう調整する。

### 釣り合った東西風

JW の $\eta_v,A(\varphi),B(\varphi)$ を用い、

$$
c(\eta)=\cos^{3/2}\eta_v,\qquad
c_s=\cos^{3/2}\eta_{v,s},\qquad
\eta_{v,s}=\frac{\pi}{2}(1-\eta_0)
$$

と置く。JW の基本場の東西風と、その地表 $\eta=1$ での値を

$$
u_{\mathrm{JW}}(\varphi,\eta)=u_0c(\eta)\sin^2(2\varphi),\qquad
u_{\mathrm{JW},s}(\varphi)=u_0c_s\sin^2(2\varphi)
$$

とする。JW のジオポテンシャルの緯度依存部分を

$$
G(\varphi,\eta)=u_0c(\eta)\left[u_0c(\eta)A(\varphi)+a\Omega B(\varphi)\right]
$$

と書くと、JW の基本場は $\Phi=\overline\Phi(\eta)+G(\varphi,\eta)$ であり、$\Phi_s=G(\varphi,1)$ である。温度を変えずに $\Phi_s=0$ とすると、静水圧の関係から全層のジオポテンシャルが一様に $G(\varphi,1)$ だけずれ、

$$
\Phi(\varphi,\eta)=\overline\Phi(\eta)+G(\varphi,\eta)-G(\varphi,1)
$$

となる。$p_s$ が一様なので $\eta$ 面は等圧面であり、$v=0$ の東西対称な流れが定常であるための条件は、曲率項を含めた傾度風バランス

$$
\left(2\Omega\sin\varphi+\frac{u\tan\varphi}{a}\right)u
=-\frac{1}{a}\frac{\partial\Phi}{\partial\varphi}
$$

である。JW の東西風は $G(\varphi,\eta)$ に対してこの関係を厳密に満たすので、

$$
F(u)=\left(2\Omega\sin\varphi+\frac{u\tan\varphi}{a}\right)u
$$

と置けば、求める東西風 $u_{\mathrm{bal}}$ は

$$
F(u_{\mathrm{bal}})=F(u_{\mathrm{JW}})-F(u_{\mathrm{JW},s})
$$

を満たす。$a/\tan\varphi$ を掛けると $u$ の 2 次方程式

$$
\left(u_{\mathrm{bal}}+a\Omega\cos\varphi\right)^2
=\left(a\Omega\cos\varphi\right)^2+D,\\
\begin{aligned}
D(\varphi,\eta)&=\left(u_{\mathrm{JW}}-u_{\mathrm{JW},s}\right)
\left(u_{\mathrm{JW}}+u_{\mathrm{JW},s}+2a\Omega\cos\varphi\right)\\
&=u_0^2\left[c(\eta)-c_s\right]\sin^2(2\varphi)
\left\{\left[c(\eta)+c_s\right]\sin^2(2\varphi)+\frac{2a\Omega}{u_0}\cos\varphi\right\}
\end{aligned}
$$

になる。$\Phi_s=0$ のときに風速 0 へ連続につながる根を選び、極付近での桁落ちを避けるため有理化して

$$
u_{\mathrm{bal}}(\varphi,\eta)
=\frac{D(\varphi,\eta)}
{a\Omega\cos\varphi+\sqrt{\left(a\Omega\cos\varphi\right)^2+D(\varphi,\eta)}}
$$

とする。$|\eta_v|$ は $0\le\eta\le1$ で $\eta=1$ のとき最大なので $c(\eta)\ge c_s$ であり、$D\ge0$ となって根号の中は常に正である。また $\sin^2(2\varphi)$ の因子により、$u_{\mathrm{bal}}$ は赤道と極で 0 になる。

$u_{\mathrm{bal}}$ は地表 $\eta=1$ で 0 となり、上空ほど JW の東西風に近づく。$N=12$ の full level では、最大値はおおむね緯度 $45^\circ$ にあり、最上層 $\eta=0.002$ で約 $23.2\,\mathrm{m\,s^{-1}}$（JW では約 $31.1$）、$\eta=0.23$ で約 $27.2\,\mathrm{m\,s^{-1}}$（約 $35.0$）、最下層 $\eta=0.94$ で約 $3.0\,\mathrm{m\,s^{-1}}$（約 $11.3$）となる。

東西風には、局所的な擾乱を加えず

$$
u=u_\mathrm{bal},v=0
$$

を初期値とし、$\zeta,\delta$ をスペクトル法で求める。温度は JW の式をそのまま用いる。また、$T_s,T_d$ は $T_N$ と等しくする。

## 時間・出力

$5$ 年間、すなわち $5\times360$ 太陽日シミュレーションする。タイムステップは $\Delta t=1200\,\mathrm{s}$ とする。このため 1 日は 72 ステップ、1 ヶ月は 2160 ステップ、1 年は 25920 ステップにそれぞれ厳密に一致する。

1日ごとに次の日平均量を出力する。

- 大気温度の全球平均 $\braket{T}$
- 地表温度の全球平均 $\braket{T_s}$
- 地面下層温度の全球平均 $\braket{T_d}$
- 大気の全球平均運動エネルギー
- 地表面気圧の全球平均 $\braket{p_s}$
- 大気上端における全球の短波入射
- 大気上端から出る全球の反射短波
- 全球の出ていった長波

これらの量を求めるために別でスペクトル法を用いるとパフォーマンスが悪化するので計算結果を流用してオンラインで集計する。

日平均は各太陽日の左端を含み右端を含まない全タイムステップ（72 サンプル）を等重みで平均し、時刻の列にはその日の開始時刻を書く。

30 日ごとに次の月平均量を出力する。最初の月は 1 年目の 4 月 1 日から 4 月 30 日までである。

- 各格子点の地表温度の月平均 $T_s(\lambda,\varphi)$
- 各格子点の地表面気圧の月平均 $p_s(\lambda,\varphi)$
- 経度方向に平均した温度・流速の月平均 $[T](\varphi,\eta),[u](\varphi,\eta),[v](\varphi,\eta)$
- 上の平均からの偏差をプライムをつけて表したときの $[u'v'](\varphi,\eta),[v'T'](\varphi,\eta)$

これらもオンラインで効率的に計算する。

月平均は各 30 日区間の左端を含み右端を含まない全タイムステップ（2160 サンプル）を等重みで平均する。`m0001` は 1 年目の 4 月、`m0002` は 5 月に対応し、以後も30日ごとに連番とする。

毎年の 4/1 に次の瞬時値を出力する。

- スペクトル $\zeta,\delta,T,\ln p_s$
- 格子 $\zeta,\delta,T,u,v,\ln p_s,T_s,T_d$

初期値を含め、最終的に6つの瞬時値が出るようにする。

瞬時値のファイル名には `yearly_` を付け、開始時刻を `y0001`、以後の各 4 月 1 日を `y0002` から `y0006` として保存する。
