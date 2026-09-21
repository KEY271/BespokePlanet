# 地形

[乾燥大気](./dry.md)・[湿潤大気](./moist.md)の下端境界条件である地表ジオポテンシャル $\Phi_s=gz_s$ と、陸の分布を表す陸面率 $f_L$ を、解析的な形状の重ね合わせとして格子上に生成する方法を定める。生成は切断波数 $T$ に依存しない角度（度）で定義し、$T=31$ と $T=63$ の両方で同じ惑星を作れるようにする。$\Phi_s$ は[乾燥大気](./dry.md#初期化)のとおり積分中ずっと固定し、超粘性や RAW フィルターは適用しない。$f_L$ は[陸と海の混合](../tendency/land-sea-surface.md)の地表収支だけが使う格子量で、スペクトル変換はしない。

## 形状の定義

経度 $\lambda$、緯度 $\varphi$、形状の中心 $(\lambda_c,\varphi_c)$、向き $\theta$（東向きから反時計回り）を度で与える。中心を原点とする局所座標を

$$
\Delta\lambda=\bigl((\lambda-\lambda_c+180^\circ)\bmod360^\circ\bigr)-180^\circ,\qquad
x=\cos\varphi\,\Delta\lambda,\qquad
y=\varphi-\varphi_c,
$$

$$
x'=x\cos\theta+y\sin\theta,\qquad
y'=-x\sin\theta+y\cos\theta
$$

とする。$x$ に評価点の $\cos\varphi$ を掛けるので、形状の東西幅は緯度によらず同じ実長さになる。局所座標は正距円筒図法上の近似であり、極を含む形状には後述の極冠を使う。

### 楕円大陸

半軸 $a,b$（度）の楕円大陸は、規格化した距離

$$
d=\sqrt{\left(\frac{x'}{a}\right)^2+\left(\frac{y'}{b}\right)^2}
$$

と、海岸の幅 $w$（度）を使って

$$
s(\lambda,\varphi)=\frac12\left[1-\tanh\frac{(d-1)\min(a,b)}{w}\right]
$$

とする。内部で $s\to1$、外部で $s\to0$ となり、$d=1$ の海岸線を横切る $10\%$–$90\%$ の遷移幅は短軸方向で約 $2.2w$ である。

### 極冠

北極または南極を中心とする陸は、縁の緯度 $\varphi_0>0$ と符号 $\varsigma=+1$（北）、$\varsigma=-1$（南）で

$$
s(\lambda,\varphi)=\frac12\left[1+\tanh\frac{\varsigma\varphi-\varphi_0}{w}\right]
$$

とする。

### 山脈

山脈は中心 $(\lambda_c,\varphi_c)$、稜線の向き $\theta$、長さ $L$、半幅 $\sigma$（いずれも度）、高さ $H$（m）を持ち、稜線に沿う座標 $x'$ と直交する座標 $y'$ で

$$
z_r(\lambda,\varphi)=H\exp\left[-\left(\frac{y'}{\sigma}\right)^2\right]\frac12\left[1-\tanh\frac{|x'|-L/2}{\sigma}\right]
$$

とする。

## 陸面率と地表高度

大陸と極冠の形状 $s_i$ の和集合を陸面率

$$
f_L=1-\prod_i\left(1-s_i\right),\qquad0\le f_L\le1
$$

とする。重なりがあっても $1$ を超えない。地表高度は、陸の基準高度 $h_0$ と山脈の和に陸面率を掛けて

$$
z_s=f_L\left(h_0+\sum_rz_r\right),\qquad
\Phi_s=gz_s
$$

とする。海上（$f_L=0$）では $z_s=0$ であり、海岸付近の山脈は $f_L$ によって自然に裾を引く。

$z_s$ は格子で計算してから $\Phi_s$ をスペクトル変換し、他の変数と同じ切断波数で切った値を用いる（[Jablonowski–Williamson](./Jablonowski-Williamson.md#地表ジオポテンシャル)と同じ扱い）。以下で「格子の $\Phi_s$」と言うときは、この切断後のスペクトル $\Phi_s$ を格子に戻した値を指す。切断により海上でも $\Phi_s$ は厳密には $0$ でなく、数 m 程度のリップルを持つ。$f_L$ は切断しないので、格子の $f_L$ は解析形そのままである。

## 解像度との関係

切断波数 $T$ で表せる最短波長は $360^\circ/T$ で、$T=31$ では約 $11.6^\circ$、$T=63$ では約 $5.7^\circ$ である。これより細かい構造は切断で落ち、その分だけ Gibbs リップルが出る。形状のパラメータは次を満たすようにとる。

- 海岸の幅と山脈の半幅は最短波長の半分以上、$w,\sigma\ge180^\circ/T$（$T=31$ で約 $5.8^\circ$、$T=63$ で約 $2.9^\circ$）。
- 楕円の短軸と山脈の長さは最短波長以上、$\min(a,b),L\ge360^\circ/T$。

同じ惑星を $T=31$ と $T=63$ で計算する場合は $T=31$ の条件で定義しておけば $T=63$ でも使える。$T=63$ だけで使う地形には $w,\sigma$ を $3^\circ$ まで細くできる。生成時には次の量を求めてメタデータに記録し、切断の影響を確認する。

- 全球の陸面率 $\braket{f_L}$。Gauss 重み $w_j$ と緯度リング $j$ の経度数 $N_j$ で $\braket{f_L}=\frac12\sum_jw_j\frac{1}{N_j}\sum_if_L(\lambda_i,\varphi_j)$ とする。
- 格子の $z_s=\Phi_s/g$ の最大値と最小値。最小値は切断による海上の負の行き過ぎを表す。
- 解析形の $z_s$ と格子の $z_s$ の差の全球 RMS。

## 既定の惑星

既定では次の形状を使う。角度は度、高さは m である。全球の陸面率は約 $27\%$ になる。

| 種類 | 中心 $(\lambda_c,\varphi_c)$ | 大きさ | 向き $\theta$ | その他 |
| --- | --- | --- | --- | --- |
| 楕円大陸 A | $(60^\circ,45^\circ)$ | $a=55,\ b=28$ | $0$ | |
| 楕円大陸 B | $(300^\circ,-10^\circ)$ | $a=28,\ b=40$ | $0$ | |
| 楕円大陸 C | $(150^\circ,-30^\circ)$ | $a=22,\ b=16$ | $0$ | |
| 南極冠 | — | $\varphi_0=70$ | — | $\varsigma=-1$ |
| 山脈 1（大陸 A） | $(75^\circ,35^\circ)$ | $L=50,\ \sigma=6$ | $0$（東西） | $H=2500$ |
| 山脈 2（大陸 B） | $(295^\circ,-10^\circ)$ | $L=60,\ \sigma=6$ | $90$（南北） | $H=2500$ |

海岸の幅は $w=6^\circ$、陸の基準高度は $h_0=300\,\mathrm{m}$ とする。$w=\sigma=6^\circ$ は $T=31$ の条件 $180^\circ/31\approx5.8^\circ$ をぎりぎり満たす値で、$T=63$ ではより余裕がある。形状の一覧とパラメータは他の設定と同じく Fortran のデフォルト初期化子に一か所で書き、ケースが上書きする。

## 地形上の初期状態

[放射ケース](../cases/radiation.md)は $\Phi_s=0$ と一様な $p_s=p_0$ のもとで、JW の温度に釣り合う東西風 $u_{\mathrm{bal}}$ を使っている。地形を入れると $p_s$ は一様でなくなり、$\eta$ 面は等圧面でなくなるので、同じ手順では釣り合わない。そこで放射ケースの基本場を気圧の関数として読み替え、地表面気圧を静水圧の関係から決める。

### 気圧の関数としての基本場

$\eta_p=p/p_0$ と置き、[JW](./Jablonowski-Williamson.md) の水平平均温度 $\overline T$ と $A(\varphi),B(\varphi)$、放射ケースの $G(\varphi,\eta)$ を使って

$$
\Phi^\ast(\varphi,p)=\overline\Phi(\eta_p)+G(\varphi,\eta_p)-G(\varphi,1),\qquad
T^\ast(\varphi,p)=T_{\mathrm{JW}}(\varphi,\eta_p),\qquad
u^\ast(\varphi,p)=u_{\mathrm{bal}}(\varphi,\eta_p)
$$

とする。ここで $\overline\Phi$ は $\overline\Phi(1)=0$ を満たす水平平均のジオポテンシャルで、対流圏 $\eta_t\le\eta$ では

$$
\overline\Phi(\eta)=\int_\eta^1R_d\overline T(\eta')\,d\ln\eta'
=\frac{gT_0}{\Gamma}\left(1-\eta^{R_d\Gamma/g}\right)
$$

である。成層圏の項は $\eta<\eta_t=0.2$ にしか現れず、地表面気圧の決定には関わらない。$\Phi^\ast,T^\ast,u^\ast$ は $p_s=p_0$ のときの放射ケースの場そのもので、気圧座標で見た傾度風バランスを満たす。静水圧の関係 $\partial\Phi^\ast/\partial\ln p=-R_dT^\ast$ も満たすので、これを地形の上に乗せた状態は、離散化誤差を除いて釣り合っている。

### 地表面気圧

各格子点で、格子の $\Phi_s$ に対して

$$
\Phi^\ast(\varphi,p_s)=\Phi_s(\lambda,\varphi)
$$

を満たす $p_s$ を求める。$\Phi^\ast$ は $p$ について単調減少なので根は一意である。反復

$$
p_s^{(m+1)}=p_s^{(m)}\exp\left[\frac{\Phi^\ast(\varphi,p_s^{(m)})-\Phi_s}{R_dT^\ast(\varphi,p_s^{(m)})}\right],\qquad
p_s^{(0)}=p_0\exp\left(-\frac{\Phi_s}{R_dT_0}\right)
$$

は等温大気では 1 回で収束し、JW の温度でも数回で $|\Phi^\ast(p_s)-\Phi_s|<10^{-6}\,\mathrm{m^2\,s^{-2}}$ に達する。$\Phi_s=0$ の格子点では $\overline\Phi(1)=0$、$G(\varphi,1)-G(\varphi,1)=0$ により補正の指数が厳密に $0$ となるので $p_s=p_0$ が厳密に得られ、切断のリップルで $\Phi_s<0$ となった点では $p_s>p_0$ になる。$z_s=2000\,\mathrm{m}$ では $p_s\approx786\,\mathrm{hPa}$ である。

### 各層の値

求めた $p_s$ から [hybrid-sigma 座標](./hybrid-sigma.md)の half level の気圧 $p_{k+1/2}=A_{k+1/2}+B_{k+1/2}p_s$ を作り、full level の気圧を

$$
\eta_{p,k}=\frac{p_{k-1/2}+p_{k+1/2}}{2p_0}
$$

とする。$p_s=p_0$ ではこれは JW の $\eta_k$ に一致する。各層の温度と風を

$$
T_k=T^\ast(\varphi,p_0\eta_{p,k}),\qquad
u_k=u^\ast(\varphi,p_0\eta_{p,k}),\qquad v_k=0
$$

とし、$\zeta,\delta$ をスペクトル法で求める。$\ln p_s$ は格子で計算してからスペクトル変換する。地表面温度と地中温度の初期値は放射ケースと同じく最下層の温度 $T_N$ にとる。地形の上では $T_N$ が JW の減率 $\Gamma=5\,\mathrm{K\,km^{-1}}$ の分だけ低くなり、$2000\,\mathrm{m}$ で約 $10\,\mathrm{K}$ である。[湿潤ケース](../cases/moist.md#初期値)の比湿の初期値は各格子点の $p_s$ から作る full level の気圧を使っているので変更はいらない。

$\Phi_s\equiv0$ のときこの手順は放射ケースの初期化にそのまま戻る。したがって既存の放射・slab ocean・湿潤ケースの初期値をこの手順で作り直しても、結果はビット単位で変わらない。この一致は実装後の確認項目にする。

### 残る不釣り合い

- $\ln p_s$ と $\Phi_s$ はそれぞれ切断されており、力学の気圧傾度項 $\nabla\Phi+R_dT\nabla\ln p_s$ の打ち消しは切断後の場では厳密ではない。
- 力学の離散化した静水圧の関係（$\alpha_k,L_k$）と、連続な $\overline\Phi$ で決めた $p_s$ の間には $O(\Delta\eta^2)$ の差がある。これは平坦な放射ケースにもある離散化誤差と同じ種類のものである。
- hybrid の $A_{k+1/2}\ne0$ の層では $\eta$ 面が地形に完全には追従しないが、力学の気圧傾度項は座標に依らない形なので、初期の釣り合いは連続の極限で保たれる。

これらの不釣り合いは重力波として数日で調整され、[湿潤ケース](../cases/moist.md#初期値)が仮想温度の分の不釣り合いを許容しているのと同じく、数年の積分では無視できる。地形を入れた乾燥・物理過程なしの短い積分で、平坦な場合と比べて発散や南北風が桁で大きくならないことを確認する。

## 力学への影響

- 気圧傾度項の $-\nabla^2_\eta\Phi_s$ は[乾燥大気](./dry.md#重力波)のとおり右辺 $\mathcal R(X)$ にだけ現れ、重力波の陰的演算子は変わらない。
- 重力波演算子の参照状態は $p_s=p_0$、$\widetilde T=300\,\mathrm{K}$ のままとする。高原の上では実際の $p_s$ が参照より $20\%$ 程度低く、その差は陽的に扱われる。参照温度が実際より高い側にあるので陰解法の安定性は保たれる。
- 超粘性は $\eta$ 面に沿って作用するので、斜面の上では温度の超粘性が鉛直方向の混合を伴う。地形に沿った補正は行わない。
- [Held–Suarez 強制](../tendency/Held-Suarez.md#rayleigh-摩擦)の地表摩擦は $\sigma=p/p_s$ で定義されており、地形に追従する。[長波放射](../tendency/longwave-radiation.md)の層の光学的厚さは $\Delta p_k$ に比例するので、地形の上では気柱が薄い分だけ自動的に小さくなる。

## 出力

生成した地形はケースの出力に、時間に依らない格子の場として一度だけ書く。

- 陸面率 $f_L(\lambda,\varphi)$
- 格子の地表高度 $z_s(\lambda,\varphi)=\Phi_s/g$（切断後、m）

メタデータには形状の一覧とパラメータ、全球の陸面率、格子の $z_s$ の最大・最小、解析形との RMS 差を書く。
