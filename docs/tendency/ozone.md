# オゾン

オゾンは予報せず、気圧に固定した有効光学的厚さを与える。[帯域別放射](./band-radiation.md)は、下の鉛直分布の規格化した層内割合 $q_k$ をそのまま使い、気柱量 $300\,\mathrm{DU}$ と帯域ごとの質量吸収係数から光学的厚さを作る（同ページの 3 節）。以下の光学的厚さ $\tau_{\mathrm{O_3,SW}}$、$\tau_{\mathrm{O_3,LW}}$ は灰色スキームだけが使う。鉛直層の添字と界面気圧は[長波放射](./longwave-radiation.md)と同じとする。

## 鉛直分布

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

とする。積分区間が $[p_a,p_b]$ の外にある部分では被積分関数を 0 とする。この規格化により $\sum_kq_k=1$ となる。

## 光学的厚さ

UV に対するオゾン層全体の鉛直光学的厚さを

$$
\tau_{\mathrm{O_3,SW}}=1.5
$$

とし、第 $k$ 層の UV 光学的厚さを

$$
\Delta\tau^{\mathrm{O_3,SW}}_k
=\tau_{\mathrm{O_3,SW}}q_k
$$

とする。同じ $q_k$ を用いて、[長波放射](./longwave-radiation.md)のオゾン長波光学的厚さを

$$
\Delta\tau^{\mathrm{O_3,LW}}_k
=\tau_{\mathrm{O_3,LW}}q_k,\qquad
\tau_{\mathrm{O_3,LW}}=0.005
$$

とする。

## 実装

実装では、積分を誤差関数 `erf` によって解析的に計算した式を用いる。まず

$$
p_L=\max\{p_{k-1/2},p_a\},\quad p_U=\min\{p_{k+1/2},p_b\}
$$

として $p_U\le p_L$ なら $q_k=0$ とする。$p_U>p_L$ なら

$$
q_k=\frac{\mathrm{erf}[x(p_U)]-\mathrm{erf}[x(p_L)]}{\mathrm{erf}[x(p_b)]-\mathrm{erf}[x(p_a)]},\quad x(p)=\frac{\ln(p/p_\mathrm{O_3})}{\sqrt{2}s_\mathrm{O_3}}
$$

とする。

気柱のすべての層をまとめて求めるとき（`ozone_layer_fractions`）は、隣り合う層が共有する界面の $\mathrm{erf}$ を 1 回だけ計算する。$p_a$ より上と $p_b$ より下の界面は $p_a$、$p_b$ に切り詰められるので、分母で求める $\mathrm{erf}[x(p_a)]$、$\mathrm{erf}[x(p_b)]$ をそのまま使う。$\mathrm{erf}$ を呼ぶのはオゾン層の中にある界面だけになり、結果は層ごとの式とビット単位で一致する。
