# Octahedral Gaussian Grid

## Gauss-Legendre 求積

球面の緯度方向を $L=2G$ 個に分ける。また、$T=G-1$ とする。Legendre 多項式 $P_L$ の零点を小さい方から順に $\mu_j,\ j=1,\dots,L$ と置く。対応する緯度は

$$
\sin\varphi_j=\mu_j
$$

である。Gauss 重みを

$$
w_j=\frac{2}{(1-\mu_j^2)P_L'(\mu_j)^2}
$$

とすると

$$
\int_{-1}^1f(\mu)d\mu\approx\sum_{j=1}^Lw_jf(\mu_j)
$$

と近似でき、$f$ が $2L-1$ 次多項式までなら厳密に等式が成り立つ。

零点と重みは行列

$$
\beta_k=\frac{k}{\sqrt{4k^2-1}},\quad k=1,2,\dots,\\
J_L=\begin{pmatrix}
  0 & \beta_1 & & & \\
  \beta_1 & 0 & \beta_2 & & \\
  & \beta_2 & 0 & \ddots & \\
  & & \ddots & \ddots & \beta_{L-1} \\
  & & & \beta_{L-1} & 0
\end{pmatrix}
$$

の固有値・固有ベクトルから求められる。零点 $\mu_j$ は固有値そのものであり、対応する重みは固有ベクトル $v_j$ の第1成分から

$$
w_j=2(v_j^1)^2
$$

となる。

## 経度の分割

経度方向は各緯度 $\mu_j$ で

$$
M_j=4\min\{j,L+1-j\}+16
$$

個に分割し、

$$
\lambda_{j,k}=\frac{2\pi k}{M_j},\quad k=0,1,\dots,M_j-1
$$

とする。陪 Legendre 関数をここでは

$$
P^m_n(\mu)=\sqrt{(2n+1)\frac{(n-m)!}{(n+m)!}}(1-\mu^2)^{m/2}\frac{d^mP_n}{d\mu^m}
$$

と定義すると直交関係

$$
\frac{1}{2}\int_{-1}^1P^m_n(\mu)P^m_{n'}(\mu)d\mu=\delta_{nn'}
$$

が成り立つ。対称性は

$$
P^m_n(-\mu)=(-1)^{n-m}P^m_n(\mu)
$$

として成り立つ。

## 連続場の変換

スカラー場 $f(\lambda,\mu)$ を

$$
f(\lambda,\mu)=\sum_{n=0}^\infty\sum_{m=-n}^ne^{im\lambda}a_n^m P_n^{|m|}(\mu)
$$

と展開する。逆変換は

$$
a_n^m=\frac{1}{4\pi}\int_{-1}^1d\mu P^{|m|}_n(\mu)\int_0^{2\pi}d\lambda f(\lambda,\mu)e^{-im\lambda}
$$

となる。

## 離散場の変換

スカラー場を離散化して

$$
f_{j,k}=f(\lambda_{j,k},\mu_j)
$$

とする。このとき

$$
a_n^m=\frac{1}{2}\sum_{j=1}^Lw_j\chi_j^mP^{|m|}_n(\mu_j)
\frac{1}{M_j}\sum_{k=0}^{M_j-1}f_{j,k}e^{-im\lambda_{j,k}},\\
\chi_j^m=\begin{cases}
  1, & |m| \le m_j^\text{max}\\
  0, & |m| > m_j^\text{max}
\end{cases},\quad m_j^\text{max}=\frac{M_j}{2}-1
$$

として変換を定義する。逆変換は

$$
f_{j,k}=\sum_{m=-T}^T\chi_j^m e^{im\lambda_{j,k}}\sum_{n=|m|}^Ta_n^mP_n^{|m|}(\mu_j)
$$

で与えられるが、正確には互いに逆にはならない。

## スペクトル微分

球面調和展開されたスカラー場

$$
f(\lambda,\varphi)
=\sum_{n=0}^T\sum_{m=-n}^n
e^{im\lambda}a_n^mP_n^{|m|}(\mu),
\qquad \mu=\sin\varphi
$$

の微分は展開係数と基底関数から求まる。

経度微分では Fourier 基底を解析的に微分できるため、各係数に $im$ を掛ければよい：

$$
\frac{\partial f}{\partial\lambda}
=\sum_{n=0}^T\sum_{m=-n}^n
e^{im\lambda}(im a_n^m)P_n^{|m|}(\mu)
$$

したがって経度微分用の係数は $d_n^m=im a_n^m$ である。

緯度微分では $\partial/\partial\varphi=\cos\varphi\,\partial/\partial\mu$ を使う。ここで

$$
\epsilon_{n,p}
=\sqrt{\frac{n^2-p^2}{4n^2-1}},
\qquad n\geq1,\qquad 0\leq p\leq n,
$$

と置き、$\epsilon_{0,0}=0$ と定義する。正規化済み陪 Legendre 関数の漸化式から

$$
(1-\mu^2)\frac{dP_n^p}{d\mu}
=(n+1)\epsilon_{n,p}P_{n-1}^p
-n\epsilon_{n+1,p}P_{n+1}^p
$$

が成り立つ。$p=|m|$ として次数ごとに項を集め直すと、

$$
\cos\varphi\frac{\partial f}{\partial\varphi}
=\sum_{m=-T}^T e^{im\lambda}
\sum_{n=|m|}^{T+1}b_n^mP_n^{|m|}(\mu)
$$

と書ける。その係数は

$$
b_n^m
=(n+2)\epsilon_{n+1,|m|}a_{n+1}^m
-(n-1)\epsilon_{n,|m|}a_{n-1}^m,
\qquad |m|\leq n\leq T+1,
$$

である。ただし、範囲 $|m|\leq k\leq T$ の外にある入力係数は $a_k^m=0$ とする。

実装では、各 $m$ について $b_{|m|:T+1}^m$ をゼロで初期化し、各入力係数 $a_n^m$（$n=|m|,\dots,T$）に対して

$$
\begin{aligned}
b_{n-1}^m &\leftarrow b_{n-1}^m+(n+1)\epsilon_{n,|m|}a_n^m
&& (n>|m|),\\
b_{n+1}^m &\leftarrow b_{n+1}^m-n\epsilon_{n+1,|m|}a_n^m
&& (n\geq|m|)
\end{aligned}
$$

と更新すればよい。

$b_{T+1}^m=-T\epsilon_{T+1,|m|}a_T^m$ は一般にはゼロでないため、緯度微分では $n=T+1$ を使わなければならない。
