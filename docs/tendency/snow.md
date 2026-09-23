# 降雪と積雪

状態: **実装済み**。[大規模凝結](./large-scale-condensation.md)による降雪、暖かい層での融解、陸の積雪バケツ、積雪のアルベド・地表交換への効果、地表での融雪、海に落ちた雪の扱い、時間積分と出力を記述する。数値検証は `core/test/check_snow.f90` にある。

温度はすべて K で書く。$-10^\circ\mathrm C$ は $263.15\,\mathrm K$、$5^\circ\mathrm C$ は $278.15\,\mathrm K$ である。

## 1. 扱う範囲と設定

| 項目 | 仕様 |
| --- | --- |
| 雪を作る過程 | 大規模凝結だけ。[湿潤対流調節](./moist-convective-adjustment.md)の降水は常に雨 |
| 相の判定 | 凝結が生じた層の温度が $T_{\rm sn}=263.15\,\mathrm K$ より低ければ、その層の凝結はすべて雪 |
| 雪の潜熱 | 雪になる凝結は $L_s=L_v+L_f$ を放出する（2 節） |
| 大気中の融解 | 雪を上から順に落とし、$T_{\rm am}=278.15\,\mathrm K$ より暖かい層では、$T_{\rm am}$ を超える分の顕熱で雪を融かして雨にする |
| 地表の積雪 | 陸だけに積雪バケツ $S$（陸面積当たりの水当量）を置く。海氷上にはためない |
| 積雪率 | $f=S/(S+S_0)$、$S_0=0.05\,\mathrm m$ 水当量 $=50\,\mathrm{kg\,m^{-2}}$ |
| 陸のアルベド | $\alpha_L+\Delta\alpha\,f$、$\Delta\alpha=0.4$ |
| 地表交換の抑制 | 陸の顕熱と水蒸気交換（蒸発・結露）を $1-f$ 倍にする |
| 地表の融雪 | 更新後の陸温度が $T_{\rm gm}=275\,\mathrm K$ より暖かければ、超えた分の熱で積雪を融かす |
| 融雪水 | 陸の[バケツ](./bucket.md)へ入れる |
| 海に落ちた雪 | 直ちに海洋混合層の熱で融かす（6 節）。開水面か海氷かを区別しない |
| 有効化 | 陸海ケース（`land`, `land-t63`, `land-earth`, `land-earth-t63`）。湿潤 aquaplanet `moist` では無効のまま |

$\Delta\alpha$ の既定値により、陸のアルベドは雪なしの $0.2$ から、$S=S_0$ で $0.4$、$S\to\infty$ で $0.6$ に近づく。$f<1$ なので陸が完全に雪で覆われることはない。

定数は $L_v=2.5\times10^6\,\mathrm{J\,kg^{-1}}$（[飽和比湿](./saturation-specific-humidity.md)）、$L_f=3.34\times10^5\,\mathrm{J\,kg^{-1}}$、$L_s=L_v+L_f=2.834\times10^6\,\mathrm{J\,kg^{-1}}$、$c_p=1004\,\mathrm{J\,kg^{-1}\,K^{-1}}$、$g=9.80616\,\mathrm{m\,s^{-2}}$、水の密度 $\rho_w=1000\,\mathrm{kg\,m^{-3}}$ とする。$S$ の単位 $\mathrm{kg\,m^{-2}}$ は水当量の $\mathrm{mm}$ に等しく、$S_0=0.05\,\mathrm m\times\rho_w=50\,\mathrm{kg\,m^{-2}}$ である。雪の密度は持たないので、$S$ は積雪深ではなく水当量である。

次は扱わない。対流性の雪、雨の再凍結、落下中の雪の昇華と雨の再蒸発、積雪面の昇華・凝華、積雪の熱容量・密度・圧密・上限と氷床化、積雪アルベドの経時変化、雪の水平移動（吹雪・雪崩）、海氷上の積雪。飽和比湿は従来どおり水面に対する値を使い、氷面の飽和は区別しない。雪で覆われた部分と覆われていない部分に別々の温度は持たせず、陸タイルは一つの温度 $T_L$ のままとする。

## 2. 相の判定と凝結

大規模凝結は[湿潤大気](../dynamics/moist.md#物理過程)の暫定場 $T^{(2)}_k,q^{(2)}_k$（乾燥・湿潤対流調節の傾向を加えた場）で評価する。凝結が生じた層の温度とは、この凝結前の暫定温度 $T^{(2)}_k$ とする。第 $k$ 層の潜熱を

$$
L_k=
\begin{cases}
L_s=L_v+L_f,& T^{(2)}_k<T_{\rm sn},\\
L_v,& T^{(2)}_k\ge T_{\rm sn}
\end{cases}
$$

とし、凝結量 $C_k$ を[大規模凝結](./large-scale-condensation.md#凝結量)の Newton 法で $L$ を $L_k$ に置き換えて

$$
q^{(2)}_k-C_k=q_s\!\left(T^{(2)}_k+\frac{L_k}{c_p}C_k,\;p_k\right)
$$

から求める。$T^{(2)}_k<T_{\rm sn}$ の層の $C_k$ はすべて雪、それ以外の層の $C_k$ はすべて雨である。雪の層も凝結後にちょうど飽和になる。

雪の層に $L_v$ ではなく $L_s$ を与えるのは、後で雪を融かすときに $L_f$ を使うからである。$L_v$ だけで雪を作ると、雪を作った層は雨の場合と同じだけ暖まるのに、融かす層や地表は $L_f$ を失い、降雪・融雪を経るたびに系から $L_f$ が消える。$L_s$ で作れば、雪が融けて雨になった後の大気・地表の合計は、はじめから雨として凝結した場合と同じになる（4 節）。

## 3. 落下と大気中の融解

第 1 層を最上層とし、上端から第 $k$ 層に入る雪のフラックスを $F_{k-1/2}$（$\mathrm{kg\,m^{-2}\,s^{-1}}$）、$F_{1/2}=0$ とする。場を進める時間幅を $\Delta\tau$（通常 $2\Delta t$）、層の質量流量換算を $m_k=\Delta p_k/(g\Delta\tau)$ と書く。第 $k$ 層の雪の生成は

$$
G_k=
\begin{cases}
m_kC_k,& T^{(2)}_k<T_{\rm sn},\\
0,&\text{それ以外}
\end{cases}
$$

である。層自身の凝結で暖まった後の温度

$$
\widehat T_k=T^{(2)}_k+\frac{L_k}{c_p}C_k
$$

が $T_{\rm am}$ を超える分の顕熱で融かせる雪の量を

$$
H_k=\frac{c_p\max(\widehat T_k-T_{\rm am},0)}{L_f}\,m_k
$$

とし、融解量と下端のフラックスを

$$
M_k=\min(F_{k-1/2}+G_k,\;H_k),\qquad
F_{k+1/2}=F_{k-1/2}+G_k-M_k
$$

とする。$M_k$ は融けて雨になり、そのまま地表まで落ちる。融かせる熱が雪より多ければ雪はすべて融け、少なければ層は $T_{\rm am}$ ちょうどまで冷える。したがって融解で層が $T_{\rm am}$ より冷えることはない。$T_{\rm sn}\le\widehat T_k\le T_{\rm am}$ の層では雪は変化せずに通過し、雨も凍らない。

層の傾向は

$$
\left(\frac{\partial q_k}{\partial t}\right)_{\rm ls}=-\frac{C_k}{\Delta\tau},\qquad
\left(\frac{\partial T_k}{\partial t}\right)_{\rm ls}=\frac{L_kC_k}{c_p\Delta\tau}-\frac{L_fM_k}{c_pm_k\Delta\tau}
$$

である。地表に届く大規模凝結の降水、そのうちの雪、気柱で融けた雪は

$$
P_{\rm ls}=\sum_km_kC_k,\qquad
S_{\rm ls}=F_{N+1/2},\qquad
M_{\rm atm}=\sum_kM_k
$$

で、$P_{\rm ls}-S_{\rm ls}$ が雨である。$P_{\rm ls}$ の定義は雪を入れる前と同じであり、雲の診断や陸の降水の集計は従来どおり $P_{\rm conv}+P_{\rm ls}$ を使う。

雪の層で $\widehat T_k$ が $T_{\rm am}$ を超えるには、$T_{\rm sn}$ 未満から $15\,\mathrm K$ 以上の凝結加熱が必要で、実際には起こらない。それでも式は同じ層で作った雪も融解の対象に含め、エネルギーを保存する。

目安として、厚さ $\Delta p=100\,\mathrm{hPa}$ の層が $T_{\rm am}$ を $1\,\mathrm K$ 超えていると、$c_p\Delta p/(gL_f)\simeq3.07\,\mathrm{kg\,m^{-2}}$、すなわち水当量 $3\,\mathrm{mm}$ 程度の雪を 1 回の更新で融かせる。

## 4. 気柱のエネルギー

各層の $c_pT+L_vq$ の変化は、雨の層で 0、雪の層で $L_fG_k\Delta\tau$、融解の層で $-L_fM_k\Delta\tau$ である。気柱で和をとると

$$
\sum_k\frac{\Delta p_k}{g}\left[c_p\left(\frac{\partial T_k}{\partial t}\right)_{\rm ls}+L_v\left(\frac{\partial q_k}{\partial t}\right)_{\rm ls}\right]
=L_f\left(\sum_kG_k-\sum_kM_k\right)=L_fS_{\rm ls}
$$

となる。すなわち大気は、地表に届いた雪の融解熱 $L_fS_{\rm ls}$ だけ $c_pT+L_vq$ を得る。この熱は雪が地表で融けるときに地表が支払う（5, 6 節）。雪が地表に届かなければ、大規模凝結は従来どおり層ごとではなく気柱全体として $c_pT+L_vq$ を保存する。

## 5. 陸の積雪

### 5.1 予報変数と積雪率

陸面積当たりの積雪水当量を

$$
S(\lambda,\varphi,t)\ge0\quad[\mathrm{kg\,m^{-2}}]
$$

とする。バケツの $W$ と同じく陸部分 $1\,\mathrm{m^2}$ 当たりの量で、格子全体の積雪量は $f_LS$ である。純海洋格子では $S=0$ に固定する。積雪率は

$$
f=\frac{S}{S+S_0}
$$

で、$S=0$ なら $f=0$、$S=S_0$ なら $f=1/2$ である。

### 5.2 地表フラックスへの効果

[陸・開水面・海氷の表面フラックス](./sea-ice.md#4-陸開水面海氷の表面フラックス)の陸の式を次のように変える。

$$
F_L=(1-\alpha_L-\Delta\alpha f)S_\downarrow+L_\downarrow-\sigma T_L^4-(1-f)H_L(T_L)-L_vE_L,
$$

$$
E_L=\min\!\left(\beta_L(W)(1-f)E_\ast(T_L),\;\frac{W}{\Delta\tau}\right)\quad(\text{上向きのとき}),\qquad
E_L=\beta_L(W)(1-f)E_\ast(T_L)\quad(\text{下向きのとき}).
$$

$H_L$ は[地面](./ground.md)のバルク式の顕熱、$E_\ast$ は[バケツ](./bucket.md#蒸発)の湿り具合を掛ける前の蒸発である。$1-f$ は、雪に覆われた部分 $f$ が空気と熱・水蒸気を交換しないことを表す。結露も同じく $1-f$ 倍にする。上向き長波 $\sigma T_L^4$ と下向き長波、地中への伝導 $K_{sd}(T_L-T_d)$ は変えない。積雪面の昇華は扱わないので、雪で覆われた部分の水蒸気交換は 0 である。大気最下層へ加える顕熱と水蒸気、地表から引く顕熱と潜熱には同じ $(1-f)H_L$ と $E_L$ を使う。格子の反射率は $w_L(\alpha_L+\Delta\alpha f)+w_{\rm ow}\alpha_o+w_i\alpha_i$ である。

$f$ は他の地表量と同じく前時刻の $\overline S^{n-1}$ から作る。

### 5.3 積雪と融雪の更新

陸の二層の式 $C_L\,dT_L/dt=F_L-K_{sd}(T_L-T_d)$ から得た傾向を $R_L$ とし、まず融雪を考えない候補

$$
T_L^\ast=T_L+\Delta\tau R_L,\qquad
S^\ast=S+\Delta\tau S_{\rm ls}
$$

を作る。$S_{\rm ls}$ は格子内で一様とみなし、陸部分にも同じ強度で降る。$T_L^\ast>T_{\rm gm}$ かつ $S^\ast>0$ なら、$T_{\rm gm}$ を超える分の熱で融かせる量と積雪量の小さい方

$$
m=\min\!\left(S^\ast,\;\frac{C_L(T_L^\ast-T_{\rm gm})}{L_f}\right)
$$

を融かし、それ以外なら $m=0$ とする。更新後は

$$
T_L^{\rm new}=T_L^\ast-\frac{L_f}{C_L}m,\qquad
S^{\rm new}=S^\ast-m,\qquad
M_s=\frac{m}{\Delta\tau}
$$

である。熱が足りれば $T_L^{\rm new}=T_{\rm gm}$ で雪が残り、雪が足りなければすべて融けて $T_L^{\rm new}>T_{\rm gm}$ になる。$T_L^\ast\le T_{\rm gm}$ なら降雪はそのまま積もる。傾向として $(T_L^{\rm new}-T_L)/\Delta\tau$ と $(S^{\rm new}-S)/\Delta\tau=S_{\rm ls}-M_s$ を返す。深い層 $T_d$ は変えない。

$C_L=2\times10^6\,\mathrm{J\,m^{-2}\,K^{-1}}$ なので、$T_L^\ast$ が $T_{\rm gm}$ を $1\,\mathrm K$ 超えるごとに水当量 $C_L/L_f\simeq6.0\,\mathrm{mm}$ の雪が融ける。

### 5.4 バケツへの融雪水

[バケツ](./bucket.md#水収支と流出)に入る液体の水を

$$
P_{\rm liq}=P_{\rm conv}+P_{\rm ls}-S_{\rm ls}+M_s
$$

とし、$W^\ast=W+\Delta\tau(P_{\rm liq}-E_L)$ から従来どおり流出 $R$ と $W^{\rm new}$ を決める。陸面の水収支は

$$
\frac{\partial}{\partial t}(W+S)=P_{\rm conv}+P_{\rm ls}-E_L-R
$$

で閉じる。積雪は流出しないので上限を持たない。

## 6. 海に落ちた雪

海部分に落ちた雪 $S_{\rm ls}$ はためず、直ちに融かす。融解熱は[海洋混合層](./slab-ocean.md)全体（海氷の下を含む）から取る。すなわち [Q flux](./q-flux.md#海洋の熱収支への加え方) と同じ位置で、海面積当たり

$$
Q\to Q-L_fS_{\rm ls}
$$

として[海氷の一回の物理更新](./sea-ice.md#51-開水面から海へ入る熱)に渡す。開水面率 $1-A$ は掛けない。融けた水は予報しない海洋貯水へ入る。

海氷があって $T_o=T_f$ の格子では、この熱は海水を凍らせ、氷の体積を $L_fS_{\rm ls}\Delta\tau/(\rho_iL_f)=S_{\rm ls}\Delta\tau/\rho_i$ だけ増やす（雪と氷の $L_f$ は同じ値）。すなわち海氷に降った雪は、同じ質量の氷に置き換わったのと熱的に同じになる。海氷がなければ $T_o$ が下がり、$T_f$ を下回れば新氷ができる。海氷が無効の設定では $T_o$ だけが下がる。

降水強度 $10\,\mathrm{mm\,day^{-1}}$ の雪は $L_fS_{\rm ls}\simeq38.7\,\mathrm{W\,m^{-2}}$ に相当し、$C_o$ の混合層を 1 日に約 $0.027\,\mathrm K$ 冷やす。

## 7. 水とエネルギーの収支

格子全体の地表のエネルギーを

$$
\mathcal E_{\rm sfc}=f_L(C_LT_L+C_dT_d-L_fS)+f_O\mathcal E_o
$$

とする（$\mathcal E_o$ は[海水と海氷のエネルギー](./sea-ice.md#6-局所エネルギー収支)）。積雪は液体の水より $L_f$ だけエネルギーが低いものとして数える。陸では降雪が $-L_fS$ を $L_fS_{\rm ls}\Delta\tau$ だけ減らし、融雪は $C_LT_L$ から $L_fm$ を取って $-L_fS$ に同じだけ戻すので、陸は正味 $L_fS_{\rm ls}\Delta\tau$ を失う。海は 6 節のとおり混合層から同じ量を失う。これが 4 節で大気が得る $L_fS_{\rm ls}$ と打ち消し合うので、1 回の物理更新で

$$
\sum_k\frac{\Delta p_k}{g}(c_pT_k+L_vq_k)+\mathcal E_{\rm sfc}
$$

の変化は、雪がないときと同じく放射と Q flux だけで決まる。水は大気の $\sum_kq_k\Delta p_k/g$ と陸の $f_L(W+S)$ の合計が、流出 $f_LR$ と海の $E-P$ を除いて保存する。RAW フィルターと切り上げを含む全積分の厳密保存は主張しない（8 節）。

## 8. 評価の順序と時間積分

[湿潤大気](../dynamics/moist.md#物理過程)の順序を次のように変える。

1. 前時刻の大気・陸・海・海氷・バケツ・積雪から蒸発候補を作る。陸の蒸発に $1-\overline f^{n-1}$ を掛ける。
2. 乾燥対流調節 → 湿潤対流調節 → 雪を含む大規模凝結（2, 3 節）→ 雲量の診断を暫定場で逐次に行う。雲量は融解の冷却を加えた暫定場から診断する。
3. 放射と地表タイルを評価する。陸のアルベドと顕熱に $\overline f^{n-1}$ を使い、陸の積雪・融雪（5.3 節）と海の融雪熱（6 節）をここで扱う。
4. 陸のバケツを雨と融雪水で更新する（5.4 節）。

従来はバケツを放射の前に更新していたが、同じ評価の融雪水を入れるために放射の後へ移した。バケツは放射の結果を使わず、放射もバケツの更新を使わないので、雪を無効にした計算ではこの入れ替えで結果は変わらない（ビット単位で一致する）。

$S$ は $W$ と同じく水平移流しない格子量で、スペクトル変換・重力波の陰解法・超粘性は適用しない。LeapFrog で

$$
S^{n+1}=\overline S^{n-1}+\Delta\tau\,(S_{\rm ls}-M_s)
$$

と進め、大気・地表温度と同じ RAW フィルターを掛ける。$\Delta\tau$ は通常 $2\Delta t$、初期化の 2 段階ではそれぞれ実際に進める幅で、どちらも $X^0$ から評価する。LeapFrog の候補と、RAW 後の $\overline S^n$、$S^{n+1}$ は最後に $S\ge0$ へ切り上げる。バケツと同じく、この切り上げは水収支に計上しない。$T_L$ が $T_{\rm gm}$ を超えたまま $S>0$ が残る状態を RAW が作っても、次の物理更新の 5.3 節で融雪されるので、RAW 後に相平衡への補正は行わない。

## 9. 初期条件・設定・互換性

陸海ケースは $S=0$ から始める。`dry_model_physics_config` の `snow_config` が次を保持し、既定値は設定型に一度だけ置く。

| フィールド | 既定値 | 意味 |
| --- | --- | --- |
| `enabled` | `.false.` | 陸海ケースで `.true.` |
| `formation_temperature` | $263.15\,\mathrm K$ | $T_{\rm sn}$ |
| `atmospheric_melting_temperature` | $278.15\,\mathrm K$ | $T_{\rm am}$ |
| `surface_melting_temperature` | $275\,\mathrm K$ | $T_{\rm gm}$ |
| `masking_water_equivalent` | $50\,\mathrm{kg\,m^{-2}}$ | $S_0$（$0.05\,\mathrm m$ 水当量） |
| `albedo_increase` | $0.4$ | $\Delta\alpha$ |
| `latent_heat_of_fusion` | $3.34\times10^5\,\mathrm{J\,kg^{-1}}$ | $L_f$ |
| `initial_water_equivalent` | $0$ | 陸を含む格子の初期積雪 $\mathrm{kg\,m^{-2}}$ |

有限な値、正の温度・$S_0$・$L_f$、$T_{\rm sn}\le T_{\rm am}$、$\Delta\alpha\ge0$、$\alpha_L+\Delta\alpha\le1$ を検証する。雪は水蒸気・大規模凝結・放射・タイル化した地表（陸海混合または海氷）を要求し、陸海混合では融雪水の受け皿としてバケツも要求する。雪が無効なのに 0 でない積雪を与えるとエラーにする。

雪を無効にした計算は従来とビット単位で一致する。陸海ケースは雪を有効にしたので結果が変わる。湿潤 aquaplanet は雪を無効のままにしたので変わらない。

## 10. 出力と診断

陸海ケースの出力に次を加える。フラックスは水当量で、ファイルでは $\mathrm{mm\,day^{-1}}$（$\mathrm{kg\,m^{-2}\,s^{-1}}$ の 86400 倍）である。

日平均（`daily_global.csv`）の列:

| 列 | 意味 |
| --- | --- |
| `snowfall_mm_day-1` | 地表に届いた大規模凝結の雪 $S_{\rm ls}$ の全球平均 |
| `atmospheric_snow_melt_mm_day-1` | 大気中で融けた雪 $M_{\rm atm}$ の全球平均 |
| `land_snowfall_mm_day-1` | $S_{\rm ls}$ の陸面積平均 |
| `land_snow_water_kg_m-2` | $S$ の陸面積平均 $\braket{S}_L$ |
| `land_snow_fraction` | $f$ の陸面積平均 $\braket{f}_L$ |
| `land_snow_melt_mm_day-1` | 地表の融雪 $M_s$ の陸面積平均 |
| `maximum_snow_budget_residual_mm_day-1` | $(S^{\rm new}-S)/\Delta\tau-(S_{\rm ls}-M_s)$ の最大絶対値。丸め誤差程度 |

月平均の格子場は `monthly_snow_water`（$S$）、`monthly_snow_fraction`（$f$ の時間平均）、`monthly_snowfall`（$S_{\rm ls}$）、`monthly_snow_melt`（$M_s$）、毎年 4/1 の瞬時値は再開用の `yearly_snow_water` と `yearly_snow_fraction` である。$S,f,M_s$ は陸面積当たり、$S_{\rm ls}$ は格子面積当たりの値で、純海洋格子の $S,f,M_s$ は 0 である。可視化では $S$ と $f$ を $f_L>0$ の格子だけに表示する。メタデータの `snow` にパラメータと方式を、`sea_ice.snow` に海氷上に雪をためないことを、`physics_order` に 8 節の順序を記録する。

## 11. 検証

`core/test/check_snow.f90` で次を確かめる。

| 条件 | 確認する結果 |
| --- | --- |
| $-10^\circ\mathrm C$ 未満の過飽和層、$-5^\circ\mathrm C$ の過飽和層、$5^\circ\mathrm C$ をわずかに超える未飽和層の気柱 | 雪の層は $L_s$ で加熱されて飽和、雨の層は $L_v$、融解層はちょうど $T_{\rm am}$ まで冷える。落下中の雪の質量保存 $S_{\rm ls}+M_{\rm atm}=\sum G_k$、気柱の $c_pT+L_vq$ の増加が $L_fS_{\rm ls}$ |
| 雪を無効にした同じ気柱 | 従来の雨だけの凝結とビット単位で一致 |
| 暖かく厚い層の上で雪が降る | 雪はすべて融け、融解層の冷却は $L_fM_{\rm atm}$ |
| 積雪率 | $f(0)=0$、$f(S_0)=1/2$、$f(3S_0)=3/4$ |
| 冷たい陸・$275\,\mathrm K$ を超える陸・薄い積雪 | 降雪は積もる／$T_L^{\rm new}=275\,\mathrm K$ で融雪量は超過熱 $/L_f$／雪がすべて融けて残りの熱は陸を暖める |
| タイルでの雪 | 反射が $(\alpha_L+0.4f)/\alpha_L$ 倍、最下層への顕熱が $-fH_L$ だけ変わる。大気・陸・積雪を合わせたエネルギーが冷たい陸と融雪中の陸の両方で閉じる。海では混合層が $L_fS_{\rm ls}$ を失い積雪は持たない |
| ケース設定 | 陸海ケースだけで雪が有効で、既定値が本仕様どおり |
| 陸海ケースの短い積分（初期積雪 $100\,\mathrm{kg\,m^{-2}}$） | $S\ge0$、$0\le f<1$、純海洋で $S=0$、バケツの範囲、積雪の収支残差が丸め誤差程度 |
