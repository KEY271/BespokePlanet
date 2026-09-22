# Held–Suarez ケース

[乾燥大気](../dynamics/dry.md)の計算手法に [Held–Suarez 強制](../tendency/Held-Suarez.md)を加えたケースの定義である。温度には同文書の Newton 緩和を、外力 $\bm{F}_k$ には同文書の Rayleigh 摩擦を与える。ほかの物理過程は用いない。

## 初期状態

初期状態は $u=v=0$ で静止した等温大気とする。ただし対称性を破るために小さな擾乱を加え

$$
T=264\,\mathrm{K}+0.1\,\mathrm{K}\cos\varphi\sin\lambda
$$

とする。圧力は $p_s=10^5\,\mathrm{Pa}$ とし、ジオポテンシャルは $\Phi_s=0$ の平坦な地形を与える。

## 時間

$200$ 日間シミュレーションする。出力は $5$ 日ごとに保存する。タイムステップは $\Delta t=1200\,\mathrm{s}$ で他の計算と同じとする。
