# BespokePlanet

架空の惑星をシミュレーションすることを目的とするプロジェクトです。Fortran のソルバーの計算結果を、Three.js 製のビジュアライザで 3D 球面と 2D 地図に表示できます。

## 依存

| 用途 | 必要なもの |
| --- | --- |
| ソルバー（`core/`） | Fortran コンパイラ（gfortran など）、[fpm](https://fpm.fortran-lang.org/)、FFTW3、LAPACK、BLAS、pkg-config |
| ビジュアライザ（`viz/`） | Python 3（外部パッケージ不要）、Node.js / npm（Three.js の取得に使用） |
| タスクランナー | [just](https://github.com/casey/just) |
| `scripts/check_hybrid_sigma.py` | matplotlib |

macOS（Homebrew）の場合のインストール例:

```sh
brew install gcc fpm fftw lapack openblas pkg-config just node
```

## 実行方法

以下のコマンドは、リポジトリ内のどのディレクトリからでも実行できます。

```sh
# テストを実行
just test

# モデルを実行（shallow-water / barotropic / dry / held-suarez / radiation / slab-ocean / all、省略時は dry）
just run dry

# Held–Suarez 強制を 200 日間実行（5 日ごとに出力）
just run held-suarez

# 日変化・季節変化を含む放射ケースを5年間実行（日次・月平均・毎年4月1日の瞬時値を出力）
just run radiation

# 地軸傾斜0・深さ30 mのslab oceanを持つ放射ケースを5年間実行
just run slab-ocean

# ビジュアライザを起動し、http://127.0.0.1:8000 を開く
just viz
```

計算結果は `output/` に保存されます。ビジュアライザについて詳しくは [viz/README.md](./viz/README.md) を参照してください。

## ドキュメント

力学コア（`docs/dynamics/`）、傾向を与える物理過程（`docs/tendency/`）、それらを使う計算ケースの定義（`docs/cases/`）に分けて置いている。

### 力学

- [球面上の微分演算子](./docs/dynamics/spherical-derivation.md)
- [Octahedral Gaussian Grid](./docs/dynamics/octahedral-gaussian-grid.md)
- [順圧渦度方程式](./docs/dynamics/barotropic-vorticity-equation.md)
- [浅水方程式](./docs/dynamics/shallow-water-equation.md)
- [hybrid-sigma 座標](./docs/dynamics/hybrid-sigma.md)
- [乾燥大気](./docs/dynamics/dry.md)
- [Jablonowski–Williamson の初期状態](./docs/dynamics/Jablonowski-Williamson.md)

### 傾向

- [物理過程を評価する時刻](./docs/tendency/physics-time-level.md)
- [長波放射](./docs/tendency/longwave-radiation.md)
- [短波放射](./docs/tendency/shortwave-radiation.md)
- [オゾン](./docs/tendency/ozone.md)
- [地面](./docs/tendency/ground.md)
- [Slab ocean](./docs/tendency/slab-ocean.md)
- [上層の Rayleigh 摩擦](./docs/tendency/upper-rayleigh-friction.md)
- [Held–Suarez 強制](./docs/tendency/Held-Suarez.md)
- [Betts–Miller 型の乾燥対流調節](./docs/tendency/dry-convective-adjustment.md)

### 暦

- [暦と軌道](./docs/calendar.md)

### ケース定義

- [乾燥大気のケース](./docs/cases/dry.md)
- [Held–Suarez ケース](./docs/cases/held-suarez.md)
- [放射ケース](./docs/cases/radiation.md)
- [Slab ocean ケース](./docs/cases/slab-ocean.md)

## ライセンス

[MIT License](./LICENSE)
