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

# モデルを実行（shallow-water / barotropic / dry / held-suarez / all、省略時は dry）
just run dry

# Held–Suarez 強制を 200 日間実行（5 日ごとに出力）
just run held-suarez

# ビジュアライザを起動し、http://127.0.0.1:8000 を開く
just viz
```

計算結果は `output/` に保存されます。ビジュアライザについて詳しくは [viz/README.md](./viz/README.md) を参照してください。

## ドキュメント

- [球面上の微分演算子](./docs/spherical-derivation.md)
- [Octahedral Gaussian Grid](./docs/octahedral-gaussian-grid.md)
- [順圧渦度方程式](./docs/barotropic-vorticity-equation.md)
- [浅水方程式](./docs/shallow-water-equation.md)
- [hybrid-sigma 座標](./docs/hybrid-sigma.md)
- [乾燥大気](./docs/dry.md)
- [Jablonowski–Williamson の初期状態](./docs/Jablonowski-Williamson.md)
- [Held–Suarez 強制](./docs/Held-Suarez.md)

## ライセンス

[MIT License](./LICENSE)
