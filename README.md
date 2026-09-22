# BespokePlanet

架空の惑星をシミュレーションすることを目的とするプロジェクトです。Fortran のソルバーの計算結果を、Three.js 製のビジュアライザで 3D 球面と 2D 地図に表示できます。

## 依存

| 用途 | 必要なもの |
| --- | --- |
| ソルバー（`core/`） | Fortran コンパイラ（gfortran など、OpenMP 対応）、[fpm](https://fpm.fortran-lang.org/)、FFTW3、LAPACK、BLAS、pkg-config |
| ビジュアライザ（`viz/`） | Python 3（外部パッケージ不要）、Node.js / npm（Three.js の取得に使用） |
| タスクランナー | [just](https://github.com/casey/just) |
| `scripts/*.py`（座標・地形の確認図、地球地形の前処理） | numpy、scipy、matplotlib（スクリプトにより一部） |
| `scripts/analyze_*.R` | R（標準パッケージのみ） |

macOS（Homebrew）の場合のインストール例:

```sh
brew install gcc fpm fftw lapack openblas pkg-config just node
```

## 実行方法

以下のコマンドは、リポジトリ内のどのディレクトリからでも実行できます。

```sh
# テストを実行
just test

# モデルを実行（shallow-water / barotropic / dry / held-suarez / radiation / slab-ocean / moist / land / land-t63 / land-earth / land-earth-t63 / all、省略時は dry）
just run dry

# Held–Suarez 強制を 200 日間実行（5 日ごとに出力）
just run held-suarez

# 日変化・季節変化を含む放射ケースを5年間実行（日次・月平均・毎年4月1日の瞬時値を出力）
just run radiation

# 地軸傾斜0・深さ30 mのslab oceanを持つ放射ケースを5年間実行
just run slab-ocean

# slab ocean ケースに水蒸気（蒸発・湿潤対流調節・大規模凝結・水蒸気に依存する長波放射・雲量診断による短波の反射）を加え、
# 地軸傾斜を地球の値に戻した湿潤ケースを T31 で5年間実行
just run moist

# 湿潤ケースに大陸・山脈と陸海混合の地表を加えた陸海ケースを T31 / T63 で5年間実行
just run land
just run land-t63

# 陸海ケースの地形を ETOPO 2022 から作った地球の地形に置き換えて T31 / T63 で5年間実行
just run land-earth
just run land-earth-t63

# ビジュアライザを起動し、http://127.0.0.1:8000 を開く
just viz
```

計算結果の解析は `scripts/` の R スクリプトで行う。いずれも第1引数にケースの出力ディレクトリ、第2引数に解析結果の出力先を取り（省略時は `output/<ケース>/analysis/`）、日平均の時系列と最終年の月平均から図と CSV を書く。陸海ケースでは雨温図・ケッペンの気候区分・帯状平均・質量流線関数などを作る。

```sh
# 湿潤 aquaplanet ケース（5 年の全球平均の時系列、最終年の帯状平均と質量流線関数）
Rscript scripts/analyze_moist_slab_ocean.R

# 陸海ケース（解析的な大陸）
Rscript scripts/analyze_moist_land_sea.R

# 地球地形の陸海ケース（第2引数で出力先、第3引数で比較対象のケースを変えられる）
Rscript scripts/analyze_moist_land_sea_earth.R output/moist_land_sea_earth_t31
```

ソルバーはスペクトル変換の鉛直層ループと物理過程の格子列ループを OpenMP で並列化している。スレッド数は既定でコア数で、`OMP_NUM_THREADS=4 just run moist` のように環境変数で変えられる。スレッド数を変えても結果はビット単位で同一である。

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
- [湿潤大気](./docs/dynamics/moist.md)
- [Jablonowski–Williamson の初期状態](./docs/dynamics/Jablonowski-Williamson.md)
- [地形](./docs/dynamics/topography.md)
- [地球の地形](./docs/dynamics/earth-topography.md)

### 傾向

- [物理過程を評価する時刻](./docs/tendency/physics-time-level.md)
- [長波放射](./docs/tendency/longwave-radiation.md)
- [短波放射](./docs/tendency/shortwave-radiation.md)
- [オゾン](./docs/tendency/ozone.md)
- [地面](./docs/tendency/ground.md)
- [Slab ocean](./docs/tendency/slab-ocean.md)
- [陸と海の混合](./docs/tendency/land-sea-surface.md)
- [上層の Rayleigh 摩擦](./docs/tendency/upper-rayleigh-friction.md)
- [Held–Suarez 強制](./docs/tendency/Held-Suarez.md)
- [Betts–Miller 型の乾燥対流調節](./docs/tendency/dry-convective-adjustment.md)
- [飽和比湿](./docs/tendency/saturation-specific-humidity.md)
- [蒸発](./docs/tendency/evaporation.md)
- [陸面のバケツモデル](./docs/tendency/bucket.md)
- [大規模凝結](./docs/tendency/large-scale-condensation.md)
- [簡略化 Betts–Miller 型の湿潤対流調節](./docs/tendency/moist-convective-adjustment.md)
- [雲](./docs/tendency/cloud.md)

### 暦

- [暦と軌道](./docs/calendar.md)

### ケース定義

- [乾燥大気のケース](./docs/cases/dry.md)
- [Held–Suarez ケース](./docs/cases/held-suarez.md)
- [放射ケース](./docs/cases/radiation.md)
- [Slab ocean ケース](./docs/cases/slab-ocean.md)
- [湿潤ケース](./docs/cases/moist.md)
- [陸海ケース](./docs/cases/land-sea.md)
- [地球地形の陸海ケース](./docs/cases/land-sea-earth.md)

## ライセンス

[MIT License](./LICENSE)
