# 鉄筋引抜き解析

[English](README.md)

CDPM2 UMATと組込みCDPを比較するAbaqus/Standardの1/4モデル例題です。メッシュ、鉄筋、支持、接触、載荷条件を共通とし、CDPの引張・圧縮損傷上限は0.99としています。粘性はCDPが1e-4、CDPM2が1e-3です。

- [例題の説明と結果](Pullout_Example_JA.pdf)
- [英語版](Pullout_Example_EN.pdf)

## 実行

正規のAbaqusと対応Fortranコンパイラを用意し、このフォルダで実行します。

```text
abaqus job=cdp interactive
abaqus job=cdpm2 user=../../src/cdpm2.f90 interactive
```

CDPM2はリポジトリ共通のソースを使用します。付属のASCII FIL出力指定と `*DEPVAR` の89を維持してください。共通メッシュは `mesh.inp`、CDP材料表は `material.inp` です。UMATの導入方法と材料モデルは[マニュアル](../../docs/CDPM2_Manual_JA.pdf)を参照してください。

## 結果

両ケースともAbaqus/Standard 2026 HF0で載荷変位2 mmまで完了しています。全試験体換算の最大荷重はCDPが178.85 kN、CDPM2が139.95 kNです。PDFには荷重–すべり比較と最大主対数ひずみコンターを掲載しています。

解析後、Abaqus PythonでODBを抽出し、通常のPython（NumPy・Matplotlib）で描画します。

```text
abaqus python extract_curve.py cdp
abaqus python extract_curve.py cdpm2
python plot_results.py
```

抽出結果は `cdp_curve.csv`、`cdpm2_curve.csv`、描画結果は `comparison.png` です。描画スクリプトは抽出した解析曲線を表示します。荷重は1/4モデルの載荷端反力を4倍し、すべりは非載荷端鉄筋平均変位からコンクリートの (50,0,0)、(0,50,0) mmの2点平均変位を引いて求めます。

## モデルと出典

単位はN・mm・MPaです。C3D8Rの86,152要素で、200 mm立方体の1/4対称モデル、公称径20 mmの鉄筋、付着長100 mmを扱います。リブはリング状とし、コンクリート上面全域を鉛直拘束します。接触は線形ペナルティ（剛性倍率0.005）、摩擦係数0.5です。形状・材料値・出典の詳細はPDFに記載しています。

参考資料：S. Seok, *Finite Element Modeling of Bond-Zone Behavior in Reinforced Concrete*, Purdue University, 2019, Chapter 4（[書誌](https://docs.lib.purdue.edu/dissertations/AAI30502919/)）、G. Metelli・G. A. Plizzari, *Influence of the relative rib area on bond behaviour*, 2014（[DOI](https://doi.org/10.1680/macr.13.00198)）。リング状リブと支持条件は本例題のモデル化条件です。PDFの実験曲線はSeok Fig.4.8からの概略読取りで、原図・原測定データは収録していません。
