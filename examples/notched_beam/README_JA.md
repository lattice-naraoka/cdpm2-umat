# 切欠き梁：CDPとCDPM2

[English](README.md)

70（`coarse`）・280（`medium`）・1120（`fine`）要素のC3D8R、Riks法、粘性なしのAbaqus/Standard比較例、計6ケースです。各ペアのメッシュ・境界条件・ステップ制御は同じです。CDPM2では材料を置き換え、必要なアワーグラス剛性を追加し、材料出力をSDVに変更しています。ASCII形式のFIL出力指定でUMATへ初期メッシュ形状を渡します。

## 実行

[マニュアル](../../docs/CDPM2_Manual_JA.pdf)に従い、`../../src/cdpm2.f90`をリポジトリ直下の`build/`でビルドします。このフォルダで、一つずつ実行してください。

```text
abaqus job=cdp_coarse cpus=1 mp_mode=threads interactive output_precision=full
abaqus job=cdpm2_coarse cpus=1 mp_mode=threads interactive output_precision=full
```

他のメッシュは`coarse`を`medium`または`fine`へ置き換えます。形状表やメッシュ準備ツールは不要です。CDPM2入力を変更するときも、付属のASCII形式FIL指定、メッシュ全域のFILE出力指定、`*DEPVAR`の **89** を保持してください。同名ジョブを同時実行しないでください。元の形状は追加状態変数へ保存します。静的リスタートはせん断例で扱い、このRiks例のリスタートは未検証です。

完了した各ODBをAbaqus Pythonで抽出し、通常のPython（NumPy・Matplotlib）で描画します。

```text
abaqus python extract_curve.py cdp_coarse.odb
abaqus python extract_curve.py cdpm2_coarse.odb
python plot_comparison.py
```

6ケースすべてを抽出してから描画してください。新しく抽出したCDP／CDPM2の曲線を描き、`results/comparison.png`へ保存します。既存の検証結果と実験との比較はマニュアルに載せています。

## 条件と出典

この例題入力はAbaqus 2026 Benchmarks、*Notched unreinforced concrete beam under 3-point bending*、公式`notchedconcbeam_3d_{coarse,medium,fine}_std.inp`を参考に作成しています。CDPの引張軟化・引張損傷を公式2DのTS-2表へ変更しています。元の3D一線形軟化ケースそのものの再現ではありません。メッシュは明示的な節点・要素・集合で記述しています。UMATは入力キーワードではなく、Abaqusが生成したFIL形状情報を読みます。物理的な形状・節点接続は変えていません。

単位はN・m・Pa。E=30 GPa、nu=0.2、fc=30 MPa、ft=3.33 MPa。TS-2の開口・応力点は（0 mm, 3.33 MPa）、（0.025 mm, 0.833 MPa）、（0.198 mm, 0.0333 MPa）。有限表区間の面積0.12697 N/mmにCDPM2の二直線軟化面積を揃えていますが、CDPの1%残留尾部とは同一ではありません。モデル幅は20・10・5 mmで、曲線は全梁・幅50 mm相当の荷重へ換算します。Riksステップ時間は物理的な載荷時間ではありません。その他の材料入力と制約はマニュアルを参照してください。

実験：Petersson（1981）、*Crack Growth and Development of Fracture Zones in Plain Concrete and Similar Materials*、Report TVBM-1006、Lund Institute of Technology。直接の図読取り元は上記Abaqusベンチマークの図11です。点は概算読取りで、実験の原始記録ではありません。Gf=115・137 N/mは実験系列のラベルで、今回の解析入力値ではありません。第三者資料の権利は元の権利者に属し、本プロジェクトのライセンスで再許諾するものではありません。
