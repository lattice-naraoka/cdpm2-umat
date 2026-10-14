# CDPM2 UMAT

[English](README.md)

CDPM2をFortranで実装した、Abaqus/Standardの陰解法用UMATです。
CDPM2はGrasslらが開発したコンクリートの損傷塑性モデルです。本プロジェクトはOOFEMの参照実装に基づき、既存の構成則をUMATとして実装しています。

## はじめに

- [日本語マニュアル](docs/CDPM2_Manual_JA.pdf)：理論・検証・入力・使い方。
- [English manual](docs/CDPM2_Manual_EN.pdf)。
- [UMATソース](src/cdpm2.f90)。
- [使用例](examples/)。
- [切欠き梁の比較例](examples/notched_beam/)：CDP／CDPM2の6ケースと結果描画。

形状表の事前準備は不要です。初回はAbaqusのASCII形式`.fil`から初期形状を読み込みます。`*DEPVAR`は **89** を指定してください。従来の材料・仕事用62変数は変更せず、追加27変数に元の形状を保存するため、再開時は元のFILや別の形状ファイルが不要です。付属例には設定済みです。従来の62変数入力も実行できますが、その旧リスタートファイルには追加形状情報がありません。詳しくはマニュアルを参照してください。

全体の取得は **Code → Download ZIP**、またはgit cloneを使ってください。詳しい説明はマニュアルに一本化しています。利用に必要なソース、実行例、説明PDF、ライセンス・権利表示を収めています。

研究用の実装であり、SIMULIAやOOFEMの公式製品ではありません。使用前にマニュアルの制約を確認してください。

LGPL-2.1-or-later、無保証。[LICENSE](LICENSE)と[第三者の権利表示](THIRD_PARTY_NOTICES_JA.md)を参照してください。第三者の実験資料の権利は各権利者に属します。
