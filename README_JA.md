# CDPM2 UMAT

[English](README.md)

GrasslらのCDPM2をAbaqus/Standardの陰解法で利用するFortran UMATです。OOFEMの参照実装に基づきます。
**研究用プレビューです。実務検証済みの安定版や、SIMULIA・OOFEMの公式製品ではありません。**

## 資料

- [マニュアル](docs/CDPM2_Manual_JA.pdf)：現在の材料モデル、数値更新法、入力・使い方。
- [検証ガイド](docs/CDPM2_Verification_JA.pdf)：現在の配布ソースの試験条件・結果・制約。
- [解析例題集](examples/pullout/CDPM2_Examples_JA.pdf)：目次、例題1、例題2。
- [更新履歴（日本語版・英語版）](docs/CDPM2_Change_History_JA.pdf)：変更点を言語別に集約。

文書版は2026-11-09です。例題集は目次が1ページ、例題1が2～6ページ、例題2が7～11ページです。

## ソースと検証範囲

- [UMATソース](src/cdpm2.f90)
- [検証詳細](VALIDATION.md)・[機械可読の証拠](evidence/)
- [例題1の入力](examples/pullout/)・[例題2の入力](examples/pullout/case2/)・[例題2の説明](CASE2.md)
- [切欠き梁の入力](examples/notched_beam/)・[1要素試験の入力と参照データ](verification/)

312条件の基準版との回帰比較を実施し、36解析が正常終了、最大応力差は3.81149e-9 MPaでした。
これは新たなOOFEM比較や実験妥当性の確認ではありません。保存局所試行1094件は724収束・370棄却、
指定15幅の接線検査は700/724点で基準達成です。同一ソースのリスタートは開始フレームの反力差が
当初基準を超えており、全量一致は確認できていません。詳細は検証ガイドを参照してください。
既存の保存例題・参照データを現在の配布ソースの再実行結果とは扱いません。

## 使用方法

正規のAbaqusと対応Fortranコンパイラを用意し、マニュアルの手順で使用してください。
初期形状はASCII形式のFILから読み込み、基本入力は24定数・89状態変数です。
異なるソース版を混在させたリスタートの互換性は保証しません。
全体の取得は **Code → Download ZIP** またはgit cloneを使用できます。

独立スペクトル試験：`python tests/spectral_reference.py --compiler gfortran`
単独試験用スタブをAbaqusへリンクしないでください。

LGPL-2.1-or-later、無保証。[LICENSE](LICENSE)と[第三者の権利表示](THIRD_PARTY_NOTICES_JA.md)を保持してください。
