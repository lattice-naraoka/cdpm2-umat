# 原著作表示と変更履歴

[English](THIRD_PARTY_NOTICES.md)

本配布物のソース、追加ツール、試験、説明資料はGNU Lesser General Public License version 2.1 or later（LGPL-2.1-or-later）で配布します。LICENSEの全文と原著作表示を保持してください。無保証です。バイナリを再配布する場合の義務を、このソースのみの配布方針で免除するものではありません。

## OOFEM

Copyright (C) 1993 - 2013 Borek Patzak

参照版はOOFEM commit 1c1dbd2b4b8c95dcc26878bf5009d44ed535e1a6です。構成則の原本（concretedpm2.C、concretedpm2.h）、関連ソース（floatmatrixf.h、structuralmaterial.C）、ライセンス（COPYING.LGPLv2.1）は以下の参照元で確認できます。

参照実装: https://github.com/oofem/oofem/tree/1c1dbd2b4b8c95dcc26878bf5009d44ed535e1a6

CDPM2原著: CDPM2: A damage-plasticity approach to modelling the failure of concrete, International Journal of Solids and Structures, volume 50, issue 24, 2013, pages 3805–3816。CDPM2固有の構成則の根拠であり、別の損傷モデルへ置換していません。

## 本配布版の変更

変更期間は2026年9月30日から10月12日です。OOFEM参照構成則をFortranへ移植し、Abaqusの応力・工学せん断・STATEV・回転・メッシュ形状のインターフェース、直接終点更新、ハイブリッド接線、粘性緩和、失敗時の状態保持と外部増分縮小を追加しました。失敗時のNewton初期値追加、収束解の照合、損傷増加量による受入れ制御を含みます。AbaqusのASCII形式の.filから形状を自動読込みし、形状表の事前準備は不要になりました。C3D8Iの8節点形状の読込みにも対応しました。追加の単体試験は同じLGPL条件です。自動読込みによる構成則更新の変更はありません。

OOFEM、Abaqus、HYPLASは別プロジェクトです。本配布版はOOFEMやSIMULIAが提供・承認したUMATではありません。非線形有限要素法の書籍とHYPLASは一般的な更新・収束・状態確定の参考として使用し、その書籍ページやプログラムは同梱していません。Abaqusのソフトウェア、マニュアル、ライブラリも同梱していません。

第三者の実験資料の権利は各権利者に属します。比較図は解析結果と公開図から読み取った数値に基づいて独自に作成したもので、原図や原測定データの再配布ではありません。

追加修正は、粘性緩和後の損傷増分判定、kappa増分に基づく局所収束判定、4元連立方程式のスケーリングと後退誤差確認、追加STATEVによる初期形状のリスタート保存です。従来の材料状態62変数は保持し、89変数入力で自動復元を利用します。構成則を別の損傷モデルへ置き換える変更ではありません。
