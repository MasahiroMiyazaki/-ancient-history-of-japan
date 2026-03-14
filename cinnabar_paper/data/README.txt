========================================================================
Ancient History of Japan: Cinnabar Resource and State Formation
- Data and Analysis Scripts -
========================================================================

本フォルダには、論文「辰砂資源の調達構造と古代国家形成モデルの再構築」
（宮崎政宏, 2026）の再現・検証に必要な全データおよび解析コードが格納されています。

【ライセンス】
本データおよびコードは CC BY 4.0 の下で公開されています。
引用・利用・改変を問わず、適切に出典を明記することで自由に利用可能です。

【フォルダ構成】
1. ソースデータ (CSVファイル)
   - mines_ancient_honshu.csv    : 古代鉱山27箇所の座標データ
   - mononobe_shrines_honshu.csv : 物部系神社17社の座標データ
   - hyakumeizan_shrines.csv     : 百名山神社18社の座標データ
   - ichinomiya_shrines_full.csv : 一之宮51社の座標データ

2. 解析スクリプト
   - analyze_ancient_japan.py
     上記4つのCSVファイルを読み込み、Haversine式による距離算出、
     Mann-Whitney U検定、効果量計算、および図表生成を一括して行うPythonコード。

3. 年代論検証用コード
   - oxcal_origin_code.txt
     既存の「春成モデル(2011)」を最新国際標準(IntCal20)で追試するためのOxCal実行コード。
   - oxcal_modify_code.txt
     信頼性の高い試料を優先的に採用した「代替ベイズモデル」用OxCal実行コード。
     Amodel=64.9%（国際基準クリア）を達成し、年代の客観的収束を証明する。

【再現手順】

--- 1. 統計分析（Python） ---
必要ライブラリ： pandas, numpy, scipy, matplotlib
インストール:
    pip install pandas numpy scipy matplotlib

実行:
    python analyze_ancient_japan.py

出力:
    images/ フォルダに、論文内で使用した「統計分布図」「距離閾値分析図」
    「詳細箱ひげ図」「距離ランキング図」「累積分布関数(CDF)」が自動生成されます。

--- 2. 年代論検証（OxCal） ---
Web版 OxCal (https://c14.arch.ox.ac.uk/oxcal.html) を使用します。

手順:
    1. OxCal を開き、[File] > [New] で新規プロジェクトを作成。
    2. [View] > [Text] を選択し、コードを貼り付け。
    3. [Run] を実行。
    4. 解析結果（[A=XX.X%] と表示される Amodel）を確認。

【論理的意図】
本公開の目的は、春成モデルの統計的不整合を隠蔽することではなく、検証可能な形で明示し、科学的基準を満たす新たな年代論を提示することにあります。
コードとデータは全て公開されているため、著者への問い合わせなしに再現・反証可能です。

========================================================================
出典：宮崎政宏 (2026) 「辰砂資源の調達構造と古代国家形成モデルの再構築」
URL: https://masahiromiyazaki.github.io/-ancient-history-of-japan/cinnabar_paper/index.html
========================================================================