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
1. dataset.csv
   - 統計分析用データセット（神社・鉱山座標・距離データ）

2. analyze_ancient_japan.py
   - 本論文の主要な統計分析（Mann-Whitney U検定、効果量算出、図表生成）
     および統計値の算出を行うPythonコード。

3. oxcal_origin_code.txt
   - 既存の「春成モデル(2011)」を最新国際標準(IntCal20)で追試するためのOxCalコード。
     統計的不整合（Amodel < 60%）の確認用。

4. oxcal_modify_code.txt
   - 本論文の「代替ベイズモデル」用OxCalコード。
     信頼性の高い試料を優先・適正な層位順序を設定することで、
     Amodel=64.9%（基準クリア）を達成し、年代の客観的収束を証明する。

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
    2. [View] > [Text] を選択し、空のテキストエリアを表示。
    3. oxcal_origin_code.txt または oxcal_modify_code.txt の内容を貼り付け。
    4. [Run] を実行。
    5. 解析結果（[A=XX.X%] と表示される Amodel）を確認。

【論理的意図】
本公開は、論文の主張に対する第三者による独立した検証・反証を促進するためにあります。
春成モデルの統計的不整合は、科学的再検証の必要性を示す客観的事実です。
コードとデータは全て公開されているため、著者への問い合わせなしに再現・反証可能です。

========================================================================
出典：宮崎政宏 (2026) 「辰砂資源の調達構造と古代国家形成モデルの再構築」
URL: https://masahiromiyazaki.github.io/-ancient-history-of-japan/cinnabar_paper/index.html
========================================================================