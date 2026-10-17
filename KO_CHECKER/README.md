# KO_Checker — RNA-seq KO evidence explorer

公共RNA-seqの**SRXをブラウザに入力し、データ取得からHISAT2マッピング、RNA上のKO評価までを実行する研究用プロトタイプ**です。マッピング済みの自前BAMを使うCLIも含みます。

**DNA上のノックアウトやタンパク質消失を確定する検査ではありません。** `KO_like`はRNA上の暫定的な証拠を示します。

## できること

- WT群・KO群に、それぞれ複数のSRX / ERX / DRXを指定。
- 各実験に属する公開済みの全runを取得し、SRA → FASTQ → HISAT2 → BAMを自動実行。
- SRX内の複数runを統合し、既存の`rnako.py`で変異・スプライス接合・予測タンパク質・NMDを解析。
- 完了済みの参照、索引、FASTQ、BAMを再利用。
- 画面から解析を中断・再開。参照ファイルは配信元が対応していれば途中から取得を再開。
- 参照の取得・展開・索引構築、SRA取得・検査・FASTQ変換を分けて表示し、取得容量・速度・実行時間を確認。
- 工程・ログ・履歴・判定根拠の表、CPM・WT基準比・変異支持率のグラフを表示。
- CSV、実行情報、ログをダウンロード。
- エクソンごとの深度・読み取りの欠落・WTとの相対深度を表示し、深度不足と局所的な読み取り低下の疑いを区別（DNA上の欠損を確定するものではありません）。

## 起動する

Docker DesktopのLinuxコンテナ環境、またはLinux上のDocker EngineとComposeを使用します。初回の参照・RNA-seq取得にはインターネット接続が必要です。非圧縮FASTQ、一時領域、ゲノム索引を保存するため、データ量に応じたメモリと、ダウンロード量より多くの空き容量を用意してください。

### Windows：ダブルクリックで起動

Docker Desktopをインストールしたうえで、このフォルダの **[start.bat](start.bat) をダブルクリック**します。Docker Desktopが停止していれば起動を試み、アプリのビルド・起動後にブラウザで **http://localhost:8080** を自動で開きます。初回のビルドには数分かかる場合があります。

次回も同じ`start.bat`で起動できます。起動後にターミナルを閉じてもアプリと解析は動き続けます。失敗した場合はウィンドウにエラーを表示して停止します。Docker DesktopはLinuxコンテナモードで使用してください。

### ターミナルで起動（Windows / Linux）

1. リポジトリをcloneするか、GitHubの「Code → Download ZIP」で取得・展開します。
2. `Dockerfile`と`compose.yaml`があるフォルダでターミナルを開きます。
3. 次を実行します。

```bash
docker compose up -d --build
```

4. **http://localhost:8080** を開きます。
5. 対象遺伝子、参照ゲノム、WT・KOのSRXを指定し、「解析を開始する」を押します。

Windowsでは、エクスプローラーでプロジェクトのフォルダを開き、アドレス欄に`powershell`と入力すると、その場所でPowerShellを開けます。

| 入力 | 内容 |
| --- | --- |
| 対象遺伝子 | 遺伝子名（例：`Phf3`）またはGTFに記載されたgene ID |
| 参照 | マウス`mm10` / `mm39`、ヒト`hg38` |
| WT・KO | 各群1〜30件のSRX / ERX / DRX。改行・カンマ・空白区切り |

群ラベルは利用者が指定します。登録情報からWT・KOを自動で決める機能はありません。**SRR / ERR / DRRはrun番号のため、Web入力では受け付けません。** 実験番号から配下の全runを自動解決し、統合します。別々のSRXが独立した生物学的反復かどうかは、実験情報で確認してください。

初回の参照取得とHISAT2索引作成は長時間かかる場合があります。上段の%は工程の目安です。その下に参照取得量・速度とファイル単位の残り時間の目安を表示します。参照の展開・FASTA索引・遺伝子注釈確認・HISAT2索引を別の工程名で表示します。索引構築中のログに出るbucketの%は、その部分処理の進捗であり全体の転送率ではありません。

## 処理の流れ

```text
WT / KOのSRX
  → ENAで配下の全runを検索・生物種とRNA-Seq構成を確認
  → GENCODE参照とHISAT2索引を準備
  → runごとにprefetch → vdb-validate → fasterq-dump
  → HISAT2 → samtools sort / index
  → 複数runをSRX単位のBAMへ統合
  → rnako.pyでRNA上の証拠を比較
  → 判定・表・グラフ・CSV
```

ENAの一時的なHTTPエラーや通信エラーは自動再試行し、通常の取得APIで取得できない場合は検索APIへ切り替えます。複数runのBAMを統合する際はリード名にrunを区別する接頭辞を付け、異なるrunの断片が誤って1つに数えられるのを防ぎます。

参照はmm10＝GENCODE M25、mm39＝M39、hg38＝human release 47です。HISAT2は`--dta`、非ストランド指定で実行します。アダプタ除去や包括的なFASTQ QCは自動実行しません。

## 結果の読み方

| 判定 | 意味 |
| --- | --- |
| `KO_like` | 比較検体に対して濃縮された有害候補変異、または条件を満たす未注釈エクソンスキップがある |
| `WT_like` | 他検体のKO候補位置などを十分に観測でき、低い変異率や正常接合の証拠がある |
| `undetermined` | カバレッジ不足、比較検体不足、明確な証拠がないなどで判断できない |

「異常が見つからない」だけではWT_likeになりません。発現量の低下だけでKO_likeにすることもありません。WT基準比は、WTラベルかつWT_likeと判定された検体の平均CPMを基準にします。基準検体が揃わない場合は空欄です。

出力と判定閾値、座標、解釈上の制限は[解析仕様](docs/analysis.md)を参照してください。

## 保存場所と再利用

| ホスト側のフォルダ | 内容 |
| --- | --- |
| `references/` | 参照FASTA・GTF・HISAT2索引 |
| `workdata/metadata/` | SRXとrunの対応。24時間キャッシュ |
| `workdata/sra/`、`workdata/fastq/` | 取得したSRAとFASTQ |
| `workdata/bam/`、`workdata/samples/` | runごとのBAMとSRX単位の統合BAM |
| `results/web/<解析ID>/` | ジョブ状態、入力、実行情報、ログ、`output/`配下の結果 |
| `inputs/` | CLIで使う検体CSVなど |

これらのフォルダはDocker起動・解析時に作成されます。ローカルデータは`.gitignore`と`.dockerignore`でソースから除外します。完了記録、ファイルのサイズ・更新日時、参照やツールバージョンなどの条件を照合して再利用します。ファイル内容の完全なハッシュ検証ではないため、キャッシュは手動で改変しないでください。

ブラウザやチャットを閉じても解析は続きます。DockerとホストPCは稼働を継続する必要があります。PCのスリープ中は処理を継続できません。

### CPUスレッド数を指定する

新しい解析の「CPUスレッド数（オプション）」を開き、1〜256の整数を入力します。空欄の場合はサーバーの既定値（通常は認識した論理CPU数と4の小さい方）を使用します。設定は解析ごとに保存され、中断後の再開・同じ条件での再実行でも引き継ぎます。実際に使用する値は解析タイトル、ログ、`provenance.json`で確認できます。HISAT2・samtools・FASTQ変換・索引構築など、対応するツールに指定します。同時に動くツールがあるため、処理全体のCPU使用量の上限を保証する値ではありません。

### 画面から中断・再開する

実行中の解析で「中断する」を押し、表示が**「中断中」**になるまで待ってください。「中断処理中」はまだ処理終了前です。実行中の外部コマンドは子プロセスを含めて終了させます。ネットワークの応答待ち、GTF読み込み、FASTA索引作成などでは中断完了まで数十秒かかる場合があります。

「再開する」を押すと**同じ解析ID**で再開し、ログを追記します。前回までに完成した参照・索引・FASTQ・BAMは再利用します。待機中のジョブも中断でき、中断したジョブの次に待機中の解析があれば処理が移ります。

新規解析で「参照・索引が完成したら中断する」を選ぶと、参照準備が完了した時点で自動的に中断します。その後の「再開」では再度同じ場所で中断せず、データ取得へ進みます。外部の監視スクリプトは不要で、アプリ自体は稼働を続けます。

| 中断した工程 | 再開時の動作 |
| --- | --- |
| 参照FASTA/GTFの取得 | 保存した部分データとETag等を使い、HTTP Range対応・同一内容を確認できれば続きから取得。非対応・内容変更時は安全に取り直す |
| SRAの取得 | SRA Toolkitのprefetchによる再開・取得済みデータ再利用 |
| 参照の展開・索引構築 | 完成していれば再利用。未完了の工程は最初からやり直す |
| FASTQ変換・マッピング | 完成済みrunは再利用。実行途中の工程はやり直す |
| KO評価・比較表作成 | BAMを再利用し、評価工程を再実行 |

参照の部分取得は`.partial`と`.partial.json`をセットで保持します。配信元がRangeを無視した応答を返した場合に部分データへ連結しない、Content-Rangeが不一致なら取り直すなどの検証を行います。旧版で作成した部分ファイルに対応する検証情報がない場合は、安全のため先頭から取得します。計算の任意の命令位置から復帰する機能ではありません。

PCやDockerも停止したい場合は、まず画面で「中断中」になったことを確認してから実行します。他の待機ジョブが動き始めていないかも確認してください。

```bash
docker compose stop
```

次回の起動:

```bash
docker compose up -d
```

サーバー停止時に実行中・待機中だった解析は、次回起動時に「中断中」として復元します。「再開する」で再開してください。以前の失敗済みジョブは従来どおり「同じ条件で再実行」で新しい解析を作成します。データ自動削除、複数サーバーによる同一キャッシュの同時利用は未対応です。

## 自前BAMをCLIで解析する

自前FASTQ/BAMのWeb入力は未実装です。座標順・インデックス付きBAMは既存CLIで解析できます。`example_samples.csv`は形式の例で、検体番号はダミーです。

ホスト側の`Analysis/<ID>/<ID>.bam`を含むディレクトリを`.env`の`RNAKO_DATA_DIR`へ絶対パスで設定します。[.env.example](.env.example)を参照してください。既存コンテナへのマウント変更は再作成を伴うため、実行中の解析が終わってから設定してください。

```bash
docker compose up -d
docker exec rnako python /app/rnako.py --gene Phf3 --assembly mm10 --samples /inputs/samples.csv --out /results/my_analysis
```

`inputs/samples.csv`の列は`srx_id,metadata_label,bam_path`です。`bam_path`はコンテナ内のパス（例：`/data/Analysis/SRX1/SRX1.bam`）で指定します。BAMに使用した参照と同じゲノムを使ってください。独自参照の場合は`--fasta`と`--gtf`を両方指定します。染色体名の別名と長さを照合しますが、長さだけでは塩基配列の完全一致を保証できません。

## 開発・テスト

```bash
docker compose run --rm rnako python -m unittest discover -s /app -p 'test_*.py' -v
```

テストには人工BAMによる統合、外部API障害への再試行、prefetchの保存先の違い、中断・再開、HTTP Rangeや配信元の内容変更、子プロセス終了を検査するケースを含みます。取得再開はローカルのHTTPサーバーで実際に通信して検証します。外部APIや大型データ取得はテスト内でモックします。

中断・再開とエクソン比較を含む統合版の検証内容はCHANGELOG.mdを参照してください。

公開RNA-seqのDaxx 2検体で、取得からマッピング・CSV出力までの完走を確認しています。両検体のRNA判定は`undetermined`であり、これは分類精度や他の検体でのKO検出を保証する検証ではありません。

Dockerを使わずに画面だけ試す場合は`python webapp.py`（Windowsでは`py webapp.py`）を実行します。実解析にはpysam、HISAT2、samtools、SRA Toolkitが必要です。Webサーバー自体はPython標準ライブラリで実装しています。推奨の実行環境はDocker内のPython 3.12です。

GitHub ActionsはDockerビルドとテストを実行します。GitHubへの追加方法とソースのみのZIP作成は[GitHub共有手順](docs/github-sharing.md)を参照してください。

## 制約と公開範囲

通常のbulk RNA-seq向けの研究用プロトタイプです。single-cell/UMI、特殊なread構成、非公開データ、複雑な構造変異、離れた変異のフェージング、群間の統計検定は対象外です。メタデータのRNA-Seq分類だけではbulkであることは保証できません。

ログインや利用者間のアクセス分離はありません。Composeは`127.0.0.1:8080`に限定して公開します。`localhost`のリンクを他人に送っても、その人があなたの解析へ接続できるわけではありません。学内やインターネットへ公開する場合は、認証やアクセス制限、計算資源の管理を別途実装してください。

## ライセンス・参照資料

この配布物にはライセンスをまだ設定していません。利用・改変・再配布の条件はリポジトリ管理者に確認してください。依存ツールと参照データには、それぞれの利用条件が適用されます。

- [NCBI SRA Toolkit：取得とFASTQ変換](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump)
- [ENA file reports](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html)
- [ENA検索API](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/advanced-search.html)
- [HISAT2マニュアル](https://daehwankimlab.github.io/hisat2/manual/)
- [GENCODE mm10 M25](https://www.gencodegenes.org/mouse/release_M25.html)、[mm39 M39](https://www.gencodegenes.org/mouse/release_M39.html)、[hg38 release 47](https://www.gencodegenes.org/human/release_47.html)

変更履歴は[CHANGELOG.md](CHANGELOG.md)を参照してください。
