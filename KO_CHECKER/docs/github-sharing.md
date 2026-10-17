# GitHubへの共有手順

この手順はソースコードの共有用です。GitHubへ置くだけで解析サーバーが起動するわけではありません。利用者はソースを取得し、自分のPCやサーバーでDockerを起動します。Python、HISAT2、SRA Toolkitを実行するため、GitHub Pagesだけではこのアプリを動かせません。

## 共有用ZIP

プロジェクトのフォルダで次を実行します。

```powershell
py prepare_release.py
```

Linuxでは`python prepare_release.py`を使います。Python 3.9以上の標準ライブラリだけで作成できます。

`dist/rnako-github.zip`には、明示したソース・設定・説明書だけが入ります。ZIP内の`SHA256SUMS.txt`は各ファイルのSHA-256一覧です。ZIP自体とローカルデータはGit追跡から除外します。公開用に`LICENSE`を追加した場合は、次回のZIP作成時に含まれます。

作成時には、作業フォルダ直下の`SHA256SUMS.txt`も同じ内容に更新します。ZIPにはWindows用の`start.bat`・`start.ps1`、中断・再開用の`transfers.py`・`test_resume.py`も含まれます。既存環境へ更新する場合はソースファイルを更新し、実行中の解析がないことを確認して`start.bat`または`docker compose up -d --build`で再ビルドしてください。`references/`・`workdata/`・`results/`は保持してください。

## GitHubの画面から追加する

1. GitHubで新しいリポジトリを作成します。公開範囲を選択してください。
2. `rnako-github.zip`をPC上の別フォルダへ展開します。
3. リポジトリの「Add file」→「Upload files」から、**展開した中身**を追加します。ZIPだけをアップロードするのではなく、`README.md`や`Dockerfile`がリポジトリの直下になるようにします。
4. `.github`フォルダ、`.gitignore`、`.dockerignore`、`.gitattributes`、`.env.example`も含めます。隠しファイルを選択できない場合はGitHub DesktopやGitを使ってください。
5. 追加するファイル一覧を確認し、コミットします。

GitHub Actionsの設定は`.github/workflows/tests.yml`です。ワークフローを追加できる権限・認証方法でアップロードしてください。CIはDockerをビルドして人工データのテストを実行し、公開RNA-seq本体や参照ゲノムは取得しません。

## Gitを使う場合

Gitが利用できる環境で、展開したソースのフォルダから実行します。`YOUR_ACCOUNT`と`YOUR_REPOSITORY`を自分の値に置き換えてください。

```bash
git init
git add .
git status --short
git diff --cached --stat
git commit -m "Add KO_Checker web prototype"
git branch -M main
git remote add origin https://github.com/YOUR_ACCOUNT/YOUR_REPOSITORY.git
git push -u origin main
```

上記は、GitHub側でまだコミットを作っていない新規リポジトリへの例です。既存リポジトリでは、その履歴に合わせて追加してください。

## 公開用に含めないもの

- `inputs/`内の実検体CSVや個人のファイルパス
- `references/`、`workdata/`、`results/`内の参照、取得データ、解析結果、ログ
- `.env`、認証情報、仮想環境、`local_notes/`内の引き継ぎ記録
- FASTQ、SRA、BAM、ゲノム配列、HISAT2索引などの大きなデータ

`.gitignore`は新しく追加するファイルを除外する設定です。すでにGitで追跡しているファイルや過去のコミットから情報を消す機能ではありません。共有用ZIPは、除外パターンではなくファイルの明示リストで作成しています。

## ライセンス

このプロジェクトには、現時点で`LICENSE`を追加していません。公開先の管理者が利用・改変・再配布の条件を決めてから、適切なライセンスファイルを追加してください。コードを公開することと、オープンソースの利用許諾を付けることは別です。

参考: [GitHubのライセンス案内](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository)。依存ツールと参照データには、それぞれの利用条件が適用されます。
