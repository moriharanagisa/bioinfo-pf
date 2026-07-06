## 1. Activate the virtual environment

```bash
source /DATA/shared/alphafold3/apps/alphafold3/.venv/bin/activate
```

## 2. Run AlphaFold3

Replace `***` with the correct file or directory path before running the command.

```bash
python3.12 \
  /DATA/shared/alphafold3/apps/alphafold3/run_alphafold.py \
  --json_path /DATA/cfam0***/*** \
  --model_dir /DATA/cfam0***/*** \
  --db_dir /DATA/shared/alphafold3/db \
  --output_dir /DATA/cfam0***/*** \
  --jackhmmer_binary_path /DATA/shared/alphafold3/apps/hmmer/3.4/bin/jackhmmer \
  --hmmsearch_binary_path /DATA/shared/alphafold3/apps/hmmer/3.4/bin/hmmsearch \
  --hmmbuild_binary_path /DATA/shared/alphafold3/apps/hmmer/3.4/bin/hmmbuild \
  --hmmalign_binary_path /DATA/shared/alphafold3/apps/hmmer/3.4/bin/hmmalign \
  --nhmmer_binary_path /DATA/shared/alphafold3/apps/hmmer/3.4/bin/nhmmer
```

## 3. Deactivate the virtual environment

```bash
deactivate
```
