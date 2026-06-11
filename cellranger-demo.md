# Cell Ranger Test Run Tutorial

## 1. Download Cell Ranger & Reference

**Environment:** Shared Server

**Purpose:**
Download Cell Ranger and the human reference dataset.

Open the Cell Ranger download page:

[Downloads](https://www.10xgenomics.com/support/jp/software/cell-ranger/downloads?utm_source=chatgpt.com#download-links)

Cell Ranger 10.0.0
1. Under **"tar.gz compression"**, Select the **wget** tab.
2. Click **Copy** to copy the download command.
3. Paste the command into the terminal and execute it.

Reference
1. Under Human reference (GRCh38) - 2024-A, Select the **wget** tab.
2. Click **Copy** to copy the download command.
3. Paste the command into the terminal and execute it.
---
## 2. Extract Cell Ranger & Reference

```bash
tar -xzvf cellranger-10.0.0.tar.gz
```
```bash
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz
```
---
## 3. Add Cell Ranger to PATH

> This command must be executed again after each new login session.

```bash
export PATH=/DATA/cfam000*/cellranger-10.0.0:$PATH
```
---
## 4. Verify Installation

```bash
cellranger testrun --id=check_install
```
---
## 5. Download Example Dataset

Browse available datasets:

https://www.10xgenomics.com/datasets?sort=publishedAt+DESC

Example: PBMC 5K Dataset ([donor1](https://www.10xgenomics.com/datasets/5k_Human_Donor1_PBMC_3p_gem-x) and [donor2](https://www.10xgenomics.com/datasets/5k_Human_Donor2_PBMC_3p_gem-x))

```bash
wget https://cf.10xgenomics.com/samples/cell-exp/9.0.0/5k_Human_Donor1_PBMC_3p_gem-x_Multiplex/5k_Human_Donor1_PBMC_3p_gem-x_Multiplex_fastqs.tar
```
```bash
wget https://cf.10xgenomics.com/samples/cell-exp/9.0.0/5k_Human_Donor2_PBMC_3p_gem-x_Multiplex/5k_Human_Donor2_PBMC_3p_gem-x_Multiplex_fastqs.tar
```
---
## 6. Extract FASTQ Files

```bash
tar -xvf 5k_Human_Donor1_PBMC_3p_gem-x_Multiplex_fastqs.tar
```
```bash
tar -xvf 5k_Human_Donor2_PBMC_3p_gem-x_Multiplex_fastqs.tar
```
---
## 7. Run Cell Ranger Count

```bash
cellranger count \
    --id=donor1 \
    --transcriptome=refdata-gex-GRCh38-2024-A \
    --fastqs=5k_Human_Donor1_PBMC_3p_gem-x_GEX_fastqs \
    --create-bam=false \
    --localcores=32 \
    --localmem=128
```
```bash
cellranger count \
    --id=donor2 \
    --transcriptome=refdata-gex-GRCh38-2024-A \
    --fastqs=5k_Human_Donor2_PBMC_3p_gem-x_GEX_fastqs \
    --create-bam=false \
    --localcores=32 \
    --localmem=128
```
### Parameters

| Parameter            | Description                         |
| -------------------- | ----------------------------------- |
| `--id`               | Output directory name               |
| `--transcriptome`    | Reference transcriptome             |
| `--fastqs`           | FASTQ directory                     |
| `--create-bam=false` | Skip BAM generation to save storage |
| `--localcores=32`    | Number of CPU cores                 |
| `--localmem=128`     | Memory limit (GB)                   |

---

## 8. Download Results to Local Computer

Open a new PowerShell or Terminal window and run:

```bash
scp -r cfam000*@133.41.125.54:/DATA/cfam000*/donor1/outs/ donor1
```
```bash
scp -r cfam000*@133.41.125.54:/DATA/cfam000*/donor2/outs/ donor2
```

---

## 9. Review QC Metrics

Open:

```text
donor1/web_summary.html
```
```text
donor2/web_summary.html
```

in your web browser.

The report includes:

* Sequencing quality metrics
* Mapping statistics
* Cell calling results
* UMI and gene count distributions
* Cell Ranger summary plots
