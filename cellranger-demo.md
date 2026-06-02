# Cell Ranger Test Run Tutorial

## 1. Download Cell Ranger & Reference

Open the Cell Ranger download page:

[Cell Ranger Downloads](https://www.10xgenomics.com/support/jp/software/cell-ranger/downloads?utm_source=chatgpt.com#download-links)

Cell Ranger 10.0.0
1. Under **"tar.gz compression"**, Select the **wget** tab.
2. Click **Copy** to copy the download command.
3. Paste the command into the terminal and execute it.

References
1. Under Human reference (GRCh38) - 2024-A, Select the **wget** tab.
2. Click **Copy** to copy the download command.
3. Paste the command into the terminal and execute it.

## 2. Extract Cell Ranger & Reference

```bash
tar -xzvf cellranger-10.0.0.tar.gz
```
```bash
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz
```

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

https://www.10xgenomics.com/support/software/cell-ranger/downloads#download-links

Example: PBMC 10K dataset

```bash
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/4.0.0/SC3_v3_NextGem_SI_PBMC_10K/SC3_v3_NextGem_SI_PBMC_10K_fastqs.tar
```

---

## 6. Extract FASTQ Files

```bash
tar -xvf SC3_v3_NextGem_SI_PBMC_10K_fastqs.tar
```

---

## 7. Reduce Dataset Size for Testing

> Remove lanes 3 and 4 to shorten runtime and reduce resource usage.

```bash
rm -f SC3_v3_NextGem_SI_PBMC_10K_fastqs*L003* \
      SC3_v3_NextGem_SI_PBMC_10K_fastqs*L004*
```

---

## 8. Run Cell Ranger Count

```bash
cellranger count \
    --id=pbmc \
    --transcriptome=refdata-gex-GRCh38-2024-A \
    --fastqs=SC3_v3_NextGem_SI_PBMC_10K_fastqs \
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

## 9. Download Results to Local Computer

Open a new PowerShell or Terminal window and run:

```bash
scp -r cfam0001@133.41.125.54:/DATA/cfam0001/pbmc/outs/ .
```

---

## 10. Review QC Metrics

Open:

```text
outs/web_summary.html
```

in your web browser.

The report includes:

* Sequencing quality metrics
* Mapping statistics
* Cell calling results
* UMI and gene count distributions
* Cell Ranger summary plots
