# Cell Ranger Test Run Tutorial

## 1. Download Cell Ranger

```bash
wget -O cellranger-10.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.0.0.tar.gz?Expires=1779986540&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=AB1uTJWPCtPHTnkX7Ic~JdcQpL7alzgGjoUvZ42k1RkT9LzI5YCXEbFD1zVX3H6AzFUBBRiq1jNBjgVTXHBaI-KAJrbiHplNr3G5SE0hSyHASbusgADF59kDv1lgjN9vEDYb3CyreEPV~tXljb7JJ-7q5WsnIYxHXTnxyMntav5s52Xy-qRZcNLx5TdufGEcQlNr1Yr5UyvgcMbFtWI7nzvOVClK86e04nP~OQvbOuWs9u1opm2e365XxUSR03o~IVdt9AzrzoaSgPyxgn4p3h2Cwd-MoNBRjay82WCGHz69Cuj2PPbenad9ZzuayQ2X9Mbsk9Y0nbQM94HYbunnnw__"
```

## 2. Extract Cell Ranger

```bash
tar -xzvf cellranger-10.0.0.tar.gz
```

---

## 3. Download Human Reference Genome

```bash
wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
```

## 4. Extract Reference Files

```bash
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz
```

---

## 5. Add Cell Ranger to PATH

> This command must be executed again after each new login session.

```bash
export PATH=/DATA/cfam0001/cellranger-10.0.0:$PATH
```

---

## 6. Verify Installation

```bash
cellranger testrun --id=check_install
```

---

## 7. Download Example Dataset

Browse available datasets:

https://www.10xgenomics.com/support/software/cell-ranger/downloads#download-links

Example: PBMC 10K dataset

```bash
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/4.0.0/SC3_v3_NextGem_SI_PBMC_10K/SC3_v3_NextGem_SI_PBMC_10K_fastqs.tar
```

---

## 8. Extract FASTQ Files

```bash
tar -xvf SC3_v3_NextGem_SI_PBMC_10K_fastqs.tar
```

---

## 9. Reduce Dataset Size for Testing

> Remove lanes 3 and 4 to shorten runtime and reduce resource usage.

```bash
rm -f SC3_v3_NextGem_SI_PBMC_10K_fastqs*L003* \
      SC3_v3_NextGem_SI_PBMC_10K_fastqs*L004*
```

---

## 10. Run Cell Ranger Count

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

## 11. Download Results to Local Computer

Open a new PowerShell or Terminal window and run:

```bash
scp -r cfam0001@133.41.125.54:/DATA/cfam0001/pbmc/outs/ .
```

---

## 12. Review QC Metrics

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
