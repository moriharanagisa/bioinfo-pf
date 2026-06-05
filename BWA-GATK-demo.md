# GATK Variant Calling Pipeline

## Extract FASTQ files

```bash
# Login
ssh gene@10.34.21.182

# Extract tar file
cd /mnt/add2/
tar -xvf /mnt/add2/X207SC26044118-Z01-F001.tar

# Set project name
export PROJECT="260501-miyoshi"

# Create project directory
cd /DATA/centos/GATK
sudo mkdir ${PROJECT}
cd ${PROJECT}

# Copy FASTQ files
sudo cp /mnt/add2/X207SC26044118-Z01-F001/01.RawData/*/*.fq.gz ./
```

---

## Trimming

```bash
# Enter Trim Galore container
docker start trimgalore-0.6.10
docker attach trimgalore-0.6.10

export PROJECT="260501-miyoshi"
cd /home/GATK/${PROJECT}

# Trim paired-end reads
for r1 in *_1.fq.gz; do
    r2="${r1/_1.fq.gz/_2.fq.gz}"
    trim_galore --paired "$r1" "$r2"
done
```

```text
# Detach container
Ctrl + P, Ctrl + Q
```

---

## Alignment

```bash
docker attach centos-DATA
export PROJECT="260501-miyoshi"
cd /home/GATK/GATK
```

```bash
process_sample() {
    fq1="$1"
    fq2="${fq1/_1_val_1.fq.gz/_2_val_2.fq.gz}"
    sample=$(basename "$fq1" _1.fq.gz)

    ../../bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem \
        -t 96 \
        -R "@RG\tID:$sample\tSM:$sample\tPL:illumina\tLB:$sample" \
        Homo_sapiens_assembly38.fasta \
        "$fq1" "$fq2" | \
        samtools sort -@ 4 -o "${sample}.sort.bam"
}

export -f process_sample

find ../${PROJECT}/ -name "*_1_val_1.fq.gz" | \
parallel -j 4 process_sample
```

```text
# Detach container
Ctrl + P, Ctrl + Q
```

---

## Remove duplicates

```bash
cd /DATA/centos/GATK/GATK

for filepath in *.sort.bam; do

    filename=$(basename "$filepath" .sort.bam)

    docker run --rm \
        -v "$(pwd)":/data \
        broadinstitute/picard:3.1.0 \
        java -jar /usr/picard/picard.jar MarkDuplicates \
        --INPUT /data/"$filename".sort.bam \
        --OUTPUT /data/"$filename".sort.dedup.bam \
        --METRICS_FILE /data/"$filename".marked_dup_metrics.txt \
        --VALIDATION_STRINGENCY LENIENT \
        --REMOVE_DUPLICATES true

done
```

---

## Index BAM files

```bash
docker attach centos-DATA
```

```bash
for bamfile in *.sort.dedup.bam; do
    samtools index "$bamfile"
done
```

---

## Base quality recalibration

```bash
parallel -j 4 '
sample=$(basename {} .sort.dedup.bam)

./gatk-4.3.0.0/gatk BaseRecalibrator \
    --input {} \
    --reference Homo_sapiens_assembly38.fasta \
    --output ${sample}.dedup.recaltab.txt \
    --known-sites Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
' ::: *.sort.dedup.bam
```

```bash
parallel -j 4 '
sample=$(basename {} .sort.dedup.bam)

./gatk-4.3.0.0/gatk ApplyBQSR \
    --input {} \
    --reference Homo_sapiens_assembly38.fasta \
    --bqsr-recal-file ${sample}.dedup.recaltab.txt \
    --output ${sample}.recal.bam
' ::: *.sort.dedup.bam
```

---

## Variant calling

```bash
ls *.recal.bam | parallel -j 4 '
sample=$(basename {} .recal.bam)

./gatk-4.3.0.0/gatk \
    --java-options "-Xmx4g" \
    HaplotypeCaller \
    -I {} \
    -O ${sample}.g.vcf.gz \
    -R Homo_sapiens_assembly38.fasta \
    -L 50genes.bed \
    --ERC GVCF \
    --dbsnp Homo_sapiens_assembly38.dbsnp138.vcf
'
```

---

## Genotyping

```bash
ls *.g.vcf.gz | parallel -j 4 '
sample=$(basename {} .g.vcf.gz)

./gatk-4.3.0.0/gatk \
    --java-options "-Xmx4g" \
    GenotypeGVCFs \
    -V {} \
    -O ${sample}.both.vcf.gz \
    -R Homo_sapiens_assembly38.fasta \
    -L 50genes.bed
'
```

---

## Separate SNPs and INDELs

```bash
ls *.both.vcf.gz | parallel -j 4 '
sample=$(basename {} .both.vcf.gz)

./gatk-4.3.0.0/gatk \
    SelectVariants \
    -R Homo_sapiens_assembly38.fasta \
    -V {} \
    --select-type-to-include SNP \
    -O ${sample}.snv.raw.vcf.gz

./gatk-4.3.0.0/gatk \
    SelectVariants \
    -R Homo_sapiens_assembly38.fasta \
    -V {} \
    --select-type-to-include INDEL \
    -O ${sample}.indel.raw.vcf.gz
'
```

---

## Variant filtration

### SNP

```bash
ls *.snv.raw.vcf.gz | parallel -j 4 '
sample=$(basename {} .snv.raw.vcf.gz)

./gatk-4.3.0.0/gatk VariantFiltration \
    -R Homo_sapiens_assembly38.fasta \
    -V {} \
    -O ${sample}.snv.pass.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
'
```

### INDEL

```bash
ls *.indel.raw.vcf.gz | parallel -j 4 '
sample=$(basename {} .indel.raw.vcf.gz)

./gatk-4.3.0.0/gatk VariantFiltration \
    -R Homo_sapiens_assembly38.fasta \
    -V {} \
    -O ${sample}.indel.pass.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"
'
```

---

## Merge SNP and INDEL VCFs

```bash
# Detach container
Ctrl + P, Ctrl + Q
```

```bash
for filepath in *.snv.pass.vcf.gz; do

    filename=$(basename "$filepath" .snv.pass.vcf.gz)

    docker run --rm \
        -v "$(pwd)":/data \
        broadinstitute/picard:3.1.0 \
        java -jar /usr/picard/picard.jar MergeVcfs \
        --INPUT /data/"$filename".snv.pass.vcf.gz \
        --INPUT /data/"$filename".indel.pass.vcf.gz \
        --OUTPUT /data/"$filename".all.pass.vcf.gz

done
```

---

## ANNOVAR annotation

```bash
docker attach centos-DATA
```

```bash
ANNOVAR_DB="./annovar-hg38"
INFO_DB="tommo-60kjpn-20240904-af_snvindel.INFO.genericdb"
MAF_DB="tommo-60kjpn-20240904-af_snvindel.MAF.genericdb"
```

```bash
process_file() {

    local file="$1"
    local sample

    sample=$(basename "$file" .all.pass.vcf.gz)

    ./annovar/convert2annovar.pl \
        --format vcf4 \
        --outfile "${sample}.all.pass.avinput" \
        --filter PASS \
        <(gunzip -c "$file")

    ./annovar/table_annovar.pl \
        "${sample}.all.pass.avinput" \
        "$ANNOVAR_DB" \
        --buildver hg38 \
        --outfile "$sample" \
        --remove \
        --protocol cytoBand,refGene,avsnp150,generic,generic \
        --operation r,g,f,f,f \
        --genericdbfile "${INFO_DB},${MAF_DB}" \
        --argument ",--hgvs --exonicsplicing,,," \
        --nastring . \
        --polish \
        --otherinfo
}

export -f process_file
export ANNOVAR_DB INFO_DB MAF_DB

parallel -j 4 process_file ::: *.all.pass.vcf.gz
```

---

## ClinVar annotation

```bash
for vcf in *.all.pass.vcf.gz; do

    sample=${vcf%.all.pass.vcf.gz}

    SnpSift annotate \
        clinvar_20250615.vcf.gz \
        "$vcf" \
        > "${sample}_clinvar.vcf"

done
```

---

## Coverage analysis

```bash
# Detach container
Ctrl + P, Ctrl + Q
```

```bash
for bamfile in *_1_val.sort.dedup.bam; do

    filename=$(basename "$bamfile" _1_val.sort.dedup.bam)

    docker run --rm \
        -v "$(pwd)":/opt/mount \
        quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2 \
        mosdepth \
        -n \
        --by /opt/mount/50genes.bed \
        --fast-mode \
        -t 4 \
        /opt/mount/"$filename" \
        /opt/mount/"$filename"_1_val.sort.dedup.bam

done
```

---

## MultiQC report

```bash
multiqc .
```
