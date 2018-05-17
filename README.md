# Bionformatics MSc, Thesis: Viral Presence In The 1000 Genomes Database


Steps:
- [x] Establish a fast command line script for automatic extraction and merging of unmapped reads
- [ ] Filter "low complexity reads"
- [x] Create a pure viral database for BLAST/bwa, alongside a mixed database for detecting false positives
- [x] Run blastn query on one sample against the viral database
- [x] Run bwa mem query on one sample against the viral database

## 1. Command line script in bash:

This script below will be the main script that automates the entire pipeline. Elements will be added after they are executed seperately and validated. Initially, the script serves to extract samples of the 1000 genomes mapped and unmapped bam files, and extracts/merges unmapped and chimeric pair reads. Fastq/a files are created from the merged bam files in the final step. The merged bam file contains reads sorted by name, while the samtools fastq converter appends a /1 or /2 according to the read flag (forward/backward), and outputs a single interleaved file.
```bash
#!/bin/bash

parent="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/"
all=$(curl -l $parent | grep -P "HG001[012]{1,2}$")

echo $all

for sample in $all; do
	page="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/$sample/alignment/"
	mapped=$(curl -l $page | grep -P ".+\.mapped.+\.bam$");
	unmapped=$(curl -l $page | grep -P ".+\.unmapped.+\.bam$");
	echo $mapped
	echo $unmapped
	curl $page$mapped --output $sample.mtemp.bam
	curl $page$unmapped --output $sample.utemp.bam	
	samtools view -@ 64 -h -F 0x4 -f 0x8 $sample.mtemp.bam -o $sample.mchi.bam 
	samtools view -@ 64 -h -f 0x4 $sample.mtemp.bam -o $sample.uchi.bam
	samtools merge -@ 64 -n $sample.merged.bam $sample.mchi.bam $sample.uchi.bam $sample.utemp.bam
	samtools sort -@ 64 -n $sample.merged.bam -o $sample.sorted.bam
	rm $sample.mchi.bam 
	rm $sample.uchi.bam
	rm $sample.utemp.bam
	rm $sample.mtemp.bam
	samtools fastq -@ 64 -N -t $sample.sorted.bam > $sample.fastq
	samtools fasta -@ 64 -N -t $sample.sorted.bam > $sample.fasta
	
done


```
## 2. Filter out low complexity reads and low quality reads (dust-threshold=2.5/quality filter=50) 

Nothing yet

## 3. Creation of a viral database for local blast:
The viral refseq database was downloaded from here: https://www.ncbi.nlm.nih.gov/genome/viruses/ after filtering for phages:

```bash
Viruses[Organism] AND srcdb_refseq[PROP] NOT unclassified dsDNA phages[Organism] NOT unclassified virophages[organism] NOT "phg"[Division] NOT wgs[PROP] NOT cellular organisms[ORGN] NOT AC_000001:AC_999999[PACC] 
```
Potential representative genomes from the assembly database:

Fungi(14/2)
Bacteria(1,543/1636)
Archaea(140/389)


```bash
("Bacteria"[Organism] OR "Fungi"[Organism] OR "Archaea" [Organism]) AND (latest[filter] AND "complete genome"[filter] AND "representative genome"[filter]) 
```

It was then transformed into a BLAST-ready local database using blast+ with the command:

```bash
makeblastdb -in "$BLASTDB/viral.fna" -dbtype nucl -parse_seqids -out "$BLASTDB\viral.fna"
```
Where $BLASTDB is the path to the fna/fasta file containing the viral genomes (~7,000 in total)

A BWA index was constructed from the same ncbi RefSeq fna/fasta file using the command:

```bash
bwa index -p viral -a is /home/adhamkmopp/viral.fna
```

## 4. Blast/bwa on test sample HG0100.bam

Before doing a blast seach, a taxonomic ID dump was obtained from (ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz) to allow the inclusion of scientific names in the results. The custom output was as follows:

| qseqid | sseqid | evalue | bitscore | sgi | sacc | staxids | sscinames | scomnames | stitle |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | 
| Query ID | Seq ID | E-Value | Bit score | Subject GI | Subject Accession | Subject Tax ID | Subject Scientific Name | Subject Common name | Subject Title |


```bash
blastn -num_threads 64 -query HG0100.fasta -db "viraldb.fasta" -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle'  > HG0100.blast
```


```bash
bwa mem -p viral HG0100.fastq > viralbwa
```
