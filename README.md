# Bionformatics MSc, Thesis: Viral Presence In The 1000 Genomes Database


Steps:
- [x] Establish a fast command line script for automatic extraction and merging of unmapped reads
- [ ] Filter "low complexity reads"
- [x] Create a pure viral database for BLAST/bwa, alongside a mixed database for detecting false positives
- [x] Run blastn query on one sample against the viral database
- [x] Run bwa mem query on one sample against the viral database

1. Command line script in bash:

This script extracts a sample of the 1000 genomes mapped and unmapped bam files, and extracts/merges unmapped and chimeric pair reads. Fastq files are created from the merged bam files in the final step. The merged bam file contains reads sorted by name, while the samtools fastq converter appends a /1 or /2 according to the read flag (forward/backward)
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
	samtools view -@ 2048 -h -F 0x4 -f 0x8 $sample.mtemp.bam -o $sample.mchi.bam 
	samtools view -@ 2048 -h -f 0x4 $sample.mtemp.bam -o $sample.uchi.bam
	samtools merge -@ 2048 -n $sample.merged.bam $sample.mchi.bam $sample.uchi.bam $sample.utemp.bam
	rm $sample.mchi.bam 
	rm $sample.uchi.bam
	rm $sample.utemp.bam
	rm $sample.mtemp.bam
	samtools fastq -@ 2048 -N -t $sample.merged.bam > $sample.fastq
done

```
2. Filter out low complexity reads and low quality reads (dust-threshold=2.5/quality filter=50) 

3. Creation of a viral database for local blast:
The viral refseq database was downloaded from here: https://www.ncbi.nlm.nih.gov/genome/viruses/, and transformed into a BLAST-ready database using blast+ with the command:
```bash
makeblastdb -in "$BLASTDB/viral.fna" -dbtype nucl -parse_seqids -out "$BLASTDB\viral.fna"
```

But firstly, the database had to be filtered for phages. The search below through ncbi yielded the database minus phages, and unclassified phages:

```bash
Viruses[Organism] AND srcdb_refseq[PROP] NOT unclassified dsDNA phages[Organism] NOT unclassified virophages[organism] NOT "phg"[Division] NOT wgs[PROP] NOT cellular organisms[ORGN] NOT AC_000001:AC_999999[PACC] 
```


bwa index construction:
```bash
bwa index -p viral -a is /home/adhamkmopp/viral.fna
```

4. Blast/bwa on test sample HG0100.bam

```bash
blastn -query HG0100.fasta -db "$BLASTDB/viral.fna"
bwa mem -p viral HG0100.fastq > viralbwa
```
