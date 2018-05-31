# Bionformatics MSc, Thesis: Viral Presence In The 1000 Genomes Database


Steps:
- [x] Establish a fast bash script for automatic extraction and merging of unmapped reads
- [x] Filter "low complexity reads"
- [x] Create a pure viral database for BLAST/bwa, alongside a mixed database for detecting false positives
- [x] Run blastn query on one sample against the viral database
- [x] Run bwa mem query on one sample against the viral database

## 1. Command line script in bash:

This script below will be the main script that automates the entire pipeline. Elements will be added after they are executed seperately and validated. Initially, the script serves to extract samples of the 1000 genomes mapped and unmapped bam files, and extracts/merges unmapped and chimeric pair reads. Fastq/a files are created from the merged bam files in the final step. The merged bam file contains reads sorted by name, while the samtools fastq converter appends a /1 or /2 according to the read flag (forward/backward), and outputs a two files split on reads /1 or /2. Finally, blast and bwa are run on the locally built databases (as outlined below) and the results are stored in a seperate folder.

The blast results are aggregated by counts for the top 20 hits found, and the corresponding reads for each of these hits are extracted and placed into a seperate file.
```bash
#!/bin/bash

parent="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/"
brits=$(curl -l $parent | grep -P "HG001[1-5]{1,3}$") # obtains and stores a list of 20 or so british individuals

echo $brits


for sample in $brits; do
	mkdir $sample # sets up the structure of the directory such that each sample is placed alone by itself
	mkdir viralblast
	mkdir viralbwa
	mkdir mixedblast
	mkdir /home/wbf326/viralblast/$sample
	mkdir /home/wbf326/viralbwa/$sample
	mkdir /home/wbf326/mixedblast/$sample
	cd $sample 
	
	page="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/$sample/alignment/" # sets up the ftp directory path by inserting the sample
	mapped=$(curl -l $page | grep -P ".+\.mapped.+\.bam$"); # sets up the ftp path to the mapped bam file and unmapped right below
	unmapped=$(curl -l $page | grep -P ".+\.unmapped.+\.bam$");
	echo $mapped
	echo $unmapped
	curl $page$mapped --output $sample.mtemp.bam # downloads the mapped bam file from the ftp list
	curl $page$unmapped --output $sample.utemp.bam	 # downloads the unmapped bam file
	samtools view -@ 64 -h -F 0x4 -f 0x8 $sample.mtemp.bam -o $sample.mchi.bam # extracts definetly mapped with unmapped pair reads into .mchi
	samtools view -@ 64 -h -f 0x4 $sample.mtemp.bam -o $sample.uchi.bam # extracts unmapped reads into .uchi
	samtools merge -@ 64 -n $sample.merged.bam $sample.mchi.bam $sample.uchi.bam $sample.utemp.bam # merges ALL mappped and unmapped reads
	samtools fastq -@ 64 -N -t -1 $sample.1.fastq -2 $sample.2.fastq $sample.sorted.bam # places forward/backward reads into respective fastq files

	# filters reads by quality and repeats as set in the paper and hopefully outputs a report in the correct place
	sga preprocess -p 1 --quality-filter=50 --no-primer-check --dust-threshold=2.5 --out=$sample.dusted.fastq $sample.1.fastq $sample.2.fastq > 		$sample.dustreport 
	cat $sample.dusted.fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > $sample.fasta # transforms fastq into fasta using sed and trim

	# outputs blast results on the total viral database with a custom column format specified elsewhere
	blastn -num_threads 64 -query $sample.fasta -db "/home/wbf326/viraldb.fasta" -outfmt '6 qseqid sseqid evalue bitscore sgi sacc qstart qend sstart send stitle' >/home/wbf326/viralblast/$sample/$sample.blast
	bwa mem -t 64 -p /home/wbf326/viraldb $sample.dusted.fastq >/home/wbf326/viralbwa/$sample/$sample.bwa

	# cleanup!except for the merged bam
	rm $sample.mchi.bam 
	rm $sample.uchi.bam
	rm $sample.utemp.bam
	rm $sample.mtemp.bam

	samtoools -@ 64 flagstat $sample.sorted.bam > $sample.flagstat # sanity check on the number of unmapped and mapped reads
	cd /home/wbf326/viralblast/$sample
	cut -f 11 $sample.blast |sort -r |  uniq -c | sort -n | tail -20 | sed -r 's/([0-9]) /\1\t/' >$sample.vHits # collects top 20 hits with counts into a single file
	cut -f 2 $sample.vHits | while read hit; do awk "/$hit/" sorted > "$hit" ; done # loops over hits and gets reads that mapped to that particular viral hit
	cd /home/wbf326

done
```
## 2. Filter out low complexity reads and low quality reads (dust-threshold=2.5/quality filter=50) 

SGA filters low complexity reads by counting the number of triplet runs in a 64 base window (given that all possible triplets are 64). So a single occurance of a triplet contributes 0 to an added score which increments with each new find. The formula then adds all runs, divided by the window, and multiplies by 10. The signficance of the 2.5 thredhold is still not clear. The quality filter removes reads where 50% of the bases have a phred score of 3 or below.

## 3. Creation of a viral database for local blast:
The viral refseq database was downloaded from here: https://www.ncbi.nlm.nih.gov/genome/viruses/ after filtering for phages:

```bash
Viruses[Organism] AND srcdb_refseq[PROP] NOT wgs[PROP] NOT cellular organisms[ORGN] NOT AC_000001:AC_999999[PACC] 
```
The above custom search could easily be modified to exclude certain viruses like phages, if the need arises.

The perl script used to download the sequences based on the custom search above:
```perl
use LWP::Simple;

# Download PubMed records that are indexed in MeSH for both asthma and 
# leukotrienes and were also published in 2009.


$db = 'nuccore';
$query = 'Viruses[Organism]+AND+srcdb_refseq[PROP]+NOT+wgs[PROP]+NOT+cellular+organisms[ORGN]+NOT+AC_000001:AC_999999[PACC]';

#assemble the esearch URL
$base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";

#post the esearch URL
$output = get($url);

#parse WebEnv and QueryKey
$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);

### include this code for ESearch-EFetch
#assemble the efetch URL
$url = $base . "efetch.fcgi?db=$db&query_key=$key&WebEnv=$web";
$url .= "&rettype=fasta";

#post the efetch URL
$data = get($url);
print "$data";
```
Potential representative genomes from the assembly database:

Fungi(14/2)
Bacteria(1,543/1636)
Archaea(140/389)

The second number is the one reported in the paper.

```bash
("Bacteria"[Organism] OR "Fungi"[Organism] OR "Archaea" [Organism]) AND (latest[filter] AND "complete genome"[filter] AND "representative genome"[filter]) 
```

The database was then transformed into a BLAST-ready local database using blast+ with the command:

```bash
makeblastdb -in "viraldb.fasta" -dbtype nucl -parse_seqids -out "viraldb.fasta"
```

A BWA index was constructed from the same ncbi RefSeq fna/fasta file using the command:

```bash
bwa index -p viral -a is /home/adhamkmopp/viral.fna
```

## 4. Blast/bwa on test samples

Before doing a blast seach, a taxonomic ID dump was obtained from (ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz) to allow the inclusion of scientific names in the results. The custom output was as follows:

| qseqid | sseqid | evalue | bitscore | sgi | sacc | qstart | qend | sstart | send | stitle |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Query ID | Seq ID | E-Value | Bit score | Subject GI | Subject Accession | Query Start| Query End | Subject Start | | Subject End | Subject Title |


```bash
blastn -num_threads 64 -query HG0100.fasta -db "viraldb.fasta" -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle'  > HG0100.blast
```


```bash
bwa mem -p viral HG0100.fastq > viralbwa
```
