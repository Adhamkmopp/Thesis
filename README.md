# Bionformatics MSc, Thesis: Viral Presence In The 1000 Genomes Database


Steps:
- [x] Establish a fast bash script for automatic extraction and merging of unmapped reads
- [x] Filter "low complexity reads"
- [x] Create a pure viral database for BLAST/bwa, alongside a mixed database for detecting false positives
- [x] Run blastn query on one sample against the viral database
- [x] Run bwa mem query on one sample against the viral database
- [x] Obtain viral genome sizes and concatenate with blast result

## 1. Command line script in bash:

This script(s) below describe the main process that automates the entire pipeline. Elements will be added after they are executed seperately and validated. Initially, the script serves to extract samples of the 1000 genomes mapped and unmapped bam files by populations, and extracts/merges unmapped and chimeric pair reads. Fastq/a files are then created from the merged bam files in the final step. The merged bam file contain reads sorted by name, while the samtools fastq converter appends a /1 or /2 according to the read flag (forward/backward), and outputs two files split on reads /1 or /2. Finally, blast and bwa are run on the locally built databases (as outlined below) and the results are stored in a seperate folder.


#### Downloading The Files

The script below downloads 20 bam files of the British population, extracts unmapped reads, performs filtering of low complexity reads and reads where 50% of bases have a quality score (Phred) of 3 and below and finally converts the output to FASTA for blasting. the files are stored in a seperate folder by the following order: /samples/population/sample_ID/sample.fasta.
```bash
#!/bin/bash

parent="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/"
brits=$(curl -l $parent | grep -P "HG001[1-5][0-9]$") # obtains and stores a list of 20 or so british individuals

echo $brits


for sample in $brits; do
if [ ! -e "/isdata/common/wbf326/samples/GBR/$sample/$sample.fasta" ];then
	mkdir /isdata/common/wbf326/samples/GBR/$sample # sets up the structure of the directory such that each sample is placed alone by itself
	cd /isdata/common/wbf326/samples/GBR/$sample 
	page="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/$sample/alignment/" # sets up the ftp directory path by inserting the sample
	mapped=$(curl -l $page | grep -P ".+\.mapped.+\.bam$"); # sets up the ftp path to the mapped bam file and unmapped right below
	unmapped=$(curl -l $page | grep -P ".+\.unmapped.+\.bam$");
	echo "Mapped file is: $mapped"
	echo "Unmapped file is: $unmapped"

	echo "Downloading mapped and unmapped bam files.."
	curl $page$mapped --output $sample.mtemp.bam # downloads the mapped bam file from the ftp list
	curl $page$unmapped --output $sample.utemp.bam	 # downloads the unmapped bam file
	echo "Collecting mapped read info.."
	samtools view -@ 64 -c -F 0x4 $sample.mtemp.bam > $sample.stats

	echo "extracting mapped singletons.."
	samtools view -@ 64 -h -F 0x4 -f 0x8 $sample.mtemp.bam -o $sample.mchi.txt # extracts definetly mapped with unmapped pair reads into .mchi

	echo "extracting unmapped singletons.."
	samtools view -@ 64 -h -f 0x4 $sample.mtemp.bam -o $sample.uchi.bam # extracts unmapped reads into .uchi

	echo "merging unmapped singletons with unmapped pairs.."
	samtools merge -@ 64 -n $sample.merged.bam $sample.uchi.bam $sample.utemp.bam # merges ALL mappped and unmapped reads

	echo "creating fastq files from unmapped reads.."
	samtools fastq -@ 64 -N -t -1 $sample.1.fastq -2 $sample.2.fastq $sample.merged.bam # places forward/backward reads into respective fastq files
	
	
	# filters reads by quality and repeats as set in the paper and hopefully outputs a report in the correct place
	sga preprocess -p 0 --quality-filter=50 --no-primer-check --dust-threshold=2.5 --out=$sample.dusted.fastq $sample.1.fastq $sample.2.fastq

	echo "creating fasta from dusted fastq.."
	cat $sample.dusted.fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > $sample.fasta # transforms fastq into fasta using sed and trim

fi

done

```

#### Running blastn On All Samples and Getting Candidate Reads

The purpose of this step is to run blast on the sample.fasta files on the locally built viral database (see: Creation Of A Local Viral Database), filtering out all phage results and creating a set of candidate reads from the filtered blast results. The key steps here are firstly removing duplicate reads from the sample.fasta file (there can be anywhere between 0 and 250 but no more from what I can tell). This is necessary as it interferes with faidx indexing, without which subsetting the sample.fasta file with an ID list would take much too long. The second step (alluded to already) is subsetting the fasta file with an ID list created from the viral blast results, thereby creating the sample.candidateReads.fasta for the second round of comprehensive blast under /viralblast/samples/population/sample_ID/sample.candidateReads.fasta 

Another key point deserve seperate mention; the fasta file is split into 32 seperate parts using fasplit, and 32 blasts are run in parallel to speed up the process significantly.
```bash
#!/bin/bash

parent="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/"
brits=$(curl -l $parent | grep -P "HG001[1-5][0-9]$") # obtains and stores a list of 20 or so british individuals

echo $brits

for sample in $brits; do
	if [[ -d "/isdata/common/wbf326/samples/GBR/$sample" && ! -e "/isdata/common/wbf326/viralblast/samples/GBR/$sample/$sample.candidateReads.fasta" ]]; then
	mkdir /isdata/common/wbf326/viralblast/samples/GBR/$sample
	cd /isdata/common/wbf326/samples/GBR/$sample
	# removing duplicates from fasta
	cat $sample.fasta | /isdata/common/wbf326/./seqkit rmdup -o $sample.clean.fasta 
	# split the candidate read fasta files into 64 seperate files, create a list of all the files and feed it into GNU parallel to run the blast command on
	# the seperate files per available core
	/isdata/common/wbf326/./faSplit sequence $sample.clean.fasta 32 $sample.split
	ls *.fa > falist
	

	echo "blasting...$sample"
	
	# outputs blast results on the total viral database with a custom column format specified elsewhere
	cat falist | parallel --no-notice 'FASTA={};blastn -evalue 1e-10 -query $FASTA -db "/isdata/common/wbf326/viralblast/viraldb.fasta" -outfmt "6 qseqid sseqid evalue bitscore sgi sacc slen qstart qend sstart send stitle" -out $FASTA.blasted'
	cat *.blasted > $sample.blast
	rm *.blasted
	mv $sample.blast /isdata/common/wbf326/viralblast/samples/GBR/$sample/$sample.blast
	
	echo "cleaning up.."
	# append each blast result into one big file
	# cleanup!except for the merged bam
	rm *fa
	rm falist

	echo "collecting viral reads counts.."
	cd /isdata/common/wbf326/viralblast/samples/GBR/$sample
	# collects viral hit counts for each sample split by virus hit
	cut -f 6,12 $sample.blast |sort -r |  uniq -c | sort -n | sed -r 's/([0-9]) /\1\t/' >$sample.vHits
	# collects lengths of genomes found
	#cut -f 2 $sample.vHits | parallel /isdata/common/wbf326/./pcregrep {} /isdata/common/wbf326/vLengths | cut -f 2 > $sample.gSizes
	#paste $sample.gSizes $sample.vHits > $sample.complete
	# removing Phage results (for reads that ONLY matched phages) and creating backup
	echo "removing phage results.."
	sed -i.bak '/phage/d' $sample.blast

	# create ID list of candidate viral reads with regex for extracting full fasta sequence
	echo "collecting candidate reads.."

	cut -f 1 $sample.blast | sort | uniq > IDs

	xargs faidx -f /isdata/common/wbf326/samples/GBR/$sample/$sample.clean.fasta < IDs >> /isdata/common/wbf326/viralblast/samples/GBR/$sample/$sample.candidateReads.fasta
 
	rm IDs
fi 
done



```

#### Running blastn On Candidate Reads Against Comprehensive Database

Candidate reads are split once again and a number of blasts are run in parallel, resulting in a final sample.mixedblast output under mixed/samples/population/sample_ID/sample.mixedblast. Between this step and the previous one however, there is a manual check on the sizes of candidate reads. If a candidate read fasta file is greater than 2MBs, the entire folder is deleted. The script initially checks for the existence of directories within /viralblast and runs blast accordingly. Sizes over 2MBs for a sample take longer than 3 hours to complete, and have often caused the server to run out of memory entirely. This decision was made as a final compromise between blast time, and the memory usage of parallel blast runs.
```bash
#!/bin/bash

parent="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/"
brits=$(curl -l $parent | grep -P "HG001[1-5][0-9]$") # obtains and stores a list of 20 or so british individuals

echo $brits

for sample in $brits; do
	if [[ -d "/isdata/common/wbf326/viralblast/samples/GBR/$sample" &&  ! -e "/isdata/common/wbf326/mixed/samples/GBR/$sample/$sample.mixedblast" ]]; then

	mkdir /isdata/common/wbf326/mixed/samples/GBR/$sample
	cd /isdata/common/wbf326/viralblast/samples/GBR/$sample
	rm *.fa
	rm falist

	/isdata/common/wbf326/./faSplit sequence $sample.candidateReads.fasta 20 $sample.split
	ls *.fa > falist

	echo "blasting candidate reads..$sample"
	# blasting each of the 64 split files seperately and cleanup!
	cat falist | parallel --no-notice 'FASTA={};blastn -evalue 1e-20 -query $FASTA -db "/isdata/common/wbf326/mixed/mixeddb.fasta" -outfmt "6 qseqid sseqid evalue bitscore sgi sacc slen qstart qend sstart send stitle" -out $FASTA.mixed'
	cat *.mixed > $sample.mixedblast
	rm *.mixed
	
	mv $sample.mixedblast /isdata/common/wbf326/mixed/samples/GBR/$sample/$sample.mixedblast

fi
done 


```

#### Post-Processing

This final step concatenate blast results, generates hit counts for quick inspection and collect stats (number of mapped reads to the human genome) in seperate files under their respective home directories. There is nothing complex happening here aside from the R script used to filter out false positives (See: False Positive Tests). Just basic string/column manipulation. 

```bash
#!/bin/bash

parent="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/"
brits=$(curl -l $parent | grep -P "HG001[1-5][0-9]$") # obtains and stores a list of 20 or so british individuals

echo $brits


#ls -d /isdata/common/wbf326/mixed/samples/GBR/*/ > GBRlist # generate list of arguments for the R script in the form of directories


cd /isdata/common/wbf326/scripts
for sample in $brits; do
if [ -d "/isdata/common/wbf326/mixed/samples/GBR/$sample" ]; then
	echo "Running RScript on $sample.."
	Rscript --vanilla analysis.R  "/isdata/common/wbf326/mixed/samples/GBR/$sample/"
	awk -F'\t' -vOFS='\t' -v sample="$sample" '{ $1 = sample"\t" $1 }1' < /isdata/common/wbf326/samples/GBR/$sample/$sample.stats >> /isdata/common/wbf326/mixed/samples/GBR/all.GBR.stats

fi
done

cd /isdata/common/wbf326/mixed/samples/GBR

echo "gathering mixed viral hit count.."
find .  -name '*.mixedblast' -exec cat {} + > all.GBR.mixed
cut -f 6,12 all.GBR.mixed  |sort -r |  uniq -c | sort -n | sed -r 's/([0-9]) /\1\t/' > all.GBR.mHits

echo "gathering mixed filtered viral hit count.."
find .  -name '*.filtered' -exec cat {} + > all.GBR.mixedf
cut -f 6,12 all.GBR.mixedf  |sort -r |  uniq -c | sort -n | sed -r 's/([0-9]) /\1\t/' > all.GBR.mHitsf

cd /isdata/common/wbf326/viralblast/samples/GBR

echo "gathering blast viral hit count.."
find .  -name '*.blast' -exec cat {} + > all.GBR.vblast
cut -f 6,12 all.GBR.vblast  |sort -r |  uniq -c | sort -n | sed -r 's/([0-9]) /\1\t/' > all.GBR.vHits

echo "gathering filtered (no phage) blast viral hit count.."
find .  -name '*.blast.bak' -exec cat {} + > all.GBR.vblastf
cut -f 6,12 all.GBR.vblastf  |sort -r |  uniq -c | sort -n | sed -r 's/([0-9]) /\1\t/' > all.GBR.vHitsf

cd /isdata/common/wbf326/mixed/samples/GBR
echo "tarring stats and all mixed filtered blast hits.."
find . -type f -name '*.filtered' | cat > flist
tar cvf GBR.mixedf.tar -T flist
```


## 2. The Sundries


#### Filtering Out Low Complexity Reads (dust-threshold=2.5/quality filter=50) 

SGA filters low complexity reads by counting the number of triplet runs in a 64 base window (given that all possible triplets of four bases is 64). So a single occurance of a triplet contributes 0 to an added score which increments with each new occurance. The formula then adds all runs, divided by the window, and multiplies by 10. The signficance of the 2.5 thredhold is still not clear. The quality filter also removes reads where 50% of the bases have a phred score of 3 or below.

#### Creation Of A Local Viral Database:
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

The database was then transformed into a BLAST-ready local database using blast+ with the command:

```bash
makeblastdb -in "viraldb.fasta" -dbtype nucl -parse_seqids -out "viraldb.fasta"
```

#### Creation of Mixed/Comprehensive Database for Blasting Candidate Reads:
The comprehensive database consisted of representaitve genomes for Fungi, Bacteria and Archae, alongwith UCSC human, chimp, mouse and chicken genomes, mixed with the refseq viral database.

Potential representative genomes from the assembly database:

Fungi(14/2)
Bacteria(1,543/1636)
Archaea(140/389)

The second number is the one reported in the paper. The search criteria for assembled genomes as listed below was used to extract the representative sequences:

```bash
("Bacteria"[Organism] OR "Fungi"[Organism] OR "Archaea" [Organism]) AND (latest[filter] AND "complete genome"[filter] AND "representative genome"[filter]) 
```

Later, the UCSC genomes were downloaded individually, converted from 2bit format to fasta and screen for non-unique read names before concatenating all the above as one large 17GB fasta file. To illustrate the first step, several entries are listed in different genomes as "chr(1 to 21) giving non-unique names for entires were concatenating (making it impossible to create a local database). To get around that issue, entries were modified with 'sed -i 's/old_entry/new_entry/g' where chrX for example, was renamed to chrX_chimp for all five genomes.



A BWA index was constructed from the same ncbi RefSeq fna/fasta file using the command:

```bash
bwa index -p viral -a is /home/adhamkmopp/viral.fna
```

#### Blast Output Format

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
#### False Positives Tests


This R script filters false positives on the final blast results. It removes hits where there are hits to organisms other than viruses with a bit score higher than or equal to that hit. So in principle, a read with multiple hits, several to viral genomes and to other organisms, may have some viral hits removed if it has a lower bit score than ANY of the hits to non-viral organisms. The input is the directory where the sample.mixedblast file is, the output is placed in the same directory as sample.mixed.filtered.



```r
#!/usr/bin/env Rscript
library(data.table)
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
files <- list.files(pattern = "*.mixedblast") # find all files with mixedblast extension
out <-sub(pattern = "(.*)\\..*$", replacement = "\\1", files) # get the sample name only


vlist <- scan("/isdata/common/wbf326/scripts/vNames", what="", sep="\n") # create virus list of names
dt<- fread(files, header = F, sep = "\t")
colnames(dt) <- c("query", "sequence", "evalue", "bitscore", "sgi", "sacc", "sLength", "qstart", "qend", "sstart",
                  "send", "stitle")


ss <- split(dt,dt$query) # split the dataframe into a list of dataframes each with a single query ID and all of its hits
z <- data.frame()
for(i in 1:length(ss)){
  boo=NULL
  x<-as.data.frame(ss[[i]]) # extract read i
  idx <- which (is.element(x$stitle, vlist)) # find where read i hits match viruses and subset in t below
  t <- x[idx,]
  for (i in nrow(t)){
    boo <- c(boo, all(t[i,4]> x[-idx,4])) # add boolean checks where each viral hit scored higher than all other non-viral hits
  }
  z<-rbind(z, t[boo,]) # concatenate result to empty dataframe and increase its size and the loop continues and finds more viral hits with scores higher than non-viral hits
}
 
z<-z[complete.cases(z),] # remove NA values
z<-cbind(z,out)
  

write.table(z, paste(args,out,".filtered", sep=""), sep = "\t", col.names = F, row.names = F, quote = F) # output in the same directory with  sample.filtered extension

```
