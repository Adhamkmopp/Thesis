# Bionformatics MSc, Thesis: Viral Presence In The 1000 Genomes Database


Steps:
- [x] Establish a fast bash script for automatic extraction and merging of unmapped reads
- [x] Filter "low complexity reads"
- [x] Create a pure viral database for BLAST/bwa, alongside a mixed database for detecting false positives
- [x] Run paralell blastn query by population against the viral database
- [x] Run paralell blastn query by population against the mixed viral database to detect false positives
- [x] Run bwa mem query on one sample against the viral database
- [x] Obtain viral genome sizes and concatenate with blast result

# Pre-Processing


This set of script(s) below describe the main process that automates the entire pipeline. Elements will be added after they are executed seperately and validated. Initially, the plan was to run on a second round of blast on a mixed database comprised of the RefSeq viral genomes plus representative genomes of Archae, Bacteria, and some fungi, commercial vectors and plasmids, and some other full genomes (human, chimp, chicken and fruitfly) to detect false positives (17GB in total, while the average set of candidate reads was around 2.5MBs in size and upwards to 15MB in some minor instances). Due to time contraints however, bwa was chosen instead. To have a reasonable measure of comparison, "mixed" blast was performed on 100 samples only to ascertain a criteria for filtering non-significant reads from bwa on the level of criteria set in the paper (bit score of 190). The details are discussed below.

### Extraction and Merging of Unmapped Reads Into FASTA/Q Files

This script serves to extract samples of the 1000 genomes mapped and unmapped bam files by populations, and extracts/merges unmapped and chimeric pair reads. Fastq/a files are then created from the merged bam files in the final step. The merged bam file contain reads sorted by name, while the samtools fastq converter appends a /1 or /2 according to the read flag (forward/backward), and outputs two files split on reads /1 or /2.

Low complexity reads and reads where 50% of bases have a quality score (Phred) of 3 and below are filtered out before finally, the FASTQ file is converted to FASTA for blasting.


/samples/population/sample_ID/sample.fasta.
```bash
#!/bin/bash

pop=$1
parent="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/$pop/"
brits=$(curl -l $parent)
echo $pop
if [[ ! -e "/isdata/common/wbf326/samples/$pop" ]];then
	mkdir /isdata/common/wbf326/samples/$pop
fi

for sample in $brits; do
if [ ! -e "/isdata/common/wbf326/samples/$pop/$sample/$sample.fasta" ];then
	mkdir /isdata/common/wbf326/samples/$pop/$sample # sets up the structure of the directory such that each sample is placed alone by itself
	cd /isdata/common/wbf326/samples/$pop/$sample 
	page="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/$sample/alignment/" # sets up the ftp directory path by inserting the sample
	mapped=$(curl -l $page | grep -P ".+\.mapped.+\.bam$"); # sets up the ftp path to the mapped bam file and unmapped right below
	unmapped=$(curl -l $page | grep -P ".+\.unmapped.+\.bam$");
	echo "Mapped file is: $mapped"
	echo "Unmapped file is: $unmapped"

	echo "Downloading mapped and unmapped bam files.."
	curl $page$unmapped --output $sample.utemp.bam	 # downloads the unmapped bam file
	curl $page$mapped --output $sample.mtemp.bam # downloads the mapped bam file from the ftp list


	echo "Collecting mapped read info.."

	samtools view -@ 64 -c -F 0x4 $sample.mtemp.bam > $sample.stats

	echo "extracting mapped singletons.."
	#samtools view -@ 64 -h -F 0x4 -f 0x8 $sample.mtemp.bam -o $sample.mchi.txt # extracts definetly mapped with unmapped pair reads into .mchi

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

	rm $sample.mtemp.bam
	rm $sample.utemp.bam
	rm $sample.merged.bam
  	rm $sample.1.fastq
	rm $sample.2.fastq

fi
done

```

#### Running blastn On All Samples and Getting Candidate Reads

The purpose of this step is to run blast on the fasta files on the locally built viral database (see: Creation Of A Local Viral Database), filtering out all phage results and creating a set of candidate reads from the filtered blast results. The key steps here firstly removing duplicate reads from the sample.fasta file. This is necessary as it interferes with faidx indexing, without which subsetting the sample.fasta file with an ID list would have taken much longer than necessary. The second step is subsetting the fasta file with an ID list created from the viral blast results, thereby creating the \*.candidateReads.fasta file for each sample for the second round of analysis (mixed blast or bwa) 

Another key point deserve seperate mention; the fasta file is split into 56 seperate parts using fasplit, and 56 blasts are run in parallel which speeds up the process significantly, as it makes full use of avaiable cores on the server. Ironically, the slowest step was in extracting candidate reads from a list of IDs. Several attempts to speed up using grep and pyfaidx, and even at attempt at parallelizing pyfaidx proved unfruitful, but they are included and commented out in the script for the sake of propriety. Seqtk , a lightweight C++ program was atleast several magnitude faster than any of the aforementioned.

```bash
#!/bin/bash

pop=$1
parent="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/$pop/"
brits=$(curl -l $parent)
echo $brits
if [[ ! -e "/isdata/common/wbf326/viralblast/samples/$pop" ]];then
	mkdir /isdata/common/wbf326/viralblast/samples/$pop
fi



for sample in $brits; do
	if [[ -e "/isdata/common/wbf326/samples/$pop/$sample/$sample.fasta" && ! -e "/isdata/common/wbf326/viralblast/samples/$pop/$sample/$sample.candidateReads.fasta" ]]; then
	mkdir /isdata/common/wbf326/viralblast/samples/$pop/$sample
	cd /isdata/common/wbf326/samples/$pop/$sample
	# removing duplicates from fasta
	echo "removing duplicat reads.."
	cat $sample.fasta | /isdata/common/wbf326/./seqkit rmdup -o $sample.clean.fasta 
	# split the candidate read fasta files into 56 seperate files, create a list of all the files and feed it into GNU parallel to run the blast command on
	# the seperate files per available core
	/isdata/common/wbf326/./faSplit sequence $sample.clean.fasta 56 $sample.split
	ls *.fa > falist
	

	echo "blasting...$sample"
	
	# outputs blast results on the total viral database with a custom column format specified elsewhere
	cat falist | parallel --no-notice 'FASTA={};blastn -evalue 1e-10 -query $FASTA -db "/isdata/common/wbf326/viralblast/viraldb.fasta" -outfmt "6 qseqid sseqid evalue bitscore pident qcovhsp gapopen gaps positive length mismatch sacc slen qstart qend sstart send stitle" -out $FASTA.blasted'
	cat *.blasted > $sample.blast
	rm *.blasted
	mv $sample.blast /isdata/common/wbf326/viralblast/samples/$pop/$sample/$sample.blast
	
	echo "cleaning up.."
	# append each blast result into one big file
	# cleanup!except for the merged bam
	#rm *fa
	#rm falist

	cd /isdata/common/wbf326/viralblast/samples/$pop/$sample

	echo "removing phage results.."
	sed -i.bak '/phage/d' $sample.blast

	# create ID list of candidate viral reads with regex for extracting full fasta sequence
	echo "collecting candidate reads.."
	# splitting the file into 32 pieces, sorting, getting unique hits in parallel and merging into an intermediate file before sorting to find unique hits one and final addtional time
	# speeds up sorting considerably
	split --number=l/32 $sample.blast
	ls x* > xlist
	echo "collecting ID list.."
	cat xlist | parallel --no-notice cut -f 1 {} | sort | uniq >> intermediate
	cat intermediate | sort | uniq > IDs
 	rm intermediate
	rm x*
	echo "faidx using ID list.."
	#cat IDs | awk '{gsub("_","\\_",$0);$0="(?s)^>"$0".*?(?=\\n(\\z|>))"}1' > IDs2
	#xargs faidx -f /isdata/common/wbf326/samples/$pop/$sample/$sample.clean.fasta >> $sample.candidateReads.fasta
	#cat falist | parallel --no-notice -j32 'FASTA={};xargs faidx -f $FASTA < /isdata/common/wbf326/viralblast/samples/$pop/'"$sample"'/IDs' >> /isdata/common/wbf326/viralblast/samples/$pop/$sample/$sample.candidateReads.fasta
	rm *fa
	/isdata/common/wbf326/seqtk/./seqtk subseq /isdata/common/wbf326/samples/$pop/$sample/$sample.clean.fasta IDs > $sample.candidateReads.fasta
	#parallel --no-notice -j41 /home/wbf326/bin/pcre2grep -oM {} /isdata/common/wbf326/samples/$pop/$sample/$sample.clean.fasta :::: IDs2 >> $sample.candidateReads.fasta

 
	
fi 
done

```

#### Running blastn On Candidate Reads Against a Comprehensive Database

Candidate reads are split once again and a number of blasts are run in parallel, resulting in a final \*.mixedblast output. The difference between the parallelization process of this process and the previous one, is that the candidate reads are split into 1200 pieces, and blastn is run in batches of 56, once again making full use of all available cores. It was a necessary step as each individual blast run could take up as much as 20GB of space, if the files were split into 56 exact pieces instead, effectively using up available space and stalling. Running as much as 67200 blast runs on a sample did not increase overall running time of the whole process, leading one to surmise that the run time is linearly proportional to the query size.

```bash
#!/bin/bash

pop=$1
parent="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/$pop/"
brits=$(curl -l $parent)
echo $brits
if [[ ! -e "/isdata/common/wbf326/mixedblast/samples/$pop" ]];then
	mkdir /isdata/common/wbf326/mixedblast/samples/$pop
fi

for sample in $brits; do
	if [[ -d "/isdata/common/wbf326/viralblast/samples/$pop/$sample" &&  ! -e "/isdata/common/wbf326/mixed/samples/$pop/$sample/$sample.mixedblast" ]]; then

	mkdir /isdata/common/wbf326/mixed/samples/$pop/$sample
	cd /isdata/common/wbf326/viralblast/samples/$pop/$sample
	rm *.fa
	rm falist

	/isdata/common/wbf326/./faSplit sequence $sample.candidateReads.fasta 1200 $sample.split
	ls *.fa > falist

	echo "blasting candidate reads..$sample"
	# blasting each of the 64 split files seperately and cleanup!
	cat falist | parallel -j52 --no-notice 'FASTA={};blastn -evalue 1e-20 -query $FASTA -db "/isdata/common/wbf326/mixed/mixeddb.fasta" -outfmt "6 qseqid sseqid evalue bitscore sgi sacc slen qstart qend sstart send stitle" -out $FASTA.mixed'
	cat *.mixed > $sample.mixedblast
	rm *.mixed
	
	mv $sample.mixedblast /isdata/common/wbf326/mixed/samples/$pop/$sample/$sample.mixedblast
	rm *.fa
	rm falist
fi
done 

```

#### Post-Processing

This final step concatenate blast results, generates hit counts for quick inspection and collect stats (number of mapped reads to the human genome) in seperate files under their respective home directories. There is nothing noteworthy of mention aside from the R script used to filter out false positives (See: False Positive Tests), for the mixed blast run and appending sample names, and a simpler script for only appending sample names to the blast results for the viral blast run.


For the second round of blast (mixed blast):
```bash
#!/bin/bash
pop=$1
parent="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/$pop/"
brits=$(curl -l $parent)
echo $brits


#ls -d /isdata/common/wbf326/mixed/samples/$pop/*/ > $poplist # generate list of arguments for the R script in the form of directories


cd /isdata/common/wbf326/scripts
for sample in $brits; do
if [[ -d "/isdata/common/wbf326/viralblast/samples/$pop/$sample" && -e "/isdata/common/wbf326/viralblast/samples/$pop/$sample/$sample.blast" ]]; then
	echo "Running RScript on $sample.."
	Rscript --vanilla viralanalysis.R  "/isdata/common/wbf326/viralblast/samples/$pop/$sample/"

fi
done

for sample in $brits; do
if [[ -d "/isdata/common/wbf326/viralblast/samples/$pop/$sample" && -e "/isdata/common/wbf326/viralblast/samples/$pop/$sample/$sample.blast" ]]; then
	awk -F'\t' -vOFS='\t' -v sample="$sample" '{ $1 = sample"\t" $1 }1' < /isdata/common/wbf326/samples/$pop/$sample/$sample.stats >> /isdata/common/wbf326/mixed/samples/$pop/all.$pop.stats
fi
done


cd /isdata/common/wbf326/viralblast/samples/$pop

echo "gathering blast viral hit count.."
find .  -name '*.blast' -exec cat {} + > all.$pop.vblast
cut -f 18 all.$pop.vblast  |sort -r |  uniq -c | sort -n | sed -r 's/([0-9]) /\1\t/' > all.$pop.vHits

echo "gathering filtered (no phage) blast viral hit count.."
find .  -name '*.blast.bak' -exec cat {} + > all.$pop.vblastf
cut -f 18 all.$pop.vblastf  |sort -r |  uniq -c | sort -n | sed -r 's/([0-9]) /\1\t/' > all.$pop.vHitsf

```

And for the first round of blast (viral blast):

```bash
#!/bin/bash
pop=$1
parent="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/$pop/"
brits=$(curl -l $parent)
echo $brits


#ls -d /isdata/common/wbf326/mixed/samples/$pop/*/ > $poplist # generate list of arguments for the R script in the form of directories


cd /isdata/common/wbf326/scripts
for sample in $brits; do
if [[ -d "/isdata/common/wbf326/mixed/samples/$pop/$sample" && ! -e "/isdata/common/wbf326/mixed/samples/$pop/$sample/$sample.filtered" ]]; then
	echo "Running RScript on $sample.."
	Rscript --vanilla analysis.R  "/isdata/common/wbf326/mixed/samples/$pop/$sample/"

fi
done

for sample in $brits; do
if [[ -d "/isdata/common/wbf326/mixed/samples/$pop/$sample" && -e "/isdata/common/wbf326/mixed/samples/$pop/$sample/$sample.filtered" ]]; then
	awk -F'\t' -vOFS='\t' -v sample="$sample" '{ $1 = sample"\t" $1 }1' < /isdata/common/wbf326/samples/$pop/$sample/$sample.stats >> /isdata/common/wbf326/mixed/samples/$pop/all.$pop.stats
fi
done

cd /isdata/common/wbf326/mixed/samples/$pop

echo "gathering mixed filtered viral hit count.."
find .  -name '*.filtered' -exec cat {} + > all.$pop.mixedf
cut -f 6,12 all.$pop.mixedf  |sort -r |  uniq -c | sort -n | sed -r 's/([0-9]) /\1\t/' > all.$pop.mHitsf

cd /isdata/common/wbf326/viralblast/samples/$pop
```

## 2. Post Processing & Others



#### Creation Of A Local Viral Database:

Refseq viral genomes were downloaded from the NCBI database using the search criteria below:

```bash
Viruses[Organism] AND srcdb_refseq[PROP] NOT wgs[PROP] NOT cellular organisms[ORGN] NOT AC_000001:AC_999999[PACC] 
```

A perl script was then used to download the sequences based on that search criteria:

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

#### Creation of Mixed/Comprehensive Database for BLAST & BWA Indexing

The comprehensive database consisted of representaitve genomes for Fungi, Bacteria and Archae, alongwith UCSC human, chimp, mouse and chicken genomes, mixed with the refseq viral database.

Potential representative genomes from the assembly database:

Fungi(14/2)
Bacteria(1,543/1636)
Archaea(140/389)

The second number is the one reported in the paper. The search criteria for use on the BLAST database to extract the representative sequences is given here:

```bash
("Bacteria"[Organism] OR "Fungi"[Organism] OR "Archaea" [Organism]) AND (latest[filter] AND "complete genome"[filter] AND "representative genome"[filter]) 
```

Later, the UCSC genomes were downloaded individually, converted from 2bit format to fasta and screened for non-unique read names before concatenating all the above as one 17GB fasta file. To illustrate the first step, several entries are listed in different genomes as "chr(1 to 21) giving non-unique names for several conflicting entires (making it impossible to create a local database). To get around that issue, entries were modified with the command 'sed -i 's/old_entry/new_entry/g' where chrX for example, was renamed to chrX_chimp for all five genomes.

Finally, the BWA index was constructed for the mixed database using the bwtsw algorithm.

#### False Positives Tests


This R script filters false positives on the mixed blast run. It removes viral hits where there were hits to other organisms, with a bit score higher than or equal to the viral hit. So in principle, a read with multiple hits, one several to viral genomes and to other organisms, may have some viral hits removed if it has a lower bit score than any of the hits to non-viral organisms. The input is the directory where the \*.mixedblast file is, the output is placed in the same directory as \*.mixed.filtered.



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
  if (length(idx) ==0){ 
    next
  }  
  if (length(idx) == nrow(x)){
    z <- rbind(z,x[idx,])
    next
  }
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

#### Getting A Table of Viral Genome Sizes

This was a simple 3-step process: i. get viral genome IDs from the viral genome database, ii. get sequence lengths for each virus and iii. concatenate into a two column files. The second step is just a tiny bit more involved as it uses samtools faidx to generate index files from which the sequence lengths are extracted


```bash
grep "^>" viraldb.fasta | cut -d'|' -f3 | cut -d ' ' -f 2- > vNames
samtools faidx viraldb.fasta
cut -f2 viraldb.fasta.fai > vSizes
paste vNames vSizes

```

