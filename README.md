# Bionformatics MSc, Thesis: Viral Presence In The 1000 Genomes Database



# Introduction


The set of script below describe the main process that automates the entire pipeline. It does not cover the entire project in detail nor does it effectively document the various changes, fixes and evolutions that went on, but it captures the general spirit of it. For example, there were some unexpected and nasty complications such as BWA removing the forward and backward read identifiers (/1 and /2) which had to be corrected for, and some other fun workarounds such as the case in making a plot of integration sites over all human chromosomes in base-r graphics (a simple for loop and a custom points( ) function was a great fix). Both of these scripts are not included here for the sake of brevity.

Initially, the plan was to run on blast on a pure viral database for all 800 samples, followed by a second round of blast on a mixed database comprised of the RefSeq viral genomes plus representative genomes of Archae, Bacteria, and some fungi, commercial vectors and plasmids, and some other full genomes (human, chimp, chicken and fruitfly) to detect false positives (17GB in total, while the average set of candidate reads was around 2.5MBs in size and upwards to 15MB in some minor instances). Due to time contraints however, bwa was chosen instead in place of the second blast run. To have a reasonable measure of comparison, "mixed" blast was performed on 100 samples only to ascertain a criteria for filtering non-significant reads from bwa on the level of criteria set in the paper (bit score of 190). The details are discussed below, and the main script for generating the plots, figures and tables is included at the very end.

### Extraction and Merging of Unmapped Reads Into FASTA/Q Files

This script served as the official starting point of the project. Its main task is to download mapped and unammped reads in bam format from the 1000 genomes repository, and ultimately generate FASTA files for blasting, while doing general housekeeping like making directories for separate samples if needed. It downloads samples of the 1000 genomes mapped and unmapped bam files by population (this one is generically set to GBR), and extracts and merges unmapped and chimeric pair reads. Fastq/a files are then created from the merged bam files in the final step for later use. The merged bam file contains reads sorted by name, while the samtools fastq converter appends a /1 or /2 according to the read flag (forward/backward), and outputs two files split on reads /1 or /2.

Low complexity reads and reads where 50% of bases have a quality score (Phred) of 3 and below are filtered out before finally, the FASTQ file is converted to FASTA for blasting.
```bash

pop=$1
parent="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/$pop/"
brits=$(curl -l $parent)
echo $pop
if [[ ! -e "/isdata/common/wbf326/samples/GBR2" ]];then
	mkdir /isdata/common/wbf326/samples/GBR2
fi

for sample in $brits; do
if [ ! -e "/isdata/common/wbf326/samples/GBR2/$sample/$sample.fasta" ];then
	mkdir /isdata/common/wbf326/samples/GBR2/$sample # sets up the structure of the directory such that each sample is placed alone by itself
	cd /isdata/common/wbf326/samples/GBR2$sample 
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
	samtools view -@ 64 -h -F 0x4 -f 0x8 $sample.mtemp.bam -o $sample.mchi.txt # extracts definetly mapped with unmapped pair reads into .mchi

	echo "extracting unmapped singletons.."
	samtools view -@ 64 -h -f 0x4 $sample.mtemp.bam -o $sample.uchi.bam # extracts unmapped reads into .uchi

	echo "merging unmapped singletons with unmapped pairs.."
	samtools merge -@ 64 -n $sample.merged.bam $sample.uchi.bam $sample.utemp.bam # merges ALL mappped and unmapped reads

	echo "creating fastq files from unmapped reads.."
	samtools fastq -@ 64 -N -t -1 $sample.1.fastq -2 $sample.2.fastq $sample.utemp.bam # places forward/backward reads into respective fastq files
	
	
	# filters reads by quality and repeats as set in the paper and hopefully outputs a report in the correct place
	sga preprocess -p 0 --quality-filter=33 --no-primer-check --dust-threshold=2.5 --out=$sample.dusted.fastq $sample.1.fastq $sample.2.fastq

	echo "creating fasta from dusted fastq.."
	cat $sample.dusted.fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > $sample.fasta # transforms fastq into fasta using sed and trim

	rm $sample.utemp.bam
  	rm $sample.1.fastq
	rm $sample.2.fastq

fi

done
```

## BLAST and GNUParallel
#### Running blastn On All Samples and Getting Candidate Reads

With FASTA files in hand, the purpose of this step is to BLAST on the locally built viral database (see: Creation Of A Local Viral Database), filtering out all phage results and to generate a set of candidate reads from the filtered blast results. 

The key steps here are firstly removing duplicate reads from each sample FASTA file. The second step is generating a simple list of read IDs from the viral BLAST hits to be used in extracting FASTQ entries (the candidate reads) for use in the second round of BLAST/BWA (the "*.candidateReads.fasta" files for each sample) 

Another key point deserve separate mention; the fasta file is split into 48 separate parts using fasplit, and 48 blasts are run in parallel which speeds up the process significantly. Ironically, the slowest step was in extracting candidate reads using the list of IDs. Several attempts to speed up using grep and pyfaidx, and even at attempt at parallelizing pyfaidx proved unfruitful, but they are included and commented out in the script for the sake of propriety. Seqtk , a lightweight C++ program was atleast several magnitude faster than any of the aforementioned.

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
	echo "removing duplicate reads.."
	cat $sample.fasta | /isdata/common/wbf326/./seqkit rmdup -o $sample.clean.fasta 
	# split the candidate read fasta files into 48 seperate files, create a list of all the files and feed it into GNU parallel to run the blast command on
	# the seperate files per available core
	/isdata/common/wbf326/./faSplit sequence $sample.clean.fasta 48 $sample.split
	ls *.fa > falist
	
	echo "blasting...$sample"
	# outputs blast results on the total viral database with a custom column format specified elsewhere
	cat falist | parallel -j48--no-notice 'FASTA={};blastn -evalue 1e-10 -query $FASTA -db "/isdata/common/wbf326/viralblast/viraldb.fasta" -outfmt "6 qseqid sseqid evalue pident bitscore score gapopen gaps positive mismatch length qstart qend sstart send slen sacc stitle" -out $FASTA.blasted'
	cat *.blasted > $sample.blast
	rm *.blasted
	mv $sample.blast /isdata/common/wbf326/viralblast/samples/$pop/$sample/$sample.blast
	
	echo "cleaning up.."
	cd /isdata/common/wbf326/viralblast/samples/$pop/$sample

	echo "removing phage results.."
	sed -i.bak '/phage/d' $sample.blast
	
	echo "collecting candidate reads.."
	# splitting the file into 32 pieces, sorting, getting unique hits in parallel and merging into an intermediate file before sorting again to find unique hits one and final addtional time
	# speeds up sorting considerably
	split --number=l/32 $sample.blast
	ls x* > xlist
	echo "collecting ID list.."
	cat xlist | parallel --no-notice cut -f 1 {} | sort | uniq >> intermediate
	cat intermediate | sort | uniq > IDs
 	rm intermediate
	rm x*
	echo "Getting candidate reads FASTA/Q.."
	rm *fa
	/isdata/common/wbf326/seqtk/./seqtk subseq /isdata/common/wbf326/samples/$pop/$sample/$sample.dusted.fastq IDs > $sample.candidateReads.fastq 
	/isdata/common/wbf326/seqtk/./seqtk subseq /isdata/common/wbf326/samples/$pop/$sample/$sample.clean.fasta IDs > $sample.candidateReads.fasta
	
fi 
done


```

#### Running blastn On Candidate Reads Against a Comprehensive Database

The second round of BLAST was against a much larger database, so now it became an issue of space rather than speed, as each individual BLAST could take up to 20GBs of RAM.
Candidate reads are split once again and a number of blasts are run in parallel, resulting in a final output (" *.mixedblast"). The difference between the parallelization process of this process and the previous one, is that the candidate reads are split into 1200 pieces each, and blastn is run in batches of 48. That was, space was conserved without sacrificing speed, as BLAST run time is linearly proportional to the query size.


```bash
#!/bin/bash

pop=$1
parent="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/$pop/"
brits=$(curl -l $parent)
echo $brits
if [[ ! -d "/isdata/common/wbf326/mixed/samples/$pop" ]];then
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
	# blasting each of the 48 split files seperately and cleanup!
	cat falist | parallel -j48 --no-notice 'FASTA={};blastn -evalue 1e-20 -query $FASTA -db "/isdata/common/wbf326/bwaplus/mixedblastdb.fasta" -outfmt "6 qseqid sseqid evalue pident bitscore score gapopen gaps positive mismatch length qstart qend sstart send slen sacc stitle" -out $FASTA.mixed'
	cat *.mixed > $sample.mixedblast
	rm *.mixed
	
	mv $sample.mixedblast /isdata/common/wbf326/mixed/samples/$pop/$sample/$sample.mixedblast
	rm *.fa
	rm falist
fi
done 



```
## 2. Post Processing & Others

#### False Positives Tests

This final step concatenate blast results, generates hit counts for quick inspection and collect stats (number of mapped reads to the human genome) in seperate files under their respective home directories. There is nothing noteworthy of mention aside from the R script used to filter out false positives (See: False Positive Tests) for the mixed blast run and appending sample names, and a simpler script for only appending sample names to the blast results for the viral blast run.


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

cd /isdata/common/wbf326/viralblast/samples/$pop

echo "gathering filtered (no phage) blast viral hit count.."
find .  -name '*.blast' -exec cat {} + > all.$pop.vblastf
cut -f 18 all.$pop.vblastf  |sort -r |  uniq -c | sort -n | sed -r 's/([0-9]) /\1\t/' > all.$pop.vHitsf

echo "gathering full blast viral hit count.."
find .  -name '*.blast.bak' -exec cat {} + > all.$pop.vblast
cut -f 18 all.$pop.vblast  |sort -r |  uniq -c | sort -n | sed -r 's/([0-9]) /\1\t/' > all.$pop.vHits

```

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

#### Manipulating BWA Strings

The script below was made solely to manipulate and translate the CIGAR string and TAGs field from BWA results to ascertain the number of gaps, matches and mismatches for every hit, and to subsequently translate that into a score that is for all intents and purposes, BLAST-equivalent.

```r
#!/usr/bin/env Rscript
library(data.table)

files <- list.files(pattern = "*.bwamemplus$") # find all files with mixedblast extension
out <-sub(pattern = "(.*)\\..*$", replacement = "\\1", files) # get the sample name only


v <- fread("vRecords", header = F)

vlist <- v$V1
dt<- fread(files, header = F, sep = "\t")
colnames(dt) <-  c("query", "flag", "subject", "sstart", "Mapq", "cigar", "mismatches", "MD", "AS")
dt$subject <-gsub( ":.*", "",dt$subject, perl = T) # remove tags from secondary hits
dt$AS <- as.numeric (regmatches(dt$AS,regexpr("[0-9]+",dt$AS))) # convert score to a number
dt$mismatches <- as.numeric (regmatches(dt$mismatches,regexpr("[0-9]+",dt$mismatches))) # convert mismatches to a number
dt$MD <- gsub("MD:Z:", "", dt$MD) # remove MD tag from MD column
splitnums <- strsplit(dt$MD, split="[A-Z]") # get a vector of matches by removing any letters and splitting on them
dt$matches <-sapply(splitnums, function(x) sum(as.numeric (regmatches(x,regexpr("[0-9]+",x))))) # add up number of matches (including insertions/deletions)

# the following one-liners remove counts of Soft and Hard clips as well as Matches, then replaces Insertions and Deletions with a comma
# and splits the remaining string (of pure numbers of insertions and deletions) by the comma, before summing all the numbers to get a true
# count of gaps and real mismatches

dt$gaps <- sapply(dt$cigar, function(x) sum(as.numeric(strsplit (gsub("I|D",",",gsub("[0-9]+M|[0-9]+S|[0-9]+H", "", x)), split=",")[[1]]))[[1]])
dt$gapO <- sapply(dt$MD, function(x) sum(strsplit(x, split="")[[1]] == "^"))
dt$gapE <- dt$gaps - dt$gapO
dt$mismatches <- dt$mismatches - dt$gaps
dt$AS <- (dt$matches - (dt$mismatches*2 + dt$gaps*2.5))


relationship <- 1.835
lowerbound <- (190/relationship)-1
lowerbound <- (lowerbound/150) *100
dt<- dt[dt$matches>=lowerbound & ((dt$matches - lowerbound) - (dt$gaps*2.5 + dt$mismatches*2) >=0),]

ss <- split(dt,dt$query) # split the dataframe into a list of dataframes each with a single query ID and all of its hits
z <- data.frame(matrix(NA, nrow = nrow(dt), ncol = ncol(dt)))
cm <- 1
print("The number of reads is:")
print(length(ss))
for(i in 1:length(ss)){
  cat("\rRead ", i, " of ", length(ss)) 
  flush.console()
  
  x<-as.data.frame(ss[[i]]) # extract read i
  idx <- which (is.element(x$subject, vlist)) # find where read i hits match viruses and subset in t below
  if (length(idx) ==0){ 
    next
  }  
  if (length(idx) == nrow(x)){
    end <- cm + nrow(x[idx,]) - 1
    z[cm:end,] <- x[idx,]
    cm <- cm + nrow(x[idx,])
    
    next
  }
  t <- x[idx,]
  boo <- sapply(t$AS, function(z) all(z>x[-idx,9])) # add boolean checks where each viral hit scored higher than all other non-viral hits
  if(all(boo==F)){
    next
  }
  end <- cm + nrow(t[boo,]) - 1
  z[cm:end,] <- t[boo,]
  cm <- cm + nrow(t[boo,])
}


colnames(z) <-  c("query", "flag", "subject", "sstart", "Mapq", "cigar", "mismatches", "MD", "AS")
z<-z[complete.cases(z),] # remove NA values
z<-cbind(z,out)
z$subject <- v$V2[match(z$subject, v$V1)]  

write.table(z, paste(args,out,".filteredplus", sep=""), sep = "\t", col.names = F, row.names = F, quote = F) # output in the same directory with  sample.filtered extension
```
## The Main Script
The R script below does everything all in all. It does some heavy housekeeping firstly (like transforming problematic viral names into simpler names and such) and performs various other transformations for generating tabular summary statitsics and the plots and figures (base-r, NOT ggplot2). But the main bulk of the computation occurs in the bwaGitme and bwaReadme functions that calculated viral abundances (it was also the most enjoyable to code). The output was then collected and presented as part of the final report.

```r
library(scales)
library(RColorBrewer)
library(data.table)
library(stringr)
setwd("/home/adhamkmopp/thesis/results")

# loading sample/viral statistics


df <- read.table("stats/vNames3", sep = "\t", quote = "")
df$V2 <- sapply ( df$V2, function(x) gsub(", complete.*", "", x))
vsIdx <- grep("Influenza A.*segment\ [0-9]", df$V2)
a <-str_match(df$V2,"Influenza A.*segment\ [0-9]")
a <- a[complete.cases(a)]
df$V2[vsIdx] <-a

vsIdx <-grep("1934", df$V2)
df$V2[vsIdx[1]] <-"Influenza A virus (A/Puerto Rico/8/1934(H1N1))"
df$V1[vsIdx[1]] <- sum (df$V1[vsIdx])
df <- df[-vsIdx[-1],]

vsIdx <- grep("2009\\(H1N1\\)\\)", df$V2)
df$V2[vsIdx[1]] <-"Influenza A virus (A/California/07/2009(H1N1))"
df$V1[vsIdx[1]] <- sum (df$V1[vsIdx])
df <- df[-vsIdx[-1],]


vsIdx <- grep("A/goose/", df$V2)
df$V2[vsIdx[1]] <-"Influenza A virus (A/Guangdong/1/1996(H5N1))"
df$V1[vsIdx[1]] <- sum (df$V1[vsIdx])
df <- df[-vsIdx[-1],]

vsIdx <- grep("A/Korea/", df$V2)
df$V2[vsIdx[1]] <-"Influenza A virus (A/Korea/426/1968(H2N2))"
df$V1[vsIdx[1]] <- sum (df$V1[vsIdx])
df <- df[-vsIdx[-1],]

vsIdx <- grep("A/Hong Kong/.* segment", df$V2)
df$V2[vsIdx[1]] <-"Influenza A virus (A/Hong Kong/1073/99(H9N2))"
df$V1[vsIdx[1]] <- sum (df$V1[vsIdx])
df <- df[-vsIdx[-1],]


vsIdx <- grep("A/Goose/", df$V2)
df <- df[-vsIdx,]

vsIdx <- grep("Human papillomavirus", df$V2)
b <-str_match(df$V2,"Human papillomavirus")
b <- b[complete.cases(b)]
df$V2[vsIdx] <-b
vsIdx <- (-1* vsIdx[-1])
df <- df[vsIdx,]

vsIdx <- grep("Torque teno virus", df$V2)
c <-str_match(df$V2,"Torque teno virus")
c <- c[complete.cases(c)]
df$V2[vsIdx] <-c
vsIdx <- (-1* vsIdx[-1])
df <- df[vsIdx,]

vnames <- df$V2
vsizes <-df$V1
vgenomes <- split(vsizes,vnames)


df <- read.table("stats/humanv", sep = "\t", quote = "")

df$V2 <- sapply ( df$V2, function(x) gsub(", complete.*", "", x))
vsIdx <- grep("Influenza A.*segment\ [0-9]", df$V2)
a <-str_match(df$V2,"Influenza A.*segment\ [0-9]")
a <- a[complete.cases(a)]
df$V2[vsIdx] <-a

vsIdx <-grep("1934", df$V2)
df$V2[vsIdx[1]] <-"Influenza A virus (A/Puerto Rico/8/1934(H1N1))"
df <- df[-vsIdx[-1],]

vsIdx <- grep("2009\\(H1N1\\)\\)", df$V2)
df$V2[vsIdx[1]] <-"Influenza A virus (A/California/07/2009(H1N1))"
df <- df[-vsIdx[-1],]


vsIdx <- grep("A/goose/", df$V2)
df$V2[vsIdx[1]] <-"Influenza A virus (A/Guangdong/1/1996(H5N1))"
df <- df[-vsIdx[-1],]

vsIdx <- grep("A/Goose/", df$V2)
df <- df[-vsIdx,]

vsIdx <- grep("Human papillomavirus", df$V2)
b <-str_match(df$V2,"Human papillomavirus")
b <- b[complete.cases(b)]
df$V2[vsIdx] <-b
vsIdx <- (-1* vsIdx[-1])
df <- df[vsIdx,]

vsIdx <- grep("Torque teno virus", df$V2)
c <-str_match(df$V2,"Torque teno virus")
c <- c[complete.cases(c)]
df$V2[vsIdx] <-c
vsIdx <- (-1* vsIdx[-1])
df <- df[vsIdx,]

vsIdx <- grep("A/Korea/", df$V2)
df$V2[vsIdx[1]] <-"Influenza A virus (A/Korea/426/1968(H2N2))"
df$V1[vsIdx[1]] <- sum (df$V1[vsIdx])
df <- df[-vsIdx[-1],]

vsIdx <- grep("A/Hong Kong/.* segment", df$V2)
df$V2[vsIdx[1]] <-"Influenza A virus (A/Hong Kong/1073/99(H9N2))"
df$V1[vsIdx[1]] <- sum (df$V1[vsIdx])
df <- df[-vsIdx[-1],]


GBR.stats <- read.table("stats/all.GBR.stats")
colnames(GBR.stats) <- c("sampleIDs", "stats")
GBR.stats$sampleIDs <- as.character(GBR.stats$sampleIDs)


DK.stats <- read.table("stats/all.DK.stats")
colnames(DK.stats) <- c("sampleIDs", "stats")
DK.stats$sampleIDs <- as.character(DK.stats$sampleIDs)


PEL.stats <- read.table("stats/all.PEL.stats")
colnames(PEL.stats) <- c("sampleIDs", "stats")
PEL.stats$sampleIDs <- as.character(PEL.stats$sampleIDs)

YRI.stats <- read.table("stats/all.YRI.stats")
colnames(YRI.stats) <- c("sampleIDs", "stats")
YRI.stats$sampleIDs <- as.character(YRI.stats$sampleIDs)

CHB.stats <- read.table("stats/all.CHB.stats")
colnames(CHB.stats) <- c("sampleIDs", "stats")
CHB.stats$sampleIDs <- as.character(CHB.stats$sampleIDs)

IBS.stats <- read.table("stats/all.IBS.stats")
colnames(IBS.stats) <- c("sampleIDs", "stats")
IBS.stats$sampleIDs <- as.character(IBS.stats$sampleIDs)

CLM.stats <- read.table("stats/all.CLM.stats")
colnames(CLM.stats) <- c("sampleIDs", "stats")
CLM.stats$sampleIDs <- as.character(CLM.stats$sampleIDs)

JPT.stats <- read.table("stats/all.JPT.stats")
colnames(JPT.stats) <- c("sampleIDs", "stats")
JPT.stats$sampleIDs <- as.character(JPT.stats$sampleIDs)

LWK.stats <- read.table("stats/all.LWK.stats")
colnames(LWK.stats) <- c("sampleIDs", "stats")
LWK.stats$sampleIDs <- as.character(LWK.stats$sampleIDs)

ALL.stats <- read.table("stats/all.ALL.stats")
colnames(ALL.stats) <- c("sampleIDs", "stats")
ALL.stats$sampleIDs <- as.character(ALL.stats$sampleIDs)
ALL.stats$Pop <- rep("Nothing", length(ALL.stats$sampleIDs))
pops <- c("GBR", "CLM", "PEL", "IBS", "LWK", "YRI", "JPT", "CHB")

for (pop in pops){
  pop.stats <- get(paste(pop, ".stats",sep=""))
  ALL.stats$Pop[ALL.stats$sampleIDs %in% pop.stats$sampleIDs] <- pop
  
}

gotem <-fread("gotem")
# setting up a ggplot-like color function to rotate hues based on the number of colors queried
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

bwaReadme <- function(pop){
  # this function reads the bwa file and adds column names
  pop.bwaf <- read.table(paste("bwa/all.",pop, ".bwaplusf", sep=""), sep="\t")
  
  
  colnames(pop.bwaf) <- c("query", "flag", "sequence", "sstart", "Mapq", "cigar", "mismatches", "MD", "AS", 
                          "matches", "gaps", "gapO", "gapE", "sampleID")
  
  pop.bwaf$sampleID <- as.character(pop.bwaf$sampleID)
  pop.bwaf$sequence <- sapply(pop.bwaf$sequence, function(x) gsub(", complete.*", "", x)) 
  pop.bwaf$sstart <-as.numeric(pop.bwaf$sstart)
  pop.bwaf$matches <- as.numeric(pop.bwaf$matches)
  pop.bwaf$mismatches <- as.numeric(pop.bwaf$mismatches)
  
  
  vsIdx <- grep("Influenza A.*segment\ [0-9]", pop.bwaf$sequence)
  a <-str_match(pop.bwaf$sequence,"Influenza A.*segment\ [0-9]")
  a <- a[complete.cases(a)]
  pop.bwaf$sequence[vsIdx] <-a
  
  vsIdx <- grep("2009\\(H1N1\\)\\)",  pop.bwaf$sequence)
  pop.bwaf$sequence[vsIdx] <-"Influenza A virus (A/California/07/2009(H1N1))"
  
  vsIdx <-grep("A/goose/|A/Goose/",  pop.bwaf$sequence)
  pop.bwaf$sequence[vsIdx] <-"Influenza A virus (A/Guangdong/1/1996(H5N1))"
  
  vsIdx <- grep("A/Korea/", pop.bwaf$sequence)
  pop.bwaf$sequence[vsIdx] <-"Influenza A virus (A/Korea/426/1968(H2N2))"
  
  
  vsIdx <- grep("A/Hong Kong/.* segment", pop.bwaf$sequence)
  pop.bwaf$sequence[vsIdx] <- "Influenza A virus (A/Hong Kong/1073/99(H9N2))"
  
  
  vsIdx <-grep("1934", pop.bwaf$sequence)
  pop.bwaf$sequence[vsIdx] <-"Influenza A virus (A/Puerto Rico/8/1934(H1N1))"
  
  vsIdx <- grep("Human papillomavirus", pop.bwaf$sequence)
  b <-str_match(pop.bwaf$sequence,"Human papillomavirus")
  b <- b[complete.cases(b)]
  pop.bwaf$sequence[vsIdx] <-b
  
  
  vsIdx <- grep("Torque teno virus", pop.bwaf$sequence)
  c <-str_match(pop.bwaf$sequence,"Torque teno virus")
  c <- c[complete.cases(c)]
  pop.bwaf$sequence[vsIdx] <-c
  
  
  
  pop.bwaf$sequence[pop.bwaf$sequence=="Mus musculus mobilized endogenous polytropic provirus clone 15 truncated gag-pol polyprotein (gag) and envelope protein (env) genes"] <- "Mus musculus mobilized endogenous polytropic provirus clone 15"
  pop.bwaf$integrated <- rep(0, length(pop.bwaf$query))
  
  
  #pop.bwaf$flag <-sapply(pop.bwaf$flag, function(x) ifelse((x==0 | x==256),1,2))
  pop.bwaf$integrated[pop.bwaf$query %in% gotem$V1 & pop.bwaf$flag %in% gotem$V2] <-1
  
  
  
  assign(paste(pop,".bwaf",sep=""), pop.bwaf,envir = .GlobalEnv)
  
}

bwaGitme <- function(pop){
  
  # this function calculates adundances for viruses in the same manner as for mixed blast
  pop.stats <- get(paste(gsub(".f", "", pop), ".stats", sep=""))
  
  pop.bwa<-get(paste(pop,".bwaf",sep=""))
  
  viruslist<-by(pop.bwa, as.character(pop.bwa$sequence), function(m) return (m))
  vfound <- names(viruslist)
  hitlist <- matrix(data=0, nrow = length(pop.stats$sampleIDs), ncol=length(vfound))
  
  rownames(hitlist) <-pop.stats$sampleIDs
  colnames(hitlist) <- vfound
  final <- list(Virus=as.character(), Sample=as.character(), Abundance=as.numeric(), Reads=as.numeric(), Depth=as.numeric(), 
                Integration=as.numeric())
  
  
  for (virus in vfound){
    curr.virus <-viruslist[[virus]]
    rcounts <-table(curr.virus$sampleID)
    rfreq <- as.data.frame(rcounts)
    idx <- rcounts>0
    rset <-rcounts[idx]
    snames<-names(rset)
    for (sample in snames){
      
      curr.sample <- curr.virus[curr.virus$sampleID==sample,]
      curr.sample <-curr.sample[order(curr.sample$sstart),]
      start <-curr.sample$sstart
      
      depth <- sum(curr.sample$matches)
      
      hitlist[sample,virus]<-rset[sample]
      
      abundance <-  (2*rset[sample]/vgenomes[[virus]]) / (pop.stats[pop.stats$sampleIDs == sample,2]/3088286401)
      
      final$Virus<-c(final$Virus, virus)
      final$Sample<-c(final$Sample, sample)
      final$Abundance <- c(final$Abundance, abundance)
      final$Reads <-c(final$Reads, rset[sample])
      final$Depth <- c(final$Depth, depth/vgenomes[[virus]])
      
      final$Integration <- c(final$Integration, ifelse(any(curr.sample$integrated==1), 1, 0))
      
      
    }
    
  }
  
  hitlist <- as.data.frame(hitlist)
  totalhits <- colSums(hitlist)
  hitlist <- rbind(hitlist, totalhits) 
  
  
  row.names(hitlist) <-c(pop.stats$sampleIDs, "total" )
  
  final<-as.data.frame(final)
  print(sum(final$Integration==1))
  final<- final[-grep("phage|bacteriophage|retrovirus", final$Virus, ignore.case = T, perl = T, value = FALSE,
                      fixed = FALSE, useBytes = FALSE, invert = FALSE),]
  hitlist<- hitlist[, -grep("phage|bacteriophage|retrovirus", colnames(hitlist), ignore.case = T, perl = T, value = FALSE,
                            fixed = FALSE, useBytes = FALSE, invert = FALSE)]
  assign(paste(pop,".f.bwa.abundance",sep=""), final,envir = .GlobalEnv)
  assign(paste(pop,".f.bwa.hits",sep=""), hitlist,envir = .GlobalEnv)
  
}
pops <-c ("GBR", "PEL", "CHB", "YRI", "LWK", "JPT", "IBS", "CLM")

for (pop in pops){
  
  bwaReadme(pop)
  bwaGitme(pop)

}

bwaReadme("DK")
bwaGitme("DK")
ALL.bwaf <- DK.bwaf
ALL.f.bwa.abundance <- DK.f.bwa.abundance
ALL.f.bwa.hits <- DK.f.bwa.hits

ALL.bwaf <-fread("ALL.bwaplusf", sep = "\t")
ALL.bwaf$integrated <- rep(0, length(ALL.bwaf$query))
ALL.bwaf$flag <-sapply(ALL.bwaf$flag, function(x) ifelse((x==0 | x==256),1,2))
ALL.bwaf$integrated[ALL.bwaf$query %in% gotem$V1 & ALL.bwaf$flag %in% gotem$V2] <-1

goodfellas <-unique(ALL.bwaf$sampleID[ALL.bwaf$integrated ==1])
bwaGitme("ALL")
ALL.f.bwa.abundance$Pop<- ALL.stats$Pop[match(ALL.f.bwa.abundance$Sample,ALL.stats$sampleIDs)]


proportion <- rep(0,8)
names(proportion) <-c("GBR", "PEL", "CHB", "YRI", "LWK", "JPT", "IBS", "CLM")
for (pop in pops){
  
  x <-get(paste(pop, ".bwaf", sep=""))
  proportion[pop] <- sum(x[x$sequence=="Human adenovirus C",15]==1)/length(x$sequence)
}

######################################### PLINK AND SUCH ###################################################
getPalette = colorRampPalette( brewer.pal(n = 11, name = "Spectral"))

adeno <- ALL.f.bwa.abundance [ALL.f.bwa.abundance$Virus=="Human adenovirus C",]
ebv <- ALL.f.bwa.abundance [ALL.f.bwa.abundance$Virus=="Human gammaherpesvirus 4",]

ebv$Pop <-factor(ebv$Pop)
adeno$Pop <-factor(adeno$Pop)
acounts <-table(adeno$Pop)
ecounts <-table(ebv$Pop)
names(acounts) <- levels(adeno$Pop)
names(ecounts) <- levels(ebv$Pop)



png(file="adenobox.png",width=12,height=6,res=300,units="in", pointsize = 10)

clrs <- getPalette(4)
clrs <-rep(clrs,2)
clrsout <- sapply(as.numeric(adeno$Pop), function(x) clrs[x])
layout(matrix(c(1, 2), nrow=1, byrow=TRUE))
par(mar=c(4,4,4,4), family="Helvetica", xpd=F)
legendstuff <- sapply(levels(adeno$Pop), function(x) paste(x, acounts[x]))
boxplot(log(adeno$Abundance) ~ adeno$Pop ,whisklty = 1, staplelwd = 2,
        col=NULL, border = alpha(clrs,0.9), outbg=clrs, outpch=16,outline=F, 
        ylim=c(-8.5,7),horizontal=F, xlab="Population", main="Adenovirus C", ylab="Log-Abundances per 100,000/cells")

symbols(x=jitter(as.numeric(adeno$Pop), factor=0.75), y=log(adeno$Abundance), circles= rep(1,length(adeno$Pop)), inches= 0.035,
        ann=F,  fg="white", bg=clrsout, add=T)
legend("top",legendstuff, fill=clrs[1:8] ,ncol=8, 
       text.width=0.5, cex=0.65)


par(mar=c(4,4,4,4), family="Helvetica", xpd=F)
clrsout <- sapply(as.numeric(ebv$Pop), function(x) clrs[x])
legendstuff <- sapply(levels(ebv$Pop), function(x) paste(x, ecounts[x]))
boxplot(log(ebv$Abundance) ~ ebv$Pop ,whisklty = 1, staplelwd = 2,
        col=NULL, border = alpha(clrs,0.9), outbg=clrs, outpch=16,outline=F, 
        ylim=c(-8.5,7),horizontal=F, xlab="Population",  main="EBV",ylab="Log-Abundances per 100,000/cells")

symbols(x=jitter(as.numeric(ebv$Pop), factor=0.75), y=log(ebv$Abundance), circles= rep(1,length(ebv$Pop)), inches= 0.035,
        ann=F,  fg="white", bg=clrsout, add=T)
legend("top",legendstuff, fill=clrs[1:8] ,ncol=8, text.width=0.55, cex=0.65)

dev.off()

######################################### TABLES AND SUCH ###################################################
viruses <- colnames(ALL.f.bwa.hits)
individuals <- apply(ALL.f.bwa.hits, 2 , function(x) sum(x!=0))
percentages <- round((individuals/797)*100, digits=2)
reads.median <- apply(ALL.f.bwa.hits, 2, function(x) median(x))
reads.max <- apply(ALL.f.bwa.hits, 2, function(x) max(x))
abundance.median <- sapply(viruses, function(x) round(median (ALL.f.bwa.abundance[ALL.f.bwa.abundance$Virus==x,3]),digits=2))
abundance.max <- sapply(viruses, function(x) round(max(ALL.f.bwa.abundance[ALL.f.bwa.abundance$Virus==x,3]),digits=3))
depth.min <-sapply(viruses, function(x) round(min (ALL.f.bwa.abundance[ALL.f.bwa.abundance$Virus==x,5]),digits=3))
depth.max <-sapply(viruses, function(x) round(max (ALL.f.bwa.abundance[ALL.f.bwa.abundance$Virus==x,5]),digits=3))

table.bwa <- data.frame(viruses,individuals,percentages,reads.median, reads.max,abundance.median,abundance.max,depth.min,
                        depth.max) 
table.bwa <- table.bwa[order(table.bwa$percentages, decreasing = TRUE), ]
colnames(table.bwa) <- c("Virus", "Number of Individuals", "Percentage of Individuals", "Median", "Maximum","Median", 
                         "A", "B", "C")
rownames(table.bwa) <- NULL
write.table(table.bwa,file="table.bwa", sep="\t", quote=F, row.names = F)


matches <- which(viruses %in% df$V2 )
ALL.f.bwa.hits1 <- ALL.f.bwa.hits[, matches]
viruses <- colnames(ALL.f.bwa.hits1)
individuals <- apply(ALL.f.bwa.hits1, 2 , function(x) sum(x!=0))
percentages <- round((individuals/797)*100, digits=2)
reads.median <- apply(ALL.f.bwa.hits1, 2, function(x) median(x))
reads.max <- apply(ALL.f.bwa.hits1, 2, function(x) max(x))
abundance.median <- sapply(viruses, function(x) round(median (ALL.f.bwa.abundance[ALL.f.bwa.abundance$Virus==x,3]),digits=2))
abundance.max <- sapply(viruses, function(x) round(max(ALL.f.bwa.abundance[ALL.f.bwa.abundance$Virus==x,3]),digits=3))
depth.min <-sapply(viruses, function(x) round(min (ALL.f.bwa.abundance[ALL.f.bwa.abundance$Virus==x,5]),digits=3))
depth.max <-sapply(viruses, function(x) round(max (ALL.f.bwa.abundance[ALL.f.bwa.abundance$Virus==x,5]),digits=3))

table.bwa1 <- data.frame(viruses,individuals,percentages,reads.median, reads.max,abundance.median,abundance.max,depth.min,
                         depth.max) 
colnames(table.bwa1) <- c("Virus", "Number of Individuals", "Percentage of Individuals", "Median", "Maximum","Median", 
                          "A", "B", "C")

write.table(table.bwa1,file="table.human.bwa", sep="\t", quote=F, row.names = F)
clrs <- brewer.pal(n = 3, name = "Spectral")
getPalette = colorRampPalette( brewer.pal(n = 11, name = "Spectral"))


######################################### PREVALENCE AND BOXPLOTS ###################################################
png(file="bwaplot.png",width=12,height=10,res=300,units="in", pointsize = 12)
X11()
layout(matrix(c(1, 2), nrow=1, byrow=TRUE),widths=c(1,0.8), heights=c(1,1))
par(mar=c(6,16,3,0), family="Helvetica", xpd=F)
viruses <-as.character(table.bwa$Virus)
clrs <- getPalette(length(viruses))

b <-barplot(sort(table.bwa$`Percentage of Individuals`), cex.main=0.90,
            yaxt="n", horiz=T, col=clrs, xlab="Percentage of Individuals", cex.axis = 0.8, plot = F)
plot(0, type="n", xaxt="n", yaxt="n",xlab="", ylab="", xlim=c(0, 100), ylim=range(b), bty="n")
grid (NULL,NULL, lty = 6, col = "cornsilk2") 

barplot(sort(table.bwa$`Percentage of Individuals`), cex.main=0.90,
        yaxt="n", horiz=T, col=clrs, xlab="Percentage of Individuals", cex.axis = 0.6, add=T)
ylab <- as.character(table.bwa$Virus[order(table.bwa$`Percentage of Individuals`)])
text(par("usr")[3], seq(0.7,length(ylab)*1.2,by=1.2), adj= 1, xpd = TRUE,
     labels = ylab, cex=0.8)

axis(1, labels=seq(from=0,to=100, by=5), at=seq(from=0,to=100, by=5), cex.axis=0.6)

viruses <- levels(ALL.f.bwa.abundance$Virus)
amatrix <- matrix(0,ncol=1898, nrow = length(viruses))
rownames(amatrix) <- viruses
colnames(amatrix) <-row.names(ALL.f.bwa.hits)

vmatrix <- matrix(0,ncol=1898, nrow = length(viruses))
rownames(vmatrix) <- viruses
colnames(vmatrix) <-row.names(ALL.f.bwa.hits)


for (i in 1:nrow(ALL.f.bwa.abundance)) {
  amatrix[as.character(ALL.f.bwa.abundance[i,1]), as.character(ALL.f.bwa.abundance[i,2])] = as.numeric(ALL.f.bwa.abundance[i,3])
}

for (i in 1:nrow(ALL.f.bwa.abundance)) {
  vmatrix[as.character(ALL.f.bwa.abundance[i,1]), as.character(ALL.f.bwa.abundance[i,2])] = as.numeric(ALL.f.bwa.abundance[i,6])
}

allviruses <- as.vector(sapply(viruses, function(x) rep(x, ncol(amatrix))))
meltedmatrix <- (as.vector(t(amatrix)))
#meltedgotem <- (as.vector(t(vmatrix)))

idx <-which(meltedmatrix==0)
meltedmatrix <- meltedmatrix[-idx]
#meltedgotem <- meltedgotem[-idx]

meltedmatrix <- log(meltedmatrix)
allviruses <- allviruses[-idx]
allviruses <- factor(allviruses)


vfound <- length(unique(ALL.f.bwa.abundance$Virus))
matches <- which(ALL.f.bwa.abundance$Virus %in% viruses  )
ALL.f.bwa.abundance <- ALL.f.bwa.abundance[matches,]


ALL.f.bwa.abundance$Virus <-factor(ALL.f.bwa.abundance$Virus,levels(ALL.f.bwa.abundance$Virus)[match(table.bwa$Virus,levels(ALL.f.bwa.abundance$Virus))]) 

allviruses <-factor(allviruses,levels(allviruses)[match(table.bwa$Virus,levels(allviruses))])

ALL.f.bwa.abundance$Virus <- factor(ALL.f.bwa.abundance$Virus, 
                                    levels(ALL.f.bwa.abundance$Virus)[order(table.bwa$`Percentage of Individuals`)])
ALL.f.bwa.abundance$Virus <- factor(ALL.f.bwa.abundance$Virus) 

allviruses <- factor(allviruses, 
                                    levels(allviruses)[order(table.bwa$`Percentage of Individuals`)])
allviruses <- factor(allviruses) 

#idx1<- meltedgotem ==1
bloods <-alpha(clrs[as.numeric(factor(allviruses))],0.6)
#cs <-sapply(seq_along(allviruses), function(i) ifelse(idx1[i], "#e34a33", bloods[i]))

par(mar=c(6,0,3,1), family="Helvetica", xpd=F)


plot(0, type="n", xaxt="n", yaxt="n",xlab="", ylab="", xlim=c(min(meltedmatrix), max(meltedmatrix)), ylim=c(1,12), bty="n")
grid (NULL,NULL, lty = 6, col = "cornsilk2") 


symbols(y=jitter(as.numeric(allviruses), factor=0.75), x=meltedmatrix, circles= rep(1,length(allviruses)), inches= 0.035,
        ann=F,  fg=alpha("black",0.8), bg=bloods, add=T)


boxplot(meltedmatrix ~ allviruses , yaxt="n",whisklty = 1, staplelwd = 2,
        col=alpha(clrs, 0.8), border = alpha("black",0.8), outline=F, 
        ylim=c(min(meltedmatrix),max(meltedmatrix)), horizontal=T, xlab="Abundances per 100,000/cells", add=T)


dev.off()
gc()


###############################################################################################################################
########################################### PREVALENCE AND BOXPLOTS (HUMAN) ###################################################

png(file="bwahumanplots.png",width=12,height=10,res=300,units="in", pointsize = 12)
layout(matrix(c(1, 2), nrow=1, byrow=TRUE),widths=c(1,0.8), heights=c(1,1))


par(mar=c(6,18,3,0), family="Helvetica", xpd=F)
viruses <-as.character(table.bwa1$Virus )
clrs <- getPalette(length(viruses))

b <-barplot(sort(table.bwa1$`Percentage of Individuals`), cex.main=0.90,
            yaxt="n", horiz=T, col=clrs, xlab="Percentage of Individuals", cex.axis = 0.8, plot = F)
plot(0, type="n", xaxt="n", yaxt="n",xlab="", ylab="", xlim=c(0, 100), ylim=range(b), bty="n")
grid (NULL,NULL, lty = 6, col = "cornsilk2") 

barplot(sort(table.bwa1$`Percentage of Individuals`), cex.main=0.90,
        yaxt="n", horiz=T, col=clrs, xlab="Percentage of Individuals", cex.axis = 0.6, add=T)
ylab <- as.character(table.bwa1$Virus[order(table.bwa1$`Percentage of Individuals`)])
text(par("usr")[3], seq(0.7,length(ylab)*1.2,by=1.2), adj= 1, xpd = TRUE,
     labels = ylab, cex=0.8)

axis(1, labels=seq(from=0,to=10, by=1), at=seq(from=0,to=10, by=1), cex.axis=0.6)


matches <- which(ALL.f.bwa.abundance$Virus %in% df$V2  )
ALL.f.bwa.abundance1 <- ALL.f.bwa.abundance[matches,]
vfound <- length(unique(ALL.f.bwa.abundance1$Virus))
ALL.f.bwa.abundance1$Virus <- factor(ALL.f.bwa.abundance1$Virus)

viruses <- levels(ALL.f.bwa.abundance1$Virus)
amatrix <- matrix(0,ncol=1898, nrow = length(viruses))
rownames(amatrix) <- viruses
colnames(amatrix) <-row.names(ALL.f.bwa.hits)

vmatrix <- matrix(0,ncol=1898, nrow = length(viruses))
rownames(vmatrix) <- viruses
colnames(vmatrix) <-row.names(ALL.f.bwa.hits)



for (i in 1:nrow(ALL.f.bwa.abundance1)) {
  amatrix[as.character(ALL.f.bwa.abundance1[i,1]), as.character(ALL.f.bwa.abundance1[i,2])] = as.numeric(ALL.f.bwa.abundance1[i,3])
}

for (i in 1:nrow(ALL.f.bwa.abundance1)) {
  vmatrix[as.character(ALL.f.bwa.abundance1[i,1]), as.character(ALL.f.bwa.abundance1[i,2])] = as.numeric(ALL.f.bwa.abundance1[i,6])
}


allviruses <- as.vector(sapply(viruses, function(x) rep(x, ncol(amatrix))))
meltedmatrix <- (as.vector(t(amatrix)))
meltedgotem <- (as.vector(t(vmatrix)))

idx <-which(meltedmatrix==0)
meltedmatrix <- meltedmatrix[-idx]
meltedgotem <- meltedgotem[-idx]

meltedmatrix <- log(meltedmatrix)
allviruses <- allviruses[-idx]
allviruses <- factor(allviruses)



ALL.f.bwa.abundance1$Virus <-factor(ALL.f.bwa.abundance1$Virus,levels(ALL.f.bwa.abundance1$Virus)[match(table.bwa1$Virus,levels(ALL.f.bwa.abundance1$Virus))]) 
allviruses <-factor(allviruses,levels(allviruses)[match(table.bwa1$Virus,levels(allviruses))])

ALL.f.bwa.abundance1$Virus <- factor(ALL.f.bwa.abundance1$Virus, 
                                     levels(ALL.f.bwa.abundance1$Virus)[order(table.bwa1$`Percentage of Individuals`)])
allviruses <- factor(allviruses, 
                     levels(allviruses)[order(table.bwa1$`Percentage of Individuals`)])
allviruses <- factor(allviruses) 

idx1<- meltedgotem ==1
bloods <-alpha(clrs[as.numeric(factor(allviruses))],0.6)
cs <-sapply(seq_along(allviruses), function(i) ifelse(idx1[i], "#e34a33", bloods[i]))

par(mar=c(6,0,3,1), family="Helvetica", xpd=F)

plot(0, type="n", xaxt="n", yaxt="n",xlab="", ylab="", xlim=c(min(meltedmatrix), max(meltedmatrix)), ylim=c(1,7), bty="n")
grid (NULL,NULL, lty = 6, col = "cornsilk2") 



symbols(y=jitter(as.numeric(allviruses), factor=0.75), x=meltedmatrix, circles= rep(1,length(allviruses)), inches=0.035,
        ann=F,  fg=alpha("black",0.6), bg=cs, add=T)

boxplot(meltedmatrix ~ allviruses , yaxt="n",
        col=alpha(clrs, 0.6), border = alpha("black",0.8), outline=F, whisklty = 1, staplelwd = 2,
        ylim=c(min(meltedmatrix),max(meltedmatrix)), horizontal=T, xlab="Abundances per 100,000/cells", add=T)

dev.off()
gc()



############################################### COVERAGE PLOTS (HUMAN) #####################################################
viruses <-as.character(table.bwa1$Virus )[c(1,2,3,5,6,4)]
getPalette = colorRampPalette( brewer.pal(n = 6, name = "Paired"))
clrs<- getPalette(6)
png(file="bwacovplots",width=12,height=6,res=300,units="in", pointsize = 12)
layout(matrix(c(1, 2, 3,
                4, 5, 6
), nrow=3, ncol=2, byrow=TRUE),widths=c(1.25,1.25,1,1), heights=c(1,1,1,1))

for (i in 1:length(viruses) ){
  
  ifelse(i==1,
         par(mar=c(3,5,3,3),xpd=T, family ="Helvetica"),
         par(mar=c(3,5,3,3),family ="Helvetica"))
  idx <- ALL.bwaf$sequence == viruses[i]
  x <- ALL.bwaf$sstart[idx]
  hist(x, main= paste(viruses[i],"\nALL", sep=" "),  lty="blank", col=clrs[i],cex.main=0.95, xlim=c(0,vgenomes[[viruses[i]]]), 
       breaks=vgenomes[[viruses[i]]]/100, xlab="Genomic Position", ylab = "Read Frequency", freq=T)
  if (i==1) {
    corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
    text(x = corners[1]-20000, y = corners[4], "A", cex=1.5)
    
  }
  
}

dev.off()
viruses <-viruses[c(1,3)]

png(file="integrationcoverage.png",width=14,height=8,res=300,units="in", pointsize = 10)
clrs <- c("#E495A5", "#ABB065", "#39BEB1", "#ACA4E2")
layout(matrix(c(1, 2), nrow=2, ncol=1, byrow=TRUE),widths=c(1,1,1,1), heights=c(1,1,1,1))
par(mar=c(3,5,3,3),xpd=T, family ="Helvetica")
idx <- ALL.bwaf$sequence == viruses[1]
x <- ALL.bwaf$sstart[idx]
hist(x, main= paste(viruses[1],"\nALL", sep=" "),  lty="blank", col=clrs[3],cex.main=0.95, xlim=c(0,vgenomes[[viruses[1]]]), 
       breaks=vgenomes[[viruses[1]]]/100, xlab="Genomic Position", ylab = "Read Frequency", freq=T)

idx <- ALL.bwaf$sequence == viruses[1] & ALL.bwaf$integrated ==1
x <- ALL.bwaf$sstart[idx]
hist(x,  lty="blank", col="red",cex.main=0.95, 
     breaks=vgenomes[[viruses[1]]]/100, yaxt='n', xaxt='n' , freq=T, add=T)

idx <- ALL.bwaf$sequence == viruses[2]
x <- ALL.bwaf$sstart[idx]
hist(x, main= paste(viruses[2],"\nALL", sep=" "),  lty="blank", col=clrs[3],cex.main=0.95, xlim=c(0,vgenomes[[viruses[2]]]), 
     breaks=vgenomes[[viruses[2]]]/100, xlab="Genomic Position", ylab = "Read Frequency", freq=T)

idx <- ALL.bwaf$sequence == viruses[2] & ALL.bwaf$integrated ==1
x <- ALL.bwaf$sstart[idx]
hist(x,  lty="blank", col="red",cex.main=0.95, 
     breaks=vgenomes[[viruses[1]]]/100, yaxt='n', xaxt='n' , freq=T, add=T)

dev.off()



############################################# TRANSFORMATION FOR PLINK ############################################
g <- read.table("./GWAS/merged.fam", sep = ' ', header = F)
samples <- ALL.stats$sampleIDs
idx <-g$V1 %in% samples
headass <- g[idx,]
headass$V1 <-factor(headass$V1)
headass$V2 <-factor(headass$V2)
herpes <- amatrix["Human gammaherpesvirus 4",]
herpes <- herpes[-length(herpes)]
herpes <- herpes[herpes!=0]
herpes <- qnorm((rank(herpes) - 1/2)/length(herpes))

headass$V7 <- as.character(headass$V6)
headass$V8 <- as.character(headass$V6)
rownames(headass) <- headass$V1
headass[names(herpes),7] <-as.character(herpes)


herpes <- amatrix["Human adenovirus C",]
herpes <- sapply(herpes, function(x) ifelse(x==0, 1,2))

headass[names(herpes),8] <-as.character(herpes)

c <- headass[1:791,]

write.table(c, "sub.fam", col.names = F, row.names = F, quote = F)
ids <- sapply(seq_along(ALL.bwaf$query), function(x) ifelse((ALL.bwaf$flag[x]==0 | ALL.bwaf$flag[x]==256)
                                                            , paste(ALL.bwaf$query[x],"/1", sep="")
                                                            , paste(ALL.bwaf$query[x],"/2", sep="")))

ids <- unique(ids)
write.table(ids, "ids.txt", col.names = F, row.names = F, quote = F)





# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
ir.pca <- prcomp(amatrix,
                 center = TRUE,
                 scale. = TRUE) 
deuces <- NULL
for (i in 1:790){
  deuces<-c(deuces, (all(amatrix[,i]==0)))
}
amatrix <- amatrix[,!deuces]

ir.pca <- prcomp(t(amatrix[8:30,]),
                 center = TRUE,
                 scale. = TRUE) 

br <- factor(ALL.stats$Pop[match(names(ir.pca$x[1:772,1]), ALL.stats$sampleIDs)])

png(file="pca.ong",width=12,height=6,res=300,units="in", pointsize = 10)

plot(x=ir.pca$x[1:772,1], y=ir.pca$x[1:772,2], col=clrs[as.numeric(br)])
dev.off()

kruskal.test(ebv$Abundance, ebv$Pop)

dunn.test (x, g=NA, method=p.adjustment.methods, kw=TRUE, label=TRUE,
           wrap=FALSE, table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05, altp=FALSE)

ad <-pairwise.wilcox.test(adeno$Abundance, adeno$Pop,p.adj = "bonf")
eb <- pairwise.wilcox.test(ebv$Abundance, ebv$Pop,p.adj = "bonf")
write.table(eb$p.value, "ebtable",quote=F, row.names = T, col.names = T, sep="\t")
write.table(ad$p.value, "adtable",quote=F, row.names = T, col.names = T, sep="\t")

dunn.test (adeno$Abundance, g=adeno$Pop, method="bonferroni", kw=TRUE, label=TRUE,
           wrap=FALSE, table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05, altp=FALSE)
integrads <- data.frame(matrix(ncol = 2, nrow = 0))

colnames(integrads) <- c("ints", "pop")
for (pop in pops){
  
  x <-get(paste(pop, ".bwaf", sep=""))
  y <- as.vector(by(x$integrated, x$sampleID, function(m){
    
    return(sum(m))
  } ))
  z <- rep(pop, length(y))
  t <- cbind(y,z)
  integrads <- rbind(t, integrads)
}


integrads$z <- as.character(integrads$z)
integrads$z <- factor(as.character(integrads$z))
integrads$y <- as.numeric(levels(integrads$y))[integrads$y] 
dunn.test (integrads$y, g=integrads$z, method="bonferroni", kw=TRUE, label=TRUE,
           wrap=FALSE, table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05, altp=FALSE)
