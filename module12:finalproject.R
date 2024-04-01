library(GenomicAlignments)


chip = readGAlignments("/Users/sydneyrattey/Desktop/groupproject/SRR4105760.bam")
View(chip)
head(chip)

control = readGAlignments("/Users/sydneyrattey/Desktop/groupproject/SRR4105763.bam")

head(control)

?GenomicAlignments

available.genomes
getBSgenome()



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsamtools")

library(BiocManager)

library(Rsamtools)
browseVignettes("Rsamtools")
control <- scanBamHeader('/Users/sydneyrattey/Desktop/groupproject/SRR4105763.bam')
controlstats <- idxstatsBam('/Users/sydneyrattey/Desktop/groupproject/SRR4105763.bam')
head(control)
View(control)


chip <- scanBamHeader('/Users/sydneyrattey/Desktop/groupproject/SRR4105760.bam')
chipstats <- idxstatsBam('/Users/sydneyrattey/Desktop/groupproject/SRR4105760.bam')
chip 
chipstats

class(chip)
class(chipstats)
isS4(chipstats)
isS4(chip)

?system.file()
scanBam


open.BamFile
input <- BamFile('/Users/sydneyrattey/Desktop/groupproject/SRR4105760.bam', 
        index='/Users/sydneyrattey/Desktop/groupproject/SRR4105760.bam.bai')
inputstats <- idxstatsBam(input)
View(input)
head(input)
summary(input)
inputstats
class(inputstats)
inputstats$seqnames

seqinfo(input)
countBam(input)
quickBamFlagSummary(input)

ipfile <-BamFile('/Users/sydneyrattey/Desktop/groupproject/SRR4105763.bam', 
                 index='/Users/sydneyrattey/Desktop/groupproject/SRR4105763.bam.bai')
ipfilestats <- idxstatsBam(ipfile)
ipfilestats
ipfilestats[1:4,2:18]

class(ipfile)


ipfile$seqnames

available.genomes()
library(BSgenome)
available.genomes()

library(BiocManager)
install("BSgenome.Scerevisiae.UCSC.sacCer1")
install("BSgenome.Scerevisiae.UCSC.sacCer3")
yeastgenome <- getBSgenome("sacCer3")
seqlength = seqlengths(yeastgenome)
seqlength

seqinfo(yeastgenome)
yeastgenome$chr10
yeastgenome
class(seqlength)
seqlength

class(inputstats)

inputstats$mapped[2:18]


MappedReadsPerChromosome <- data.frame("Chromosome Length"=seqlength,"Mapped Reads" = inputstats$mapped[2:18])

MappedReadsPerChromosome

Mappedreadsperchrom <- function(x){
  data.frame("Chromosome Length"=seqlength,"Mapped Reads" = inputstats$mapped[2:18])
}

?data.frame

library(normr)
peakfit = enrichR(treatment = "/Users/sydneyrattey/Desktop/groupproject/SRR4105763.bam",
control= "/Users/sydneyrattey/Desktop/groupproject/SRR4105760.bam",
genome="sacCer3",verbose=FALSE)
summary(peakfit)

install.packages("ShortReadQ")
library(ShortRead)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ShortRead")

library(ShortRead)
?ShortRead
browseVignettes("ShortRead")


data <- FastqStreamer("/Users/sydneyrattey/downloads/SRR4105760.fastq.gz")
View(data)
countFastq("/Users/sydneyrattey/downloads/SRR4105760.fastq.gz")
readFastq("/Users/sydneyrattey/downloads/SRR4105760.fastq.gz")



alphabetByCycle("/Users/sydneyrattey/downloads/SRR4105760.fastq.gz", alphabet)

alphabetFrequency("/Users/sydneyrattey/downloads/SRR4105760.fastq.gz")
FastqSampler("/Users/sydneyrattey/downloads/SRR4105760.fastq.gz")

install.packages("seqinr")
library(seqinr)

fasta <- read.fasta("/Users/sydneyrattey/downloads/M77821.1.fasta")
write.table(fasta, "new")
head(fasta)
class(fasta)
as.data.frame(fasta)
fasta
class(fasta)
as.matrix(fasta)
fasta[2,1]
DNAString("/Users/sydneyrattey/downloads/M77821.1.fasta")

my_dataframe <- as.data.frame(fasta)

my_dataframe
fastalength <- length(my_dataframe[i,1])
my_dataframe[5,1]
my_dataframe[i=="a",1]
A <- length(my_dataframe[i=="a",1])
A
(A/fastalength)*100

my_dataframe[i=="t",1]
T <- length(my_dataframe[i=="t",1])


my_dataframe[i=="g",1]
G <- length(my_dataframe[i=="g",1])

my_dataframe[i=="c",1]
C <- length(my_dataframe[i=="c",1])

percentA <- function(x){
  my_dataframe <- as.data.frame(x)
  fastalength <- length(my_dataframe[i,1])
  A <- length(my_dataframe[i=="a",1])
  (A/fastalength)*100
}

percentA(fasta)

base(my_dataframe)

is.SeqFastadna(fasta)
s

fasta[grep("[a]", fasta)]

length(fasta)
nchar(fasta)
nchar(fasta[which(fasta = "a")])

as.SeqFastadna(fasta)
as.SeqFastadna(fasta, name="test")
summary.SeqFastadna(fasta)


?grep

practice <- data.frame(fasta)
class(practice)
practice
nchar(practice[which(practice = "a")])
length(practice)
write.table(practice, "test")

for (x in practice[x,1]){
  print(x)
}

G(practice)

base <- function(x){
  y = length(x[which(x == "a")])
  y
  z = length(x)
  z
  (y/z)*100
}
practice[4]
nrows(practice)

#have user copy and paste a fasta file
#insert commas into file

length(fasta[which(fasta == t)])
fasta[which(fasta=="a") ]

fasta2 <- data.frame(fasta)
fasta2
fasta2[i,1]
head(fasta2)
length(fasta2)
class(fasta2)
length(fasta2)
fasta2[i,1]
new2
length(new2)
class(new2)

fasta2[grep("[a]", fasta2)]


nchar(fasta ="a")

nchar(fasta == "a")
nchar(fasta == "t")
nchar(fasta == "c")
nchar(fasta == "g")

nchar(fasta)
length(fasta)
width(fasta)
count(fasta)

vector <- data.frame("a","a","a","c","c","c","t","g","g")
vector
length(vector[which(vector == "g")])
length(vector[which(vector == "t")])
length(vector)


vector2 <- data.frame("a","a","a","c","c","c","t","g","g")
vector2
length(vector2)
nchar(vector2)
#calculate the percentage of g in the values
  y = length(vector[which(vector == "g")])
  y
  z = length(vector)
  z
  (y/z)*100
  
  
  #create a function to calculate percentage of g
  G <- function(x){
  y = length(x[which(x == "g")])
  y
  z = length(x)
  z
  (y/z)*100
  }
  
  G(vector)

  fpkm <- function(x){
    x/geneLen/sum(x)*10^9
  }
  
  
  
  vector
class(vector)

random <- data.frame("GAATTCCGCTTAGAATGAACATTATTGTTACTGCTATGGAGATTGTTGTAGAGCCGGTGC
ATATTTTCCTCATCAGAACCTCGACGAGGATCATGTTCATTAGAGTGCATAATATGGTCT
AATGAGTGGTTCGGTGACACATCACTCAGGTTTCTATCTACAAAGGATTTGGGCAATGCA
TGCTTCCTCTCCGCCTCCTTGCCGTTGGGAGTTCTTCGTTTTTCTCTTTGTAAACTGCTA")

nchar(random)
nchar(random[where(random="G")])



random
nchar(random[which(random == G)])

for( i in fasta){
  print(nchar(fasta[which(fasta = "t")]))
}

?data.frame


?List
  
  for(i in 1:ncol(dfReads)){
    print(fpkm(dfReads[,2], df$Length))
  }

fpkm <- function(x, geneLen = df$Length){
  x/geneLen/sum(x)*10^9
}



