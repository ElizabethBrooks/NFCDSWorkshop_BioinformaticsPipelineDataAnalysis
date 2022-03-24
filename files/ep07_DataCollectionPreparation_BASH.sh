#export PATH=$PATH:$PWD/sratoolkit.2.11.2-mac64/bin

prefetch SRR8288561 SRR8288562 SRR8288563 SRR8288564 SRR8288557 SRR8288560 SRR8288558 SRR8288559 SRR8288565 SRR8288566 SRR8288567 SRR8288568

fastq-dump --gzip SRR8288561; fastq-dump --gzip SRR8288562; fastq-dump --gzip SRR8288563; fastq-dump --gzip SRR8288564; fastq-dump --gzip SRR8288557; fastq-dump --gzip SRR8288560
fastq-dump --gzip SRR8288558; fastq-dump --gzip SRR8288559; fastq-dump --gzip SRR8288565; fastq-dump --gzip SRR8288566; fastq-dump --gzip SRR8288567; fastq-dump --gzip SRR8288568

#brew install fastqc

fastqc SRR8288561.fastq.gz --extract

#wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh

#bash Miniconda3-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda

#source /Users/bamflappy/miniconda/bin/activate

#conda init zsh

#conda list

#conda install -c bioconda gffread

gffread -E -F -T Tribolium_castaneum.gff3 -o Tribolium.gtf

#export PATH=$PATH:$PWD/hisat2-2.2.1

hisat2-build Tribolium_castaneum.genome.fa TriboliumBuild

##cntrl samples 4h
hisat2 -q -x TriboliumBuild -U SRR8288561.fastq.gz -S SRR8288561_accepted_hits.sam --summary-file SRR8288561_alignedSummary.txt
hisat2 -q -x TriboliumBuild -U SRR8288562.fastq.gz -S SRR8288562_accepted_hits.sam --summary-file SRR8288562_alignedSummary.txt
hisat2 -q -x TriboliumBuild -U SRR8288563.fastq.gz -S SRR8288563_accepted_hits.sam --summary-file SRR8288563_alignedSummary.txt

##cntrl samples 24h
hisat2 -q -x TriboliumBuild -U SRR8288558.fastq.gz -S SRR8288558_accepted_hits.sam --summary-file SRR8288558_alignedSummary.txt
hisat2 -q -x TriboliumBuild -U SRR8288567.fastq.gz -S SRR8288567_accepted_hits.sam --summary-file SRR8288567_alignedSummary.txt
hisat2 -q -x TriboliumBuild -U SRR8288568.fastq.gz -S SRR8288568_accepted_hits.sam --summary-file SRR8288568_alignedSummary.txt

##treat samples 4h
hisat2 -q -x TriboliumBuild -U SRR8288564.fastq.gz -S SRR8288564_accepted_hits.sam --summary-file SRR8288564_alignedSummary.txt
hisat2 -q -x TriboliumBuild -U SRR8288557.fastq.gz -S SRR8288557_accepted_hits.sam --summary-file SRR8288557_alignedSummary.txt
hisat2 -q -x TriboliumBuild -U SRR8288560.fastq.gz -S SRR8288560_accepted_hits.sam --summary-file SRR8288560_alignedSummary.txt

##treat samples 24h
hisat2 -q -x TriboliumBuild -U SRR8288559.fastq.gz -S SRR8288559_accepted_hits.sam --summary-file SRR8288559_alignedSummary.txt
hisat2 -q -x TriboliumBuild -U SRR8288565.fastq.gz -S SRR8288565_accepted_hits.sam --summary-file SRR8288565_alignedSummary.txt
hisat2 -q -x TriboliumBuild -U SRR8288566.fastq.gz -S SRR8288566_accepted_hits.sam --summary-file SRR8288566_alignedSummary.txt
