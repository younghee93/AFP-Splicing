#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=24
#SBATCH --mem=64000
#SBATCH -o ./tmp/chr.rMATs.o.%j
#SBATCH --err ./tmp/chr.rMATs.e.%j


module load tophat
module load bowtie2
module load samtools/1.2
module load python/2.7.11
ml star/2.5

mkdir -p $3

echo $3
python ./Tools/rMATS.3.2.5/RNASeq-MATS.py -b1 $1 -b2 $2 -t paired -len 48 -gtf ./hg.GRCh37.75.gtf -c 0.00001 -novelSS 0 -analysis U -o $3



