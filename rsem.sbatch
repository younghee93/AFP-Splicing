#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=64000
#SBATCH -o ./05.rsem_result/tmp/rsem.o.%j
#SBATCH --err ./05.rsem_result/tmp/rsem.e.%j

ml star/2.5
ml rsem/1.3.0
ml samtools

rsem-calculate-expression  --paired-end --star $1 $2 ./5.reference/5.rsem/hg19 $3


