#!/bin/bash

# submit this script to Eddie

#$ -cwd
#$ -l h_rt=48:00:0,h_vmem=24G #we are asking for a 48h time slot and 24GB RAM
#$ -pe sharedmem 8

# get an email when job is submitted, ended or removed
#$ -M UUN@ed.ac.uk
#$ -m bes

# allow for modules to be loaded (we need this for the packages we're using)
source /etc/profile.d/modules.sh

# Load modules, i.e. muscle and iqtree
module load roslin/muscle/5.3
module load roslin/iqtree/2.4.0

####################

# run the actual code

# give the name of the file you're using
filename="DC11_filtered_sequences.noATS" #remove .fasta ending!

# align sequences
muscle -align $filename.fasta -output $filename.aligned.fasta

# construct tree
iqtree2 --safe -s $filename.aligned.fasta -m JTT+F+I+I+R10 -alrt 1000 -nt 8
