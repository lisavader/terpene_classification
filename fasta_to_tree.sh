#!/bin/sh

fasta=$1  #specify path to fasta file as command line argument
file_name=$(basename $fasta .fasta)
dir_path=$(dirname $fasta)
data_path=${dir_path%/*}

muscle -in $fasta -out "${data_path}/msa/${file_name}.afa"
#java -jar $HOME/BMGE/src/BMGE.jar -t AA -i "${data_path}/msa/${file_name}.afa" -of "${data_path}/msa/${file_name}_trim.afa"
FastTree "${data_path}/msa/${file_name}.afa" > "${data_path}/newick/${file_name}.tree"