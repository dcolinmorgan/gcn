#!/bin/bash
#PBS -l nodes=1:ppn=48
#PBS -l mem=200g
#PBS -l walltime=96:00:00
#PBS -N gcn
#PBS -q cgsd
#PBS -e .GCerr
#PBS -o .GCout
##qsub run/gcn/gcn.sh

module purge
module add miniconda3
module add parallel
source activate mypy3
# https://stackoverflow.com/questions/54245735/unable-to-run-pbs-script-on-multiple-nodes-using-gnu-parallel
# sort -u $PBS_NODEFILE > nodelist.dat
# export JOBS_PER_NODE=1

chmod +x run/gcn/gcn.py

#################PRODIGAL in parallel##############
function gcnXpat {
patt=$1
pat=$(echo "$1"|cut -d/ -f8)
eval "prodigal -i "$patt" -o run/gcn/pat/"$pat"_out -f gff -d run/gcn/pat/"$pat"_genes -a run/gcn/pat/"$pat"_proteins -p meta"

# eval "diamond blastp -b12 -c1 -d ../../groups/cgsd/gordonq/database/gunc_db_gtdb95.dmnd -q run/gcn/pat/"$pat"_proteins --threads 64 -o run/gcn/pat/"$pat"_proteins_o -f6"#--multiprocessing --tmpdir run/gcn/tmpdir --parallel-tmpdir run/gcn/ptmpdir"
}

export -f gcnXpat
# parallel -j 1 gcnXpat ::: ../../groups/cgsd/gordonq/LauG_Metagenomics_CPOS-200710-CJX-3455a/assemble/*/final.contigs.fa

#################DIAMOND in multithreading##############
# patdir='run/gcn/pat/*'
# pats=$(ls $patdir*_proteins)

# for pat in $pats
# do
# # eval "diamond makedb "
# eval "diamond blastp -b12 -c1 -d ../../groups/cgsd/gordonq/database/uniref90/uniref90_201901.dmnd -q "$pat" --threads 64 -o "$pat"_o3 -f6" # qseqid sseqid evalue bitscore staxids sscinames" #--multiprocessing --tmpdir run/gcn/tmpdir --parallel-tmpdir run/gcn/ptmpdir"
# done

#################make network in parallel##############
function dia2nk {

patt=$1
pat=$(echo "$1"|cut -d/ -f4|cut -d_ -f1)

awk '($3>95)' $patt > tmp/"$pat"_nko
cut -f2 tmp/"$pat"_nko |cut -d '|' -f 1|cut -d '_' -f 2 > tmp/"$pat"_tmpA
cut -f2 tmp/"$pat"_nko |cut -d '|' -f 2 > tmp/"$pat"_tmpB
eval "pr -mt -s, tmp/"$pat"_tmpB tmp/"$pat"_tmpA >run/gcn/pat/"$pat"_nk"
rm tmp/*tmp*

# python run/gcn/gcn.py $pat

}

export -f dia2nk

# parallel dia2nk ::: run/gcn/pat/*_proteins_o2

########################################################################
function cnvt {

patt=$1
pat=$(echo "$1"|cut -d/ -f4|cut -d_ -f1)
less "$patt" | tr , '\t' > "$patt"_t
}

export -f cnvt

# parallel cnvt ::: run/gcn/pat/*_nk

########################################################################
function humann3X {
pat=$1
humann --threads 40 --input-format fasta --input $pat --output humann3 --nucleotide-database /groups/cgsd/gordonq/database/humann_chocophlan/ --protein-database /groups/cgsd/gordonq/database/uniref90/ --remove-temp-output --metaphlan-options "--bowtie2db /home/dcmorgan/.conda/envs/mypy3/lib/python3.9/site-packages/metaphlan/metaphlan_databases/"
# https://github.com/biobakery/biobakery/wiki/humann3#22-running-humann-the-basics
}

export -f humann3X
# parallel humann3X ::: /groups/cgsd/gordonq/LauG_Metagenomics_CPOS-200710-CJX-3455a/primary_seq/merged/*.fasta


patdir='/groups/cgsd/gordonq/LauG_Metagenomics_CPOS-200710-CJX-3455a/primary_seq/merged/*.fasta'
pats=$(ls $patdir)

# for pat in $pats
# do
# humann3X $pat
# done
# humann_join_tables -i /groups/cgsd/gordonq/all_hypertension/humann3_res/ -o ~/ht_subset_genefamilies.tsv --file_name genefamilies

# humann_renorm_table -i ~/ht_subset_genefamilies.tsv -o ~/data/gcn/ht_genefamilies-cpm.tsv --units cpm


########################################################################
python run/gcn/gcn.py