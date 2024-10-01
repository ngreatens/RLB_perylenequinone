#!/bin/bash

module load blast+ seqtk
gffpath=/project/fdwsru_fungal/Nick/git_repos/gffread/
PATH="$PATH:$gffpath" #add to path

##### 
# using input of gffs and fastas for target organisms, and a text file of the the target organisms and names of genes.
# ensure that gene coordinates are short and concise. Rename fastas and gffs if necessary



while read line; do
    genome=$(echo $line | awk '{print $1}')
    gene=$(echo $line | awk '{print $2}')
    echo -en "$genome\t"
    grep $gene ${genome}.gff | head -1 | cut -f 1,4 | tr "\n" "\t"
    grep $gene ${genome}.gff | tail -1 | cut -f 5,7
done < targetgenes.txt > targetgenes_boundaries.txt


#### get coordinates for cluster. the central genes (e.g. pks) +/- some amount
## note that the output might be negative or extend beyond the length of the contig/scaffold/chromosome. check this.

space=35000 #space on side of cluster

while read line; do
    genome=$(echo $line | awk '{print $1}')
    chrom=$(echo $line | awk '{print $2}')
    pks_start=$(echo $line | awk '{print $3}')
    pks_end=$(echo $line | awk '{print $4}')
    arr=("$pks_start" "$pks_end")
    pks_start_reordered=$(sort -n <(printf "%s\n" "${arr[@]}") | sed -n "1p")
    pks_end_reordered=$(sort -n <(printf "%s\n" "${arr[@]}") | sed -n "2p")
    pks_start_minus20=$(((${pks_start_reordered} - $space)))
    pks_end_plus20=$(((${pks_end_reordered} + $space)))
    echo -e "${genome}\t${chrom}:${pks_start_minus20}-${pks_end_plus20}"
done < targetgenes_boundaries.txt > cluster_coords.txt

#now with cluster coords available (and not negative),
# use gff read to subset gffs and get protein fastas for genes fully within the coordinates

gffpath=/project/fdwsru_fungal/Nick/git_repos/gffread/
PATH="$PATH:$gffpath" #add to path

while read line; do
    id=$(echo $line | awk '{print $1}')
    coords=$(echo $line | awk '{print $2}')
    gffread \
        -g ${id}.fna \
        -y ${id}_cluster_proteins.fa \
        -R \
        -r ${coords} \
        ${id}.gff
    gffread \
        -g ${id}.fna \
        -R \
        -r ${coords} \
        ${id}.gff | grep -P "mRNA\t" > ${id}_transcripts.gff
done < cluster_coords.txt

# remove comments on protein fasta if present and rename them with organism tag
rm *renamed_cluster_proteins.fa     #first make sure file is empty
for file in *cluster_proteins.fa; do
    cat $file | bioawk -c fastx -v id="$(basename $file _cluster_proteins.fa)" '{print ">"id"_"$name"\n"$seq}' > $(basename $file _cluster_proteins.fa)_renamed_cluster_proteins.fa
done 

# get all proteins in single fa file
cat *renamed_cluster_proteins.fa > prots.fa


# now get list of genes in tsv format that is able ot be parsed by gggenes
# note that it will still require some amount of configuring in r to orient genes and align a gene of choice

echo -e "molecule\tstart\tend\tstrand\tprot\tattributes" > genes.txt
for file in *transcripts.gff; do 
    while read line; do
        id=$(basename $file _transcripts.gff)
        start_pos=$(echo $line | awk '{print $4}')
        end_pos=$(echo $line | awk '{print $5}')
        nameprot=${id}_$(echo $line | awk '{print $9}' | cut -f1 -d";" | sed 's/ID=//1')
        attributes=$(echo $line | awk '{print $9}')
        strand=$(echo $line | awk '{print $7}')
        if [ "$strand" = "+" ] ; then 
            strand="1" 
        else 
            strand="0" 
        fi
        echo -e "$id\t$start_pos\t$end_pos\t$strand\t$nameprot\t$attributes"
    done < $file >> genes.txt
done

#now get blast hits
cluster=coniothyrium_glycines_cluster.fa #e.g. in this case all genes fully contained within the boundaries of the PKS +/- 35kb

module load blast+

makeblastdb -in $cluster -dbtype prot -out ${cluster%.*}

blastp \
    -db ${cluster%.*} \
    -query prots.fa \
    -outfmt 6  \
    -evalue 1e-10 \
    -qcov_hsp_perc 30 |
awk '$3 > 30 {print ;}' > blast.out

#get top match for each hit

# first get number of hits for each query seq. 
# if multiple hits are expected, this could be left out.
while read line; do 
    id=$(echo $line | awk '{print $1}')
    num=$(grep -c $id blast.out )
    echo -e $id $num; 
done < blast.out | uniq > hits_count.txt

# now for hits with just one hit, print blast output
# and for hits with more than one hit, sort for top hit by percent match and print only this blast output

while read line; do 
    id=$(echo $line | awk '{print $1}')
    count=$(echo $line | awk '{print $2}')
    if [ $count == 1 ]; then 
        grep "$id" blast.out
    else
         grep "$id" blast.out | sort -nrk3 | head -1 
    fi
done < hits_count.txt > blast_tophits.out

# now replace names of predicted proteins in genes file with those from the top blast hit in the cluster

head -1 genes.txt > genes2.txt
while read line; do
    prot=$(echo $line | awk '{print $5}')
    match=$(grep $prot blast_tophits.out | cut -f2 | tr -d "\n")
    if [ $(echo $match| wc -c) != 1 ]; then 
        echo -e "$(echo $line | awk '{print $1"\t"$2"\t"$3"\t"$4}')\t$match\t$(echo $line | awk '{print $6}')"
    fi
done < genes.txt >> genes2.txt

###############################################################################################################################################

# some amount of manual curation may be required before the next series of steps
# for example, if another gene nearby does not appear to be in the cluster, but has a blast hit to the query cluster, it may mes up the figure.
# delete lines as needed and record

#in this case, DELETE second Viridothelium hit for nearby PKS gene

### at this point, the genes2.txt file should be good to show similarity to the cluster.
### additional work and some amount of manual creation is required to select additional genes that will also appear in the cluster

# now select genes that are notable in other clusters and make fasta.

# get attributes of genes that did not have hits cluster

cat genes.txt | awk '{print $6}'| sort > genes_attrributes.txt
cat genes2.txt | awk '{print $6}' |sort > genes2_attributes.txt
comm -3 genes_attrributes.txt genes2_attributes.txt > atts_without_hits.txt

#now use attributes to get proteins and subset these proteins from earlier protein file
while read line; do
    grep $line genes.txt | awk '{print $5}'
done < atts_without_hits.txt > prots_without_hits.txt

module load seqtk

seqtk subseq prots.fa prots_without_hits.txt > prots_without_hits.fa

#### now blast these hits against other refs

other_refs=other_refs.fa

module load blast+

makeblastdb -in $other_refs -dbtype prot -out ${other_refs%.*}

blastp \
    -db ${other_refs%.*} \
    -query prots_without_hits.fa \
    -outfmt 6  \
    -evalue 1e-10 \
    -qcov_hsp_perc 30 |
awk '$3 > 30 {print ;}' > blast_against_other_refs.out

while read line; do
    prot=$(echo $line | awk '{print $5}')
    match=$(grep $prot blast_against_other_refs.out | cut -f2 | tr -d "\n")
    if [ $(echo $match| wc -c) != 1 ]; then 
        echo -e "$(echo $line | awk '{print $1"\t"$2"\t"$3"\t"$4}')\t$match\t$(echo $line | awk '{print $6}')"
    fi
done < genes.txt >> genes2.txt

#add other genes that are within the boundaries of the comparison clusters and check for similarity with other genes 

# again get coordinates for cluster 

while read line; do
    genome=$(echo $line | awk '{print $1}')
    gene=$(echo $line | awk '{print $2}')
    echo -en "$genome\t"
    grep $gene ${genome}.gff | head -1 | cut -f 1
done < targetgenes.txt | sort > cluster_ids_chroms.txt

while read line; do
    id=$(echo $line | awk '{print $1}')
    chrom=$(echo $line | awk '{print $2}')
    start=$(grep $id genes2.txt | sort -n -k2 | head -1 | awk '{print $2}')
    end=$(grep $id genes2.txt | sort -n -k2 | tail -1 | awk '{print $3}')
    echo -e "${id}\t${chrom}:${start}-${end}"
done < cluster_ids_chroms.txt > query_hits_coords.txt

## again subset gffs to pull out genes within clusters, more narrowly

while read line; do
    id=$(echo $line | awk '{print $1}')
    coords=$(echo $line | awk '{print $2}')
    gffread \
        -g ${id}.fna \
        -y ${id}_cluster_proteins_strict.fa \
        -R \
        -r ${coords} \
        ${id}.gff
    gffread \
        -g ${id}.fna \
        -R \
        -r ${coords} \
        ${id}.gff | grep -P "mRNA\t" > ${id}_transcripts_strict.gff
done < query_hits_coords.txt

# remove comments on protein fasta if present and rename them with organism tag
rm *renamed_cluster_proteins_strict.fa     #first make sure file is empty
for file in *cluster_proteins_strict.fa; do
    cat $file | bioawk -c fastx -v id="$(basename $file _cluster_proteins_strict.fa)" '{print ">"id"_"$name"\n"$seq}' > $(basename $file _cluster_proteins_strict.fa)_renamed_cluster_proteins_strict.fa
done 

# get all proteins in single fa file again
cat *renamed_cluster_proteins_strict.fa > prots_strict.fa

# now compare names with blast out file and subset prots without match

cat blast.out blast_against_other_refs.out | awk '{print $1}' | sort -f - | uniq > prots_with_hits.txt
cat prots_strict.fa | bioawk -c fastx '{print $name}' | sort > prots_strict.txt
comm -3 prots_strict.txt prots_with_hits.txt | grep -vP "\t" > other_prots_in_cluster.txt


##get proteins 
seqtk subseq prots.fa other_prots_in_cluster.txt > other_prots_in_cluster.fa

##blast proteins against each other
other_genes=other_prots_in_cluster.fa

ml blast+
makeblastdb -in $other_genes -dbtype prot -out ${other_genes%.*}


blastp \
    -db ${other_genes%.*} \
    -query other_prots_in_cluster.fa\
    -outfmt 6  \
    -evalue 1e-10 \
    -qcov_hsp_perc 30 |
awk '$3 > 30 {print ;}' > blast_other_genes_against_other_genes.out

#Manually inspect blast hits
while read line; do 
    query=$(echo $line | awk '{print $1}')
    hit=$(echo $line | awk '{print $2}')
    if [ $query != $hit ]; then 
        echo $line
    fi
done < blast_other_genes_against_other_genes.out > blast_other_genes_against_other_genes_unique.out 

#Manually add hits to genes file if needed. 

for file in *transcripts_strict.gff; do 
    while read line; do
        id=$(basename $file _transcripts_strict.gff)
        start_pos=$(echo $line | awk '{print $4}')
        end_pos=$(echo $line | awk '{print $5}')
        nameprot=${id}_$(echo $line | awk '{print $9}' | cut -f1 -d";" | sed 's/ID=//1')
        attributes=$(echo $line | awk '{print $9}')
        strand=$(echo $line | awk '{print $7}')
        if [ "$strand" = "+" ] ; then 
            strand="1" 
        else 
            strand="0" 
        fi
        echo -e "$id\t$start_pos\t$end_pos\t$strand\tother\t$attributes"
    done < $file 
done > tmp

while read line; do
    id=$(echo $line | awk '{print $1}')
    sed -i "s/${id}_//g" other_prots_in_cluster.txt
done < targetgenes.txt

while read line; do
    grep $line tmp
done < other_prots_in_cluster.txt > addl_genes.txt

rm tmp

# add query to figure. Here just down manually

gffread -r scaff_46:93700-120242 Coniothyrium_glycines.gff | grep -P "mRNA\t" > Coniothyrium_glycines_transcripts.gff

for file in *transcripts.gff; do 
    while read line; do
        id=$(basename $file _transcripts.gff)
        start_pos=$(echo $line | awk '{print $4}')
        end_pos=$(echo $line | awk '{print $5}')
        nameprot=$(echo $line | awk '{print $9}' | cut -f1 -d";" | sed 's/ID=//1') 
        attributes=$(echo $line | awk '{print $9}')
        strand=$(echo $line | awk '{print $7}')
        if [ "$strand" = "+" ] ; then 
            strand="1" 
        else 
            strand="0" 
        fi
        echo -e "$id\t$start_pos\t$end_pos\t$strand\t$nameprot\t$attributes"
    done < $file
done > Coniothyrium_glycines_genes.txt 

#custom to eliminate a few clusters, etc. 
cat genes2.txt addl_genes.txt Coniothyrium_glycines_genes.txt  > genes3.txt
sed -i 's/ctb6/other/1' genes3.txt

#here just renaming the gene for hte genes file to make them appear first in the figure legend
sed -i 's/g653/a653/1' genes3.txt

# get final cluster coords for figure

while read line; do
    id=$(echo $line | awk '{print $1}')
    chrom=$(echo $line | awk '{print $2}')
    start=$(grep $id genes2.txt | sort -n -k2 | head -1 | awk '{print $2}')
    end=$(grep $id genes2.txt | sort -n -k2 | tail -1 | awk '{print $3}')
    echo -e "${id}\t${chrom}:${start}-${end}"
done < cluster_ids_chroms.txt >cluster_final_coords.txt




