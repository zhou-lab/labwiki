#!/bin/sh

zcat /mnt/isilon/zhou_lab/HFS10T/2019_08_28_HPC_Laird_secondary/2016_05_13_InfiniumAnnotation/20180909/EPIC/EPIC.hg19.manifest.tsv.gz | awk -F"\t" -v OFS="\t" '$47=="FALSE" && $5~/^cg/{print $1,$2,$3,$4,$5,$26,$32,$38;}' | sortbed | bedtools intersect -a - -b ~/references/hg19/annotation/cpg/cpg_noDecoy.bed.gz -sorted -wo | awk '$12==2' | gzip -c >~/references/InfiniumArray/EPIC/EPIC.hg19.manifest.cg.bed.gz

bedtools intersect -a ~/references/hg19/annotation/cpg/cpg.bed -b ~/references/InfiniumArray/EPIC/EPIC.hg19.manifest.cg.bed.gz -sorted -loj | awk -v OFS="\t" 'NR==FNR{a[$1]=$4;}NR!=FNR{if($8 in a) aa=a[$8]; else aa=-1; print $1,$2,$3,aa,$8;}' <(zcat /home/zhouw3/references/InfiniumArray/EPIC/EPIC.idx.gz) - | bgzip -c >/home/zhouw3/references/InfiniumArray/EPIC/EPIC_to_hg19.idx.gz

tabix -p bed ~/references/InfiniumArray/EPIC/EPIC_to_hg19.idx.gz
