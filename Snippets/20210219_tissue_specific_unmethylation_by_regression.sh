cd /mnt/isilon/zhoulab/labprojects/20210204_methSignature

cat <<'EOF' | tee 20210221_human_bloodRegression.pbs | qsub4
excel2tsv.R ~/samplesheets/2020/20200917_BLUEPRINTsubset_blood_references.xlsx |  awk 'NR>1&&$3=="N"{print $2;}' | tbmate view -l - -cd |  Rscript /mnt/isilon/zhoulab/labtools/Rutils/wzsearch_methsignatureLM.R <(excel2tsv.R ~/samplesheets/2020/20200917_BLUEPRINTsubset_blood_references.xlsx | awk 'NR>1&&$3=="N"{print $1,$5;}') - | gzip -c >/mnt/isilon/zhoulab/labprojects/20210204_methSignature/20210221_human_bloodRegression.tsv.gz
EOF

cat <<'EOF' | tee 20210220_human_tissueRegression.tsv_MLHL-B.pbs | qsub4
zcat /mnt/isilon/zhoulab/labprojects/20210204_methSignature/20210220_human_tissueRegression.tsv.gz | python /mnt/isilon/zhoulab/labtools/pyutils/wzsearch_methsignature2.py -u BCell,Plasma -x Testis,Placenta,Oocyte --fm 0.5 -m - | awk 'NR>1' | bedtools merge -i - | awk '{print $1,$2,$3,$3-$2}' | bedtools intersect -a - -b ~/references/hg19/annotation/cpg/cpg.bed -sorted -c >/mnt/isilon/zhoulab/labprojects/20210204_methSignature/20210220_human_tissueRegression.tsv_MLHL-B.bed
EOF

