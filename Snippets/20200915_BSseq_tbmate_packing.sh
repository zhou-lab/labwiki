cd /mnt/isilon/zhou_lab/projects/20200106_human_WGBS

mkdir -p hg19_tbk

for f in bed/*.bed.gz; do echo -ne $f"\t"; zcat $f | awk '{print NF;}' | head -1; done |  awk '$2==5' | cut -f1 >5col
for f in bed/*.bed.gz; do echo -ne $f"\t"; zcat $f | awk '{print NF;}' | head -1; done |  awk '$2==4' | cut -f1 >4col

doawk() { awk '{if($5<0)$5=0; print;}' $1; }
export -f doawk

## five column samples packed to float.int
cat 5col | parallel -j 20 'f={}; b=$(basename ${f%.bed.gz}); bedtools intersect -a /mnt/isilon/zhou_lab/projects/20191221_references/hg19/annotation/cpg/idx.gz -b {} -sorted -loj | cut -f1-3,8,9 | doawk | tbmate pack -s float.int -m /mnt/isilon/zhou_lab/projects/20191221_references/hg19/annotation/cpg/idx.gz - hg19_tbk/$b.tbk'

## four column samples packed to float
cat 4col | parallel -j 20 'f={}; b=$(basename ${f%.bed.gz}); bedtools intersect -a /mnt/isilon/zhou_lab/projects/20191221_references/hg19/annotation/cpg/idx.gz -b {} -sorted -loj | cut -f1-3,8 | doawk | tbmate pack -s float -m /mnt/isilon/zhou_lab/projects/20191221_references/hg19/annotation/cpg/idx.gz - hg19_tbk/$b.tbk'
