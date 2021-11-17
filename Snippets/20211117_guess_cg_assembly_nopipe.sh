guess_cg_assembly_nopipe() {
  sample=$1
  echo -ne $sample"\t"$(zcat $sample | awk '1;$2>10000000{exit;}' | wc -l)
  echo -ne "\thg19\t"$(zcat ~/references/hg19/annotation/cpg/cpg.bed.gz | awk '1;$2>10000000{exit;}' | bedtools intersect -a - -b <(zcat $sample | awk '1;$2>10000000{exit;}') -sorted | wc -l)
  echo -ne "\thg38\t"$(zcat ~/references/hg38/annotation/cpg/cpg.bed.gz | awk '1;$2>10000000{exit;}' | bedtools intersect -a - -b <(zcat $sample | awk '1;$2>10000000{exit;}') -sorted | wc -l)
  echo -ne "\tmm10\t"$(zcat ~/references/mm10/annotation/cpg/cpg.bed.gz | awk '1;$2>10000000{exit;}' | bedtools intersect -a - -b <(zcat $sample | awk '1;$2>10000000{exit;}') -sorted | wc -l)
  echo
}
