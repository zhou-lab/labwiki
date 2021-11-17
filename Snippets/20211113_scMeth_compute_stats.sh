compute_stats() {
  # usage compute_stats allc samples.tsv
  # samples.tsv contains one column of sample names
  # this script update samples.tsv, and add one more column

  dir=$1
  meta=$2
  output=$3
  
  mkdir -p tmp
  rm -f tmp/*
  parallel -j 20 'f={/}; echo -e ${f%.bed.gz}"\t"$(zcat {} | wc -l) >tmp/${f%.bed.gz};' ::: $dir/*.bed.gz
  cat tmp/* >tmp_count

  awk 'NR==FNR{a[$1]=$2;}NR!=FNR{if($1 in a) aa=a[$1]; else aa="NA"; print $0,aa;}' tmp_count $meta >tmp_count2
  mv tmp_count2 $output
  rm -f tmp_count
  rm -f tmp/*
}
