

function format_mapping_InfiniumI {
  mkdir -p tmp
  species=$1
  [[ -z "$2" ]] && fa=/mnt/isilon/zhou_lab/projects/20191221_references/Ensembl101/fasta/ || fa="$2"
  
  ## reformat A
  grep '^cg' sam/${species}.sam | awk '$2<256 && $1~/_1$/' | awk '{$2=and($2,0x14);print $0;}' | awk -F" " -v OFS="\t" '!and($2,0x4) && !and($2,0x10){print $3,$4+47,$4+49,substr($10,50,1),$0} !and($2,0x4) && and($2,0x10){print $3,$4-1,$4+1,substr($10,1,1),$0;} and($2,0x4){print "*",0,0,substr($10,50,1),$0;}' | awk -v OFS="\t" '{$5=substr($5,1,length($5)-2);print $0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >tmp/${species}_A
    
  ## reformat B
  grep '^cg' sam/${species}.sam | awk '$2<256 && $1~/_2$/' | awk '{$2=and($2,0x14);print $0;}'| awk -F" " -v OFS="\t" '!and($2,0x4) && !and($2,0x10){print $3,$4+47,$4+49,substr($10,50,1),$0} !and($2,0x4) && and($2,0x10){print $3,$4-1,$4+1,substr($10,1,1),$0;} and($2,0x4){print "*",0,0,substr($10,50,1),$0;}' | awk -v OFS="\t" '{$5=substr($5,1,length($5)-2);print $0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >tmp/${species}_B
  
  ## merge A,B
  awk 'NR==FNR{a[$5]=$0;}NR!=FNR && ($5 in a){$6=and($6,0x14); print $0"\t"a[$5];}' tmp/${species}_B tmp/${species}_A >tmp/${species}_AB
  
  ## add extension and target sequence
  awk '$1!="*" && $1==$15 && $2==$16 && $3==$17 && $10=="50M" && $24=="50M"' tmp/${species}_AB | awk '$6==0{print $1,$3,$3+1,$0;}$6==16{print $1,$2-1,$2,$0;}' | wzseqtk.py getfasta -i - -f ${fa}/${species}.fa | cut -f4- | awk '$6==0{print $1,$3-2,$3,$0;}$6==16{print $1,$2,$2+2,$0;}' | wzseqtk.py getfasta -i - -f ${fa}/${species}.fa | cut -f4- >tmp/${species}_AB_clean
  ## add color channel
  awk 'NR==FNR{ext[$5]=$29;tgt[$5]=$30;}NR!=FNR{if($5 in ext) {extn=ext[$5]; tgtn=tgt[$5];} else {extn="NA"; tgtn="NA";} print $0,extn,tgtn;}' tmp/${species}_AB_clean tmp/${species}_AB | awk -f wanding.awk -e '{if($6==16) $29=dnarev($29); if($29=="C") col="G"; else if ($29=="NA") col="NA"; else col="R"; print $0"\t"col;}' | sortbed >tmp/${species}_AB_final
}

function format_mapping_InfiniumII {
  mkdir -p tmp
  species=$1
  [[ -z "$2" ]] && fa=/mnt/isilon/zhou_lab/projects/20191221_references/Ensembl101/fasta/ || fa="$2"
  
  ## reformat A
  grep '^cg' sam/${species}.sam | awk '$2<256' | awk '{$2=and($2,0x14);print $0;}' | awk -F" " -v OFS="\t" '!($1~/_1$/ || $1~/_2$/) && $2==0{if($1~/_[TBN]O/){print $3,$4+49,$4+51,substr($10,50,1),$0;}else{print $3,$4+48,$4+50,substr($10,50,1),$0;}}!($1~/_1$/ || $1~/_2$/) && $2==16{if($1~/_[TBN]O/){print $3,$4-3,$4-1,substr($10,1,1),$0;}else{print $3,$4-2,$4,substr($10,1,1),$0;}}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >tmp/${species}_II
  grep '^ch' sam/${species}.sam | awk '$2<256' | awk '{$2=and($2,0x14);print $0;}' | awk -F" " -v OFS="\t" '!($1~/_1$/ || $1~/_2$/) && $2==0{print $3,$4+49,$4+50,substr($10,50,1),$0}!($1~/_1$/ || $1~/_2$/) && $2==16{print $3,$4-2,$4-1,substr($10,1,1),$0;} !($1~/_1$/ || $1~/_2$/) && $2==4{print "*",0,0,substr($10,50,1),$0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >>tmp/${species}_II
  
  ## add extension and target sequence
  awk '$1!="*" && $10=="50M"' tmp/${species}_II | wzseqtk.py getfasta -i - -f ${fa}/${species}.fa >tmp/${species}_II_clean

  ## add tgt
  awk 'NR==FNR{tgt[$5]=$15;}NR!=FNR && ($5 in tgt){print $0,tgt[$5];}' tmp/${species}_II_clean tmp/${species}_II | sortbed >tmp/${species}_II_final
}

function validate_Infinium {
  mkdir -p validation
  species=$1
  [[ -z "$2" ]] && fa=/mnt/isilon/zhou_lab/projects/20191221_references/Ensembl101/fasta/ || fa="$2"
  
  awk '$1!="*" && ($10=="50M" || $24=="50M"){$2=$2-50;$3=$3+50;print $0;}' tmp/${species}_AB | wzseqtk.py getfasta -i - -f ${fa}/${species}.fa | awk -f wanding.awk -e 'NR!=FNR{print $5, "I", s1[$5]; print joinr(1,28); print $29; if(!and($6,0x10)){printf("> ");} else {printf("<"); for(i=1;i<50;++i) printf(" "); } print $11; if(!and($6,0x10)) {printf("> "); os=s1[$5];}else {printf("<"); for(i=1;i<50;++i) printf(" "); os=dnarev(s1[$5])} print os; if(!and($6,0x10)){printf("> ");} else {printf("<"); for(i=1;i<50;++i) printf(" "); } print $25; if(!and($6,0x10)) {printf("> "); os=s2[$5];}else {printf("<"); for(i=1;i<50;++i) printf(" "); os=dnarev(s2[$5])} print os;print("\n");}NR==FNR{s1[$1]=$2;s2[$1]=$3;}' tmp/probe2originalseq.txt - >validation/${species}.txt
  awk '$1!="*" && $10=="50M"{$2=$2-50;$3=$3+50;print $0;}' tmp/${species}_II | wzseqtk.py getfasta -i - -f ${fa}/${species}.fa | awk -f wanding.awk -e 'NR!=FNR{print $5, "II", s1[$5]; print joinr(1,14); print $15; if(!and($6,0x10)){printf(">");} else {printf("<"); for(i=1;i<51;++i) printf(" "); } print $11; if(!and($6,0x10)) {printf(">"); os=s1[$5];}else {printf("<"); for(i=1;i<51;++i) printf(" "); os=dnarev(s1[$5])} print os; print("\n");}NR==FNR{s1[$1]=$2}' tmp/probe2originalseq.txt - >>validation/${species}.txt
}

