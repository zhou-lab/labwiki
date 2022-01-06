

function format_mapping_InfiniumI {
  mkdir -p tmp
  species=$1
  
  [[ -z "$2" ]] && fa=/mnt/isilon/zhou_lab/projects/20191221_references/Ensembl101/fasta/ || fa="$2"
  
  ## reformat A
  grep '^cg' sam/${species}.sam | awk '$2<256' | awk -F" " -v OFS="\t" '$1~/_1$/ && !and($2,0x4) && !and($2,0x10){print $3,$4+47,$4+49,substr($10,50,1),$0}$1~/_1$/ && !and($2,0x4) && and($2,0x10){print $3,$4-1,$4+1,substr($10,1,1),$0;} $1~/_1$/ && and($2,0x4){print "*",0,0,substr($10,50,1),$0;}' | awk -v OFS="\t" '!and($6,0x4) && ($5~/_1$/) && (($4=="A" && !and($6,0x10)) || ($4=="T" && and($6,0x10))){$5=substr($5,1,length($5)-2);print $0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >tmp/${species}_A
    
  ## reformat B
  grep '^cg' sam/${species}.sam | awk '$2<256' | awk -F" " -v OFS="\t" '$1~/_2$/ && !and($2,0x4) && !and($2,0x10){print $3,$4+47,$4+49,substr($10,50,1),$0}$1~/_2$/ && !and($2,0x4) && and($2,0x10){print $3,$4-1,$4+1,substr($10,1,1),$0;} $1~/_2$/ && and($2,0x4){print "*",0,0,substr($10,50,1),$0;}' | awk -v OFS="\t" '!and($6,0x4) && ($5~/_2$/) && (($4=="G" && !and($6,0x10)) || ($4=="C" && and($6,0x10))){$5=substr($5,1,length($5)-2);print $0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >tmp/${species}_B
  
  ## merge A,B
  awk 'NR==FNR{a[$5]=$0;}NR!=FNR && ($5 in a){$6=and($6,0x10); print $0"\t"a[$5];}' tmp/${species}_A tmp/${species}_B >tmp/${species}_AB
  ## add extension and target sequence
  awk '$1!="*" && $1==$15 && $2==$16 && $3==$17 && $10=="50M" && $24=="50M"' tmp/${species}_AB | awk '$6==0{print $1,$3,$3+1,$0;}$6==16{print $1,$2-1,$2,$0;}' | wzseqtk.py getfasta -i - -f ${fa}/${species}.fa | cut -f4- | awk '$6==0{print $1,$3-2,$3,$0;}$6==16{print $1,$2,$2+2,$0;}' | wzseqtk.py getfasta -i - -f /mnt/isilon/zhou_lab/projects/20191221_references/Ensembl101/fasta/${species}.fa | cut -f4- >tmp/${species}_AB_clean
  ## add color channel
  awk 'NR==FNR{ext[$5]=$29;tgt[$5]=$30;}NR!=FNR && ($5 in ext){print $0,ext[$5],tgt[$5];}' tmp/${species}_AB_clean tmp/${species}_AB | awk -f wanding.awk -e '{if($6==16) $29=dnarev($29); if($29=="C") col="G"; else if ($29==".") col="."; else col="R"; print $0"\t"col;}' | sortbed >tmp/${species}_AB_final
}

function format_mapping_InfiniumII {
  mkdir -p tmp
  species=$1
  ## reformat A
  grep '^cg' sam/${species}.sam | awk '$2<256' | awk '!and($2,0x4){$2=and($2,0x10);print $0;}' | awk -F" " -v OFS="\t" '!($1~/_1$/ || $1~/_2$/) && $2==0{print $3,$4+48,$4+50,substr($10,50,1),$0}!($1~/_1$/ || $1~/_2$/) && $2==16{print $3,$4-2,$4,substr($10,1,1),$0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >tmp/${species}_II
  
  ## add extension and target sequence
  awk '$1!="*" && $10=="50M"' tmp/${species}_II | wzseqtk.py getfasta -i - -f /mnt/isilon/zhou_lab/projects/20191221_references/Ensembl101/fasta/${species}.fa >tmp/${species}_II_clean
  ## add tgt
  awk 'NR==FNR{tgt[$5]=$15;}NR!=FNR && ($5 in tgt){print $0,tgt[$5];}' tmp/${species}_II_clean tmp/${species}_II | sortbed >tmp/${species}_II_final
}


