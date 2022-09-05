function test_export {
  ## ex1
  export BISCUIT_FASTA=~/references/hg19/biscuit/hg19.fa
  export BASEDIR=~/zhou_lab/projects/20220822_rebuild_Infinium_array_manifest
  export REFCODE=hg19
  export PLATFORM=EPIC+
  export CSV=~/zhou_lab/projects/20220822_rebuild_Infinium_array_manifest/csv/EPIC+/CombinedManifestEPIC.manifest.LEGXservices.csv.gz
  export TMPFDR=$BASEDIR/tmp/${PLATFORM}_${REFCODE}

  ## ex2
  export BISCUIT_FASTA=~/references/hg19/biscuit/hg19.fa
  export BASEDIR=~/zhou_lab/projects/20220822_rebuild_Infinium_array_manifest
  export REFCODE=hg19
  export PLATFORM=HM27
  export CSV=~/zhou_lab/projects/20220822_rebuild_Infinium_array_manifest/csv/HM27/GPL8490_HumanMethylation27_270596_v.1.2.csv.gz
  export TMPFDR=$BASEDIR/tmp/${PLATFORM}_${REFCODE}
  export parser=HM27

  ## ex3
  export BISCUIT_FASTA=~/references/hg38/biscuit/hg38.fa
  export BASEDIR=~/zhou_lab/projects/20220822_rebuild_Infinium_array_manifest
  export REFCODE=hg38
  export PLATFORM=Mammal40
  export CSV=~/zhou_lab/projects/20220822_rebuild_Infinium_array_manifest/csv/Mammal40/HorvathMammal40.CanonicalManifest.3.2019.manifest.csv
  export TMPFDR=$BASEDIR/tmp/${PLATFORM}_${REFCODE}
  export parser=Mammal40
}

function set_environment {
  export BISCUIT_FASTA=$1
  export BASEDIR=$2
  export REFCODE=$3             # hg38, mm10, ...
  export PLATFORM=$4            # how shall we call this new platform
  export CSV=$5                 # original csv file
  export parser=$6              # which parser to use
  export GMGTF=$7               # gene model gtf.gz
  export GMCODE=$8              # gene model name we should call the gene annotation file
  export TMPFDR=$BASEDIR/tmp/${PLATFORM}_${REFCODE}
}

function run_pipeline_build_Infinium_20220823() ( ## note that fun() (...) runs in a subshell
  set_environment $1 $2 $3 $4 $5 $6 $7 $8
  cd $BASEDIR
  rm -rf $TMPFDR

  runpipe1 | tee $BASEDIR/tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz.log
  echo "Temp folder:" $TMPFDR
  echo "Please delete this folder after checking"
)

function runpipe1 {
  csv2standard_input_tsv_$parser
  prepare_fa
  biscuit_mapping_thread5
  format_mapping_InfiniumI
  format_mapping_InfiniumII
  validate_Infinium
  ~/repo/labwiki/Snippets/20220105_Infinium_mapping_mergeV2.R $TMPFDR ref $BASEDIR/tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz NA
  count_manifest
  ~/repo/labwiki/Snippets/20220105_build_sesameAddress.R $BASEDIR/tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz RDS_ordering/${PLATFORM}_${REFCODE}.rds
  ~/repo/labwiki/Snippets/20220105_build_sesameAddressGR.R $BASEDIR/tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz RDS_GRanges/${PLATFORM}_${REFCODE}.rds decoy
  set_features_$REFCODE
  buildFeatureOverlaps
  buildFeatureGene
}

function csv2standard_input_tsv_HM27 {
  # requirement: csvtk
  # column requirement:
  # 1:cgNumber, 7:I/II, 3:AddressA_ID, 4:AlleleA_ProbeSeq, 5:AlleleB_ID, 6:AlelleB_ProbeSeq
  mkdir -p $TMPFDR/fa
  cd $TMPFDR/fa
  zcat -f $CSV | awk -F"," 'BEGIN{a=0}/^IlmnID/{a=1;next;}/^\[Controls\]/{a=0;}a{print $1,"I",$4,$5,$6,$7;}' >standard_input.tsv
  zcat -f $CSV | awk -F"," 'a{print "ctl_"$1,$1,"NA","NA","II";}/^\[Control/{a=1}' >standard_input_control.tsv
  zcat -f $CSV | awk -F"," 'a{print "ctl_"$1,$1,$2,$3,$4,"II";}/^\[Control/{a=1}' >standard_input_control_anno.tsv
  cd $BASEDIR
}

function csv2standard_input_tsv_EPIC {
  # requirement: csvtk
  # column requirement:
  # 1:cgNumber, 7:I/II, 3:AddressA_ID, 4:AlleleA_ProbeSeq, 5:AlleleB_ID, 6:AlelleB_ProbeSeq
  mkdir -p $TMPFDR/fa
  cd $TMPFDR/fa
  zcat -f $CSV | awk 'BEGIN{a=0}/^IlmnID/{a=1;next;}/^\[Controls\]/{a=0;}a' | csvtk csv2tab | awk '{print $1,$7,$3,$4,$5,$6}' >standard_input.tsv
  zcat -f $CSV | awk -F"," 'a{print "ctl_"$1,$1,"NA","NA","II";}/^\[Control/{a=1}' >standard_input_control.tsv
  zcat -f $CSV | awk -F"," 'a{print "ctl_"$1,$1,$2,$3,$4,"II";}/^\[Control/{a=1}' >standard_input_control_anno.tsv
  cd $BASEDIR
}

function csv2standard_input_tsv_Mammal40 {
  # requirement: csvtk
  # column requirement:
  # 1:cgNumber, 20:I/II, 3:AddressA_ID, 4:AlleleA_ProbeSeq, 6:AlleleB_ID, 7:AlelleB_ProbeSeq
  mkdir -p $TMPFDR/fa
  cd $TMPFDR/fa
  zcat -f $CSV | awk 'BEGIN{a=0}/^IlmnID/{a=1;next;}/^\[Controls\]/{a=0;}a' | csvtk csv2tab | awk '{print $1,$20,$3,$4,$6,$7}' >standard_input.tsv
  zcat -f $CSV | awk -F"," 'a{print "ctl_"$1,$1,"NA","NA","II";}/^\[Control/{a=1}' >standard_input_control.tsv
  zcat -f $CSV | awk -F"," 'a{print "ctl_"$1,$1,$2,$3,$4,"II";}/^\[Control/{a=1}' >standard_input_control_anno.tsv
  cd $BASEDIR
}

function csv2standard_input_tsv_EPICplus {
  # requirement: csvtk
  # column requirement:
  mkdir -p $TMPFDR/fa
  cd $TMPFDR/fa
  # 1:cgNumber, 2:I/II, 3:AddressA_ID, 4:AlleleA_ProbeSeq, 5:AlleleB_ID, 6:AlelleB_ProbeSeq
  zcat -f $CSV | awk 'BEGIN{a=0}/^IlmnID/{a=1;next;}/^[0-9]/{a=0;}a' | csvtk csv2tab | awk '{if ($6!="") $6=int($6); print $1,$9,$4,$5,$6,$7}' >standard_input.tsv
  zcat -f $CSV | awk -F"," 'a{print "ctl_"$1,$1,"NA","NA","II";}/^[0-9]/{a=1}' >standard_input_control.tsv
  zcat -f $CSV | awk -F"," 'a{print "ctl_"$1,$1,$2,$3,$4,"II";}/^[0-9]/{a=1}' >standard_input_control_anno.tsv
  cd $BASEDIR
}

function csv2standard_input_tsv_EPICv2draft {
  # requirement: csvtk
  # column requirement:
  mkdir -p $TMPFDR/fa
  cd $TMPFDR/fa
  # 1:cgNumber, 2:I/II, 3:AddressA_ID, 4:AlleleA_ProbeSeq, 5:AlleleB_ID, 6:AlelleB_ProbeSeq
  zcat -f $CSV | awk 'BEGIN{a=0}/^IlmnID/{a=1;next;}/^\[Controls\]/{a=0;}a' | csvtk csv2tab | awk '{print $1,$15,$3,$4,$5,$6}' >standard_input.tsv
  zcat -f $CSV | awk -F"," 'a{print "ctl_"$1,$1,"NA","NA","II";}/^\[Control/{a=1}' >standard_input_control.tsv
  zcat -f $CSV | awk -F"," 'a{print "ctl_"$1,$1,$2,$3,$4,"II";}/^\[Control/{a=1}' >standard_input_control_anno.tsv
  cd $BASEDIR
}

function csv2standard_input_tsv_MM285 {
  # requirement: csvtk
  # column requirement:
  mkdir -p $TMPFDR/fa
  cd $TMPFDR/fa
  # 1:cgNumber, 2:I/II inferred from AddressA_ID, 3:AddressA_ID, 4:AlleleA_ProbeSeq, 5:AlleleB_ID, 6:AlelleB_ProbeSeq
  zcat -f $CSV | awk 'BEGIN{a=0}/^IlmnID/{a=1;next;}/^\[Controls\]/{a=0;}a' | csvtk csv2tab | awk '{if($5=="") type="II"; else type="I"; print $1,type,$3,$4,$5,$6}' >standard_input.tsv
  zcat -f $CSV | awk -F"," 'a{print "ctl_"$1,$1,"NA","NA","II";}/^\[Control/{a=1}' >standard_input_control.tsv
  zcat -f $CSV | awk -F"," 'a{print "ctl_"$1,$1,$2,$3,$4,"II";}/^\[Control/{a=1}' >standard_input_control_anno.tsv
  cd $BASEDIR
}

function prepare_fa {
  cd $TMPFDR/fa
  echo "=== Prepare FASTA ==="
  cat standard_input.tsv | awk '$2=="II"{print ">"$1; print gensub(/[MWKSVBDH]/,"N","g",gensub(/R/,"A","g",gensub(/Y/,"T","g",$4)));}' >Infinium_II.fa
  cat standard_input.tsv | awk '$2=="I"{print ">"$1"_1" >"Infinium_I_1.fa"; print gensub(/[MWKSVBDH]/,"N","g",gensub(/R/,"A","g",gensub(/Y/,"T","g",$4))) >"Infinium_I_1.fa"; print ">"$1"_2" >"Infinium_I_2.fa"; print gensub(/[MWKSVBDH]/,"N","g",gensub(/R/,"A","g",gensub(/Y/,"T","g",$6))) >"Infinium_I_2.fa";}'
  cat standard_input.tsv | awk '{if($6=="") $6="NA"; print $1,$4,$6;}' >probe2originalseq.txt
  cat standard_input.tsv | awk '{if($5=="") $5="NA"; print $1,$3,$5;}' >probe2address.txt
  ## sanity check, no non-conventional bases
  grep -v '^>' Infinium_I_1.fa | grep '[^ATGCNMWKSVBDH]'
  grep -v '^>' Infinium_I_2.fa | grep '[^ATGCNMWKSVBDH]'
  grep -v '^>' Infinium_II.fa | grep '[^ATGCNMWKSVBDH]'
  cd $BASEDIR
}

function biscuit_mapping_thread5 {
  cd $TMPFDR
  echo "=== BISCUIT alignment ===="
  biscuit align -F@ 5 $BISCUIT_FASTA fa/Infinium_II.fa | grep -v "^@" >ref.sam;
  biscuit align -F@ 5 $BISCUIT_FASTA fa/Infinium_I_1.fa fa/Infinium_I_2.fa | grep -v "^@" >>ref.sam
  cd $BASEDIR
}

function format_mapping_InfiniumI {

  cd $TMPFDR
  # [[ -z "$2" ]] && fa=/mnt/isilon/zhou_lab/projects/20191221_references/Ensembl101/fasta/ || fa="$2"
  echo "=== Format Inf-I ===="
  ## reformat A
  grep '^cg' ref.sam | awk '$2<256 && $1~/_1$/' | awk '{$2=and($2,0x14);print $0;}' | awk -F" " -v OFS="\t" '!and($2,0x4) && !and($2,0x10){if($1~/_[TBN]O/){print $3,$4+48,$4+50,substr($10,50,1),$0;} else {print $3,$4+47,$4+49,substr($10,50,1),$0;}} !and($2,0x4) && and($2,0x10){if($1~/_[TBN]O/){print $3,$4-2,$4,substr($10,1,1),$0;}else{print $3,$4-1,$4+1,substr($10,1,1),$0;}} and($2,0x4){print "*",0,0,substr($10,50,1),$0;}' | awk -v OFS="\t" '{$5=substr($5,1,length($5)-2);print $0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >ref_A
  grep '^rs' ref.sam | awk '$2<256 && $1~/_1$/' | awk '{$2=and($2,0x14);print $0;}' | awk -F" " -v OFS="\t" '!and($2,0x4) && !and($2,0x10){print $3,$4+48,$4+49,substr($10,50,1),$0} !and($2,0x4) && and($2,0x10){print $3,$4-1,$4,substr($10,1,1),$0;} and($2,0x4){print "*",0,0,substr($10,50,1),$0;}' | awk -v OFS="\t" '{$5=substr($5,1,length($5)-2);print $0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >>ref_A
  ## reformat B
  grep '^cg' ref.sam | awk '$2<256 && $1~/_2$/' | awk '{$2=and($2,0x14);print $0;}' | awk -F" " -v OFS="\t" '!and($2,0x4) && !and($2,0x10){if($1~/_[TBN]O/){print $3,$4+48,$4+50,substr($10,50,1),$0;} else {print $3,$4+47,$4+49,substr($10,50,1),$0;}} !and($2,0x4) && and($2,0x10){if($1~/_[TBN]O/){print $3,$4-2,$4,substr($10,1,1),$0;}else{print $3,$4-1,$4+1,substr($10,1,1),$0;}} and($2,0x4){print "*",0,0,substr($10,50,1),$0;}' | awk -v OFS="\t" '{$5=substr($5,1,length($5)-2);print $0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >ref_B
  grep '^rs' ref.sam | awk '$2<256 && $1~/_2$/' | awk '{$2=and($2,0x14);print $0;}' | awk -F" " -v OFS="\t" '!and($2,0x4) && !and($2,0x10){print $3,$4+48,$4+49,substr($10,50,1),$0} !and($2,0x4) && and($2,0x10){print $3,$4-1,$4,substr($10,1,1),$0;} and($2,0x4){print "*",0,0,substr($10,50,1),$0;}' | awk -v OFS="\t" '{$5=substr($5,1,length($5)-2);print $0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >>ref_B
  ## merge A,B
  awk 'NR==FNR{a[$5]=$0;}NR!=FNR && ($5 in a){$6=and($6,0x14); print $0"\t"a[$5];}' ref_B ref_A >ref_AB
  
  ## add extension and target sequence
  awk '$1!="*" && $10=="50M"' ref_AB | awk '$6==0{if($5~/_[TBN]O/ && !($5~/^rs/)){pos=$3-1;}else{pos=$3;} print $1,pos,pos+1,$0;}$6==16{if($5~/_[TBN]O/ && !($5~/^rs/)){pos=$2+1;}else{pos=$2;} print $1,pos-1,pos,$0;}' | wzseqtk.py getfasta -i - -f $BISCUIT_FASTA | cut -f4- | awk '$6==0{print $1,$2,$3,$0;}$6==16{print $1,$2,$3,$0;}' | wzseqtk.py getfasta -i - -f $BISCUIT_FASTA | cut -f4- >ref_AB_clean
  ## flip ext if need be, add color channel
  awk 'NR==FNR{ext[$5]=$29;tgt[$5]=$30;}NR!=FNR{if($5 in ext) {extn=ext[$5]; tgtn=tgt[$5];} else {extn="NA"; tgtn="NA";} print $0,extn,tgtn;}' ref_AB_clean ref_AB | awk -f wanding.awk -e '{if($6==16 && $29!="NA") {$29=dnarev($29);} if($29=="C" && ($5~/_[TBN]O/)) $29="Y"; if($29=="G" && !($5~/_[TBN]O/)) $29="R"; if($29=="C" || $29=="G") col="G"; else if ($29=="NA") col="NA"; else col="R"; print $0"\t"col;}' | sortbed >ref_AB_final
  cd $BASEDIR
}

function format_mapping_InfiniumII {

  cd $TMPFDR
  echo "=== Format Inf-II ===="
  ## reformat A
  grep '^cg' ref.sam | awk '$2<256' | awk '{$2=and($2,0x14);print $0;}' | awk '!($1~/_1$/ || $1~/_2$/)' | awk -F" " -v OFS="\t" '$2==0 {if($1~/_[TBN]O/){print $3,$4+49,$4+51,substr($10,50,1),$0;}else{print $3,$4+48,$4+50,substr($10,50,1),$0;}} $2==16 {if($1~/_[TBN]O/){print $3,$4-3,$4-1,substr($10,1,1),$0;}else{print $3,$4-2,$4,substr($10,1,1),$0;}} $2==4 {print "*",0,0,substr($10,50,1),$0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >ref_II
  grep '^ch' ref.sam | awk '$2<256' | awk '{$2=and($2,0x14);print $0;}' | awk '!($1~/_1$/ || $1~/_2$/)' | awk -F" " -v OFS="\t" '$2==0 {print $3,$4+49,$4+50,substr($10,50,1),$0} $2==16 {print $3,$4-2,$4-1,substr($10,1,1),$0;} $2==4 {print "*",0,0,substr($10,50,1),$0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >>ref_II
  grep '^rs' ref.sam | awk '$2<256' | awk '{$2=and($2,0x14);print $0;}' | awk '!($1~/_1$/ || $1~/_2$/)' | awk -F" " -v OFS="\t" '$2==0 {print $3,$4+49,$4+50,substr($10,50,1),$0} $2==16 {print $3,$4-2,$4-1,substr($10,1,1),$0;} $2==4 {print "*",0,0,substr($10,50,1),$0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >>ref_II
  
  ## add target sequence
  awk '$1!="*" && $10=="50M"' ref_II | wzseqtk.py getfasta -i - -f $BISCUIT_FASTA >ref_II_clean

  ## add tgt
  awk 'NR==FNR{tgt[$5]=$15;}NR!=FNR {if($5 in tgt) {tgtn=tgt[$5];} else {tgtn="NA";} print $0,tgtn;}' ref_II_clean ref_II | sortbed >ref_II_final
  cd $BASEDIR
}

function validate_Infinium {
  mkdir -p $BASEDIR/validation
  cd $BASEDIR
  echo "=== Create validation files ===="
  ## I
  ## cg and ch
  awk '$1!="*" && ($10=="50M" || $24=="50M") && !($5~/^rs/){$2=$2-50;$3=$3+50;print $0;}' $TMPFDR/ref_AB | wzseqtk.py getfasta -i - -f $BISCUIT_FASTA | awk -f wanding.awk -e 'NR!=FNR{print $5, "I", s1[$5]; print joinr(1,28); print $29; if(!and($6,0x10)){printf(">");if(!($5~/_[TBN]O/)){printf(" ");}} else {printf("<"); for(i=1;i<50;++i) printf(" "); if($5~/_[TBN]O/){printf(" ");}} print $11; if(!and($6,0x10)) {printf(">"); if(!($5~/_[TBN]O/)){printf(" ");} os=s1[$5];}else {printf("<"); for(i=1;i<50;++i) printf(" "); os=dnarev(s1[$5]); if($5~/_[TBN]O/){printf(" ");}} print os; if(!and($6,0x10)){printf(">"); if(!($5~/_[TBN]O/)){printf(" ");}} else {printf("<"); for(i=1;i<50;++i) printf(" "); if($5~/_[TBN]O/){printf(" ");}} print $25; if(!and($6,0x10)) {printf(">"); if(!($5~/_[TBN]O/)){printf(" ");} os=s2[$5];} else {printf("<"); for(i=1;i<50;++i) printf(" "); os=dnarev(s2[$5]); if($5~/_[TBN]O/){printf(" ");}} print os;print("\n");}NR==FNR{s1[$1]=$2;s2[$1]=$3;}' $TMPFDR/fa/probe2originalseq.txt - >validation/${PLATFORM}_${REFCODE}.txt
  ## rs
  awk '$1!="*" && ($10=="50M" || $24=="50M") && $5~/^rs/{$2=$2-50;$3=$3+50;print $0;}' $TMPFDR/ref_AB | wzseqtk.py getfasta -i - -f $BISCUIT_FASTA | awk -f wanding.awk -e 'NR!=FNR{print $5, "I", s1[$5]; print joinr(1,28); print $29; if(!and($6,0x10)){printf(" ");} else {printf("<"); for(i=1;i<50;++i) printf(" "); if($5~/_[TBN]O/){printf(" ");}} print $11; if(!and($6,0x10)) {printf(" "); os=s1[$5];}else {printf("<"); for(i=1;i<50;++i) printf(" "); os=dnarev(s1[$5]); if($5~/_[TBN]O/){printf(" ");}} print os; if(!and($6,0x10)){printf(" ");} else {printf("<"); for(i=1;i<50;++i) printf(" "); if($5~/_[TBN]O/){printf(" ");}} print $25; if(!and($6,0x10)) {printf(" "); os=s2[$5];} else {printf("<"); for(i=1;i<50;++i) printf(" "); os=dnarev(s2[$5]); if($5~/_[TBN]O/){printf(" ");}} print os;print("\n");}NR==FNR{s1[$1]=$2;s2[$1]=$3;}' $TMPFDR/fa/probe2originalseq.txt - >>validation/${PLATFORM}_${REFCODE}.txt  

  ## II
  ## cg
  awk '$1!="*" && $10=="50M" && $5~/^cg/{$2=$2-50;$3=$3+50;print $0;}' $TMPFDR/ref_II | wzseqtk.py getfasta -i - -f $BISCUIT_FASTA | awk -f wanding.awk -e 'NR!=FNR{print $5, "II", s1[$5]; print joinr(1,14); print $15; if(!and($6,0x10)){printf(">");} else {printf("<"); for(i=1;i<51;++i) printf(" "); } print $11; if(!and($6,0x10)) {printf(">"); os=s1[$5];}else {printf("<"); for(i=1;i<51;++i) printf(" "); os=dnarev(s1[$5])} print os; print("\n");}NR==FNR{s1[$1]=$2}' $TMPFDR/fa/probe2originalseq.txt - >>validation/${PLATFORM}_${REFCODE}.txt
  ## ch and rs
  awk '$1!="*" && $10=="50M" && !($5~/^cg/){$2=$2-50;$3=$3+50;print $0;}' $TMPFDR/ref_II | wzseqtk.py getfasta -i - -f $BISCUIT_FASTA | awk -f wanding.awk -e 'NR!=FNR{print $5, "II", s1[$5]; print joinr(1,14); print $15; if(!and($6,0x10)){printf("")} else {printf("<"); for(i=1;i<51;++i) printf(" "); } print $11; if(!and($6,0x10)) {os=s1[$5];}else {printf("<"); for(i=1;i<51;++i) printf(" "); os=dnarev(s1[$5])} print os; print("\n");}NR==FNR{s1[$1]=$2}' $TMPFDR/fa/probe2originalseq.txt - >>validation/${PLATFORM}_${REFCODE}.txt
  gzip -f validation/${PLATFORM}_${REFCODE}.txt
  cd $BASEDIR
}

function count_manifest {
  echo -e "=== Build tsv manifest "tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz"===="
  echo "Original manifest rows: " $(zcat -f $TMPFDR/fa/standard_input.tsv | wc -l)
  echo "Final manifest rows non-ctl: " $(zcat -f tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | cut -f9 | awk 'NR>1' | grep -v 'ctl' |  wc -l)
  echo "Final manifest rows non-ctl (uniq): " $(zcat -f tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | cut -f9 | awk 'NR>1' | sort | uniq | grep -v 'ctl' |  wc -l)
  echo "Final manifest rows cg: " $(zcat -f tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | cut -f9 | awk 'NR>1' | grep 'cg' |  wc -l)
  echo "Final manifest rows cg (uniq): " $(zcat -f tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | cut -f9 | awk 'NR>1' | grep 'cg' | sort | uniq |  wc -l)
  echo "Final manifest rows ch: " $(zcat -f tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | cut -f9 | awk 'NR>1' | grep 'ch' |  wc -l)
  echo "Final manifest rows ch (uniq): " $(zcat -f tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | cut -f9 | awk 'NR>1' | grep 'ch' | sort | uniq |  wc -l)
  echo "Final manifest rows rs: " $(zcat -f tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | cut -f9 | awk 'NR>1' | grep 'rs' |  wc -l)
  echo "Final manifest rows rs (uniq): " $(zcat -f tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | cut -f9 | awk 'NR>1' | grep 'rs' | sort | uniq |  wc -l)
  echo "Final manifest rows ctl: " $(zcat -f tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | cut -f9 | awk 'NR>1' | grep 'ctl' |  wc -l)
  echo "Final manifest rows ctl (uniq): " $(zcat -f tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | cut -f9 | awk 'NR>1' | grep 'ctl' | sort | uniq |  wc -l)
  echo "Check C/O tag (only valid when using the new ID system):"
  zcat tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | awk 'NR>1' | awk '$9~/cg/{split($9,a,"_"); aa=substr(a[2],2,1); if(aa=="C") { if((and($10,0x10) && $18=="f") || (!and($10,0x10) && $18=="r")) {ok="ok";}else{ok="not ok";} } if(aa=="O") { if((!and($10,0x10) && $18=="f") || (and($10,0x10) && $18=="r")) {ok="ok";}else{ok="not ok";} } print aa,$10,$18,ok;}' | sort | uniq -c
  echo "Check 1/2 tag (only valid when using the new ID system):"
  zcat tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | awk 'NR>1&&$9~/cg/{if($5=="NA"){st="II";}else{st="I";} split($9,a,"_"); aa=substr(a[2],3,1); if(aa==1){bb="I";} else {bb="II";} if(st==bb){ok="ok";}else{ok="not ok";} print st,bb,ok;}' | sort | uniq -c
}

function set_features_hg38 {
  export FEATURES="Blacklist.20220304 ChromHMM.20220303 H3K9me3HighConfidence.20220801 ImprintingDMR.20220818 rmsk1.20220307 rmsk2.20220321 Tetranuc2.20220321 CGI.20220904"
}

function set_features_mm10 {
  export FEATURES="ChromHMM.20220414"
}

function buildFeatureOverlaps {

  echo "=== Create annotation files ===="
  mkdir -p $TMPFDR
  mkdir -p features/${PLATFORM}_${REFCODE}/
  [[ -f tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz ]] || exit 1;
  zcat tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | awk 'NR>1&&$1!="NA"&&$9~/^c[gh]/' | sortbed >$TMPFDR/${PLATFORM}_${REFCODE}.bed
  
  for f in $FEATURES; do
    echo -n Processing feature $f;
    bedtools intersect -b $TMPFDR/${PLATFORM}_${REFCODE}.bed -a ~/references/${REFCODE}/features/$f.bed.gz -sorted -wo | awk 'BEGIN{print "Probe_ID\tKnowledgebase"}{print $14,$4;}' | gzip -c >features/${PLATFORM}_${REFCODE}/$f.gz
    echo " (captured" $(zcat features/${PLATFORM}_${REFCODE}/$f.gz | wc -l) "rows)"
  done

  echo -n Processing probe type;
  zcat tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | awk 'NR==1{print $9,"Knowledgebase";}NR>1{print $9,"ProbeType;"substr($9,1,2);}' | gzip -c >features/${PLATFORM}_${REFCODE}/ProbeType.gz
  echo " (captured the following probe types)"
  zcat features/${PLATFORM}_${REFCODE}/ProbeType.gz | awk 'NR>1' | cut -f2 | uniq -c
  echo "Finished processing all features."
}

function buildFeatureGene {
  echo "=== Create gene annotation file ===="
  zcat $GMGTF | awk '(!/^#/)&&$3=="transcript"{match($9,/gene_name "([^"]*)"/, genename); match($9,/transcript_id "([^"]*)"/, transid); match($9,/transcript_type "([^"]*)"/, transcripttype); print $1,$4,$5,$7,genename[1],transcripttype[1],transid[1]}' | awk -f wanding.awk -e '$4=="+"{print $1,max($2-1500,0),$3,joinr(4,7),$2}$4=="-"{print $1,max($2,0),$3+1500,joinr(4,7),$3}' | sortbed >$TMPFDR/$GMCODE.wtss1500.bed
  
  zcat tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | awk 'NR>1 && $1!="NA"{if(and($10, 0x10))s="-"; else s="+";print $1,$2,$3,s,$9;}' | sortbed | bedtools intersect -a - -b $TMPFDR/$GMCODE.wtss1500.bed -sorted -loj | awk '$10!="."{if ($9=="+") a=$2+1-$13; else a=$13-$2; b=$6":"$7"-"$8$9; print $1,$2+1,$3,$4,$5,$10,$11,$12,a,b}' | sort -k1,1 -k2,2n -k4,7 | bedtools groupby -i - -g 1-4 -c 5,6,6,7,8,9 -o distinct,distinct,collapse,collapse,collapse,collapse -delim ";" | awk -f wanding.awk -e 'NR==FNR{gene[$5]=joinr(6,10)}NR!=FNR&&FNR>1{if($9 in gene){g=gene[$9];}else{g="NA\tNA\tNA\tNA\tNA"} if(and($10,0x10)) {st="-";}else{st="+";} print $1,$2,$3,st,$9,g;}' - <(zcat tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz) | awk -f wanding.awk -e 'BEGIN{print "CpG_chrm\tCpG_beg\tCpG_end\tprobe_strand\tprobeID\tgenesUniq\tgeneNames\ttranscriptTypes\ttranscriptIDs\tdistToTSS"}1' | gzip -c >features/${PLATFORM}_${REFCODE}/${PLATFORM}.${REFCODE}.manifest.${GMCODE}.tsv.gz
  echo "Created: "features/${PLATFORM}_${REFCODE}/${PLATFORM}.${REFCODE}.manifest.${GMCODE}.tsv.gz
}

