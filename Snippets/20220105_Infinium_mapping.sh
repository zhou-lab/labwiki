#NOTE: Run on compute node with 110GB of RAM (Srun24) if building overlaps with large feature files such as TFBSrm.20221005.bed.gz

function set_environment {
  export BISCUIT_FASTA=$1
  export BASEDIR=$2
  export REFCODE=$3             # hg38, mm10, ...
  export PLATFORM=$4            # how shall we call this new platform
  export CSV=$5                 # original csv file
  export parser=$6              # which parser to use
  export GMGTF=$7               # gene model gtf.gz
  export GMCODE=$8              # gene model name we should call the gene annotation file
  export SNPMERGER=$9
  export TMPFDR=$BASEDIR/tmp/${PLATFORM}_${REFCODE}
}

function run_pipeline_build_Infinium_20220823() ( ## note that fun() (...) runs in a subshell
  set_environment $1 $2 $3 $4 $5 $6 $7 $8 $9
  cd $BASEDIR
  rm -rf $TMPFDR

  runpipe1 | tee $BASEDIR/tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz.log
  echo "Temp folder:" $TMPFDR
  echo "Please delete this folder after checking"
)

function run_pipeline_build_Infinium_20221128() (
  set_environment $1 $2 $3 $4 $5 $6 $7 $8 $9
  cd $BASEDIR
  rm -rf $TMPFDR

  runpipe0 | tee $BASEDIR/tsv_manifest/${REFCODE}.tsv.gz.log
  echo "Temp folder:" $TMPFDR
  echo "Please delete this folder after checking"
)

function runpipe0 {
  csv2standard_input_tsv_$parser
  prepare_fa
  biscuit_mapping_thread5
  format_mapping_InfiniumI
  format_mapping_InfiniumII
  ~/repo/labwiki/Snippets/20220105_Infinium_mapping_mergeV2.R $TMPFDR ref $BASEDIR/tsv_manifest/${REFCODE}.tsv.gz NA # build manifest.tsv.gz
  count_manifest $BASEDIR/tsv_manifest/${REFCODE}.tsv.gz # validate the manifest counts
}

function runpipe1 {
  csv2standard_input_tsv_$parser
  prepare_fa
  biscuit_mapping_thread5
  format_mapping_InfiniumI
  format_mapping_InfiniumII
  validate_Infinium             # create validation drawings
  ~/repo/labwiki/Snippets/20220105_Infinium_mapping_mergeV2.R $TMPFDR ref $BASEDIR/tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz NA # build manifest.tsv.gz
  count_manifest $BASEDIR/tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz # validate the manifest counts
  
  create_masks_step1
  for SNP in ~/references/$REFCODE/annotation/snp/snp_bed_standard/*; do create_mask_snp; done
  merge_masks_$SNPMERGER

  set_features_$REFCODE
  buildFeatureGenome
  buildFeatureTechnical
  buildFeatureStudy
  buildFeatureGene
}

function create_masks_step1 {
  echo "=== Create masks ===="
  rm -f $TMPFDR/mask_*
  zcat tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | awk '$9~/^uk/' | awk '{print $9,"M_designUK";}' >$TMPFDR/mask_designUK.tsv
  zcat tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | awk '$9~/^c[gh]/' | awk '$1=="NA"||$17=="NA"||$17<35||($26!="NA"&&$26<35)||($21!="NA"&&$12!=$21)' | awk '{print $9,"M_mapping"}' >$TMPFDR/mask_mapping.tsv
  zcat tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | awk '$9~/^c[gh]/' | awk 'NR==FNR{a[$1]}NR!=FNR&&!($9 in a)' $TMPFDR/mask_mapping.tsv - | sortbed | gzip -c >$TMPFDR/cg_and_ch.bed.gz
  zcat $TMPFDR/cg_and_ch.bed.gz | awk '$13<=10||($22!="NA"&&$22<=10){print $9,"M_nonuniq"}' >$TMPFDR/mask_nonunique.tsv
}

function create_mask_snp {
  echo "=== Create SNP masks "$SNP
  mname=$(basename $SNP .bed.gz)
  zcat $TMPFDR/cg_and_ch.bed.gz | awk '{if(and($10,0x10)){st="-";} else {st="+";} print $11,$12-1,$12+49,st,$9,$18;}' | awk '$4=="+"{$2=$3-5;print;}$4=="-"{$3=$2+5;print;}' | sortbed | bedtools intersect -a - -b $SNP -sorted -wo | awk -v mname=$mname '{if($4=="-"){d=$8-$2;}else{d=$3-$9;} print $5,"M_"mname";"$10";"$11$12";"$6";"d";"$13;}' >$TMPFDR/mask_$mname.tsv
  zcat $TMPFDR/cg_and_ch.bed.gz | awk '{if(and($10,0x10)){st="-";} else {st="+";} if($5=="NA"){type="II";}else{type="I";} print $11,$12-1,$12+49,st,$9,$18,type;}' | awk '$4=="+"{$2=$3;$3=$3+1;print;}$4=="-"{$3=$2;$2=$2-1;print;}' | sortbed | bedtools intersect -a - -b $SNP -sorted -wo -sorted | awk '{print $5,$11";"$12$13";"$6";"$7";"$14,$12,$13,$6,$7;}' >$TMPFDR/tmp_maskExtention_$mname.tsv
  awk -v mname=$mname '$6=="I"{if($5=="f"){if($3=="C")$3="T";if($4=="C")$4="T";} if($5=="r"){if($3=="G")$3="A";if($4=="G")$4="A";} if(($3~/[CG]/ && $4~/[TA]/)||($4~/[CG]/ && $3~/[TA]/)) {print $1,"M_1baseSwitch"mname";"$2;}}' $TMPFDR/tmp_maskExtention_$mname.tsv >$TMPFDR/mask_${mname}_1baseSwitch.tsv
  awk -v mname=$mname '$6=="II"{print $1,"M_2extBase_"mname";"$2;}' $TMPFDR/tmp_maskExtention_$mname.tsv >$TMPFDR/mask_${mname}_2extBase.tsv
}

function merge_masks_hg38dbSNP20180418 {
  find $TMPFDR/ mask_fromStudies -name 'mask_*.tsv' | xargs cat | awk '{split($2,a,";"); print $1,$2,a[1];}' | sort -k1,1 | bedtools groupby -g 1 -c 2,3 -o collapse,distinct | awk 'BEGIN{print "Probe_ID\tmask\tmaskUniq\tM_general";}{split($3,a,",");g="FALSE";for(i=1;i<=length(a);++i){if(a[i]=="M_mapping"||a[i]=="M_nonuniq"||a[i]=="M_SNPcommon_5pt"||a[i]=="M_1baseSwitchSNPcommon_5pt"||a[i]=="M_2extBase_SNPcommon_5pt"){g="TRUE";}} print $0,g;}' | gzip -c >mask/${PLATFORM}.${REFCODE}.mask.tsv.gz
  echo $(zcat mask/${PLATFORM}.${REFCODE}.mask.tsv.gz | awk '$4=="TRUE"' | wc -l)" probes masked."
}

function merge_masks_mm10dbSNP142 {
  echo "=== mask/${PLATFORM}.${REFCODE}.mask.tsv.gz"
  find $TMPFDR/ mask_fromStudies -name 'mask_*.tsv' | xargs cat | awk '{split($2,a,";"); print $1,$2,a[1];}' | sort -k1,1 | bedtools groupby -g 1 -c 2,3 -o collapse,distinct | awk 'BEGIN{print "Probe_ID\tmask\tmaskUniq\tM_general";}{split($3,a,",");g="FALSE";for(i=1;i<=length(a);++i){if(a[i]=="M_mapping"||a[i]=="M_nonuniq"){g="TRUE";}} print $0,g;}' | gzip -c >mask/${PLATFORM}.${REFCODE}.mask.tsv.gz
  echo $(zcat mask/${PLATFORM}.${REFCODE}.mask.tsv.gz | awk '$4=="TRUE"' | wc -l)" probes masked."
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
  zcat -f $CSV | awk -F"," 'a{print "ctl_"$1"_"gensub(" ","_","g",$2),$1,"NA","NA","II";}/^\[Control/{a=1}' >standard_input_control.tsv
  zcat -f $CSV | awk -F"," 'a{print "ctl_"$1"_"gensub(" ","_","g",$2),$1,$2,$3,$4,"II";}/^\[Control/{a=1}' >standard_input_control_anno.tsv
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
  # 1:IlmnID (cgNumber with suffix), 2:Name (cgNumber), 3:AddressA_ID, 4:AlleleA_ProbeSeq, 5:AddressB_ID, 6:AlleleB_ProbeSeq
  zcat -f $CSV | awk 'BEGIN{a=0}/^IlmnID/{a=1;next;}/^\[Controls\]/{a=0;}a' | csvtk csv2tab | awk '{print $1,$15,$3,$4,$5,$6}' >standard_input.tsv
  zcat -f $CSV | awk -F"," 'a{print "ctl_"$1"_"gensub(" ","_","g",$2),$1,"NA","NA","II";}/^\[Control/{a=1}' >standard_input_control.tsv
  zcat -f $CSV | awk -F"," 'a{print "ctl_"$1"_"gensub(" ","_","g",$2),$1,$2,$3,$4,"II";}/^\[Control/{a=1}' >standard_input_control_anno.tsv
  cd $BASEDIR
}

function csv2standard_input_tsv_MSAdraft {
  # requirement: csvtk
  # column requirement:
  mkdir -p $TMPFDR/fa
  cd $TMPFDR/fa
  # 1:cgNumber, 2:I/II, 3:AddressA_ID, 4:AlleleA_ProbeSeq, 5:AlleleB_ID, 6:AlelleB_ProbeSeq
  zcat -f $CSV | awk '{if($3=="NA"){type="II"}else{type="I";}print $1,type,"12345"NR,$2,"12345"NR,$3}' >standard_input.tsv
  :>standard_input_control.tsv
  :>standard_input_control_anno.tsv
  cd $BASEDIR
}

function csv2standard_input_tsv_MSAdraft2 {
  # requirement: csvtk
  # column requirement:
  mkdir -p $TMPFDR/fa
  cd $TMPFDR/fa
  # 1:bucket, 2:Probe_ID, 3:AlleleA_ID 4:AlleleA_ProbeSeq 5:AlleleB_ID 6:AlleleB_ProbeSeq 7:NormalizationBin
  zcat -f $CSV | csvtk csv2tab | awk 'NR>1{if($6==""||$6=="NA"){type="II";$6=""}else{type="I";}print $2,type,"123"NR,$4,"123"NR,$6;}' >standard_input.tsv
  :>standard_input_control.tsv
  :>standard_input_control_anno.tsv
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

function csv2standard_input_tsv_LEGX {
  # requirement: csvtk
  # column requirement:
  mkdir -p $TMPFDR/fa
  cd $TMPFDR/fa
  # 1:ProbeID, 2:I/II inferred from AddressA_ID, 3:AddressA_ID, 4:AlleleA_ProbeSeq, 5:AlleleB_ID, 6:AlelleB_ProbeSeq
  zcat -f $CSV | awk 'NR>1 && !($1~/^ct/ || $1~/^uk/){if($16=="NA" || $16=="") type="II"; else type="I"; if ($2=="NA") $2=""; if ($16=="NA") $16=""; print $1,type,$3,$15,$2,$16}' >standard_input.tsv
  zcat -f $CSV | awk 'NR>1 &&  ($1~/^ct/ || $1~/^uk/){if($16=="NA" || $16=="") type="II"; else type="I"; if ($2=="NA") $2=""; if ($16=="NA") $16=""; print $1,$3,$2,"NA",type}' >standard_input_control.tsv
  zcat -f $CSV | awk 'NR>1 &&  ($1~/^ct/ || $1~/^uk/){if($16=="NA" || $16=="") type="II"; else type="I"; if ($2=="NA") $2=""; if ($16=="NA") $16=""; print $1,$3,"NA","NA","NA","NA";}' >standard_input_control_anno.tsv
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
  grep '^rs\|^nv' ref.sam | awk '$2<256 && $1~/_1$/' | awk '{$2=and($2,0x14);print $0;}' | awk -F" " -v OFS="\t" '!and($2,0x4) && !and($2,0x10){print $3,$4+48,$4+49,substr($10,50,1),$0} !and($2,0x4) && and($2,0x10){print $3,$4-1,$4,substr($10,1,1),$0;} and($2,0x4){print "*",0,0,substr($10,50,1),$0;}' | awk -v OFS="\t" '{$5=substr($5,1,length($5)-2);print $0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >>ref_A
  ## reformat B
  grep '^cg' ref.sam | awk '$2<256 && $1~/_2$/' | awk '{$2=and($2,0x14);print $0;}' | awk -F" " -v OFS="\t" '!and($2,0x4) && !and($2,0x10){if($1~/_[TBN]O/){print $3,$4+48,$4+50,substr($10,50,1),$0;} else {print $3,$4+47,$4+49,substr($10,50,1),$0;}} !and($2,0x4) && and($2,0x10){if($1~/_[TBN]O/){print $3,$4-2,$4,substr($10,1,1),$0;}else{print $3,$4-1,$4+1,substr($10,1,1),$0;}} and($2,0x4){print "*",0,0,substr($10,50,1),$0;}' | awk -v OFS="\t" '{$5=substr($5,1,length($5)-2);print $0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >ref_B
  grep '^rs\|^nv' ref.sam | awk '$2<256 && $1~/_2$/' | awk '{$2=and($2,0x14);print $0;}' | awk -F" " -v OFS="\t" '!and($2,0x4) && !and($2,0x10){print $3,$4+48,$4+49,substr($10,50,1),$0} !and($2,0x4) && and($2,0x10){print $3,$4-1,$4,substr($10,1,1),$0;} and($2,0x4){print "*",0,0,substr($10,50,1),$0;}' | awk -v OFS="\t" '{$5=substr($5,1,length($5)-2);print $0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >>ref_B
  ## merge A,B
  awk 'NR==FNR{a[$5]=$0;}NR!=FNR && ($5 in a){$6=and($6,0x14); print $0"\t"a[$5];}' ref_B ref_A >ref_AB
  
  ## add extension and target sequence
  awk '$1!="*" && $10=="50M"' ref_AB | awk '$6==0{if($5~/_[TBN]O/ && !($5~/^rs/||$5~/^nv/)){pos=$3-1;}else{pos=$3;} print $1,pos,pos+1,$0;}$6==16{if($5~/_[TBN]O/ && !($5~/^rs/||$5~/^nv/)){pos=$2+1;}else{pos=$2;} print $1,pos-1,pos,$0;}' | wzseqtk.py getfasta -i - -f $BISCUIT_FASTA | cut -f4- | awk '$6==0{print $1,$2,$3,$0;}$6==16{print $1,$2,$3,$0;}' | wzseqtk.py getfasta -i - -f $BISCUIT_FASTA | cut -f4- >ref_AB_clean
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
  grep '^rs\|^nv' ref.sam | awk '$2<256' | awk '{$2=and($2,0x14);print $0;}' | awk '!($1~/_1$/ || $1~/_2$/)' | awk -F" " -v OFS="\t" '$2==0 {print $3,$4+49,$4+50,substr($10,50,1),$0} $2==16 {print $3,$4-2,$4-1,substr($10,1,1),$0;} $2==4 {print "*",0,0,substr($10,50,1),$0;}' | awk -f wanding.awk -e '{if(!match(joinr(16,NF),/(NM:[^[:space:]]*)/,nm)){nm1=".";} else nm1=nm[1]; match(joinr(16,NF),/(AS:[^[:space:]]*)/,as); match($0,/(YD:[^[:space:]]*)/,yd); print joinr(1,10),$14,nm1,as[1],yd[1];}' >>ref_II
  
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
  awk '$1!="*" && ($10=="50M" || $24=="50M") && !($5~/^rs/||$5~/^nv/){$2=$2-50;$3=$3+50;print $0;}' $TMPFDR/ref_AB | wzseqtk.py getfasta -i - -f $BISCUIT_FASTA | awk -f wanding.awk -e 'NR!=FNR{print $5, "I", s1[$5]; print joinr(1,28); print $29; if(!and($6,0x10)){printf(">");if(!($5~/_[TBN]O/)){printf(" ");}} else {printf("<"); for(i=1;i<50;++i) printf(" "); if($5~/_[TBN]O/){printf(" ");}} print $11; if(!and($6,0x10)) {printf(">"); if(!($5~/_[TBN]O/)){printf(" ");} os=s1[$5];}else {printf("<"); for(i=1;i<50;++i) printf(" "); os=dnarev(s1[$5]); if($5~/_[TBN]O/){printf(" ");}} print os; if(!and($6,0x10)){printf(">"); if(!($5~/_[TBN]O/)){printf(" ");}} else {printf("<"); for(i=1;i<50;++i) printf(" "); if($5~/_[TBN]O/){printf(" ");}} print $25; if(!and($6,0x10)) {printf(">"); if(!($5~/_[TBN]O/)){printf(" ");} os=s2[$5];} else {printf("<"); for(i=1;i<50;++i) printf(" "); os=dnarev(s2[$5]); if($5~/_[TBN]O/){printf(" ");}} print os;print("\n");}NR==FNR{s1[$1]=$2;s2[$1]=$3;}' $TMPFDR/fa/probe2originalseq.txt - >validation/${PLATFORM}_${REFCODE}.txt
  ## rs
  awk '$1!="*" && ($10=="50M" || $24=="50M") && ($5~/^rs/||$5~/^nv/){$2=$2-50;$3=$3+50;print $0;}' $TMPFDR/ref_AB | wzseqtk.py getfasta -i - -f $BISCUIT_FASTA | awk -f wanding.awk -e 'NR!=FNR{print $5, "I", s1[$5]; print joinr(1,28); print $29; if(!and($6,0x10)){printf(" ");} else {printf("<"); for(i=1;i<50;++i) printf(" "); if($5~/_[TBN]O/){printf(" ");}} print $11; if(!and($6,0x10)) {printf(" "); os=s1[$5];}else {printf("<"); for(i=1;i<50;++i) printf(" "); os=dnarev(s1[$5]); if($5~/_[TBN]O/){printf(" ");}} print os; if(!and($6,0x10)){printf(" ");} else {printf("<"); for(i=1;i<50;++i) printf(" "); if($5~/_[TBN]O/){printf(" ");}} print $25; if(!and($6,0x10)) {printf(" "); os=s2[$5];} else {printf("<"); for(i=1;i<50;++i) printf(" "); os=dnarev(s2[$5]); if($5~/_[TBN]O/){printf(" ");}} print os;print("\n");}NR==FNR{s1[$1]=$2;s2[$1]=$3;}' $TMPFDR/fa/probe2originalseq.txt - >>validation/${PLATFORM}_${REFCODE}.txt  

  ## II
  ## cg
  awk '$1!="*" && $10=="50M" && $5~/^cg/{$2=$2-50;$3=$3+50;print $0;}' $TMPFDR/ref_II | wzseqtk.py getfasta -i - -f $BISCUIT_FASTA | awk -f wanding.awk -e 'NR!=FNR{print $5, "II", s1[$5]; print joinr(1,14); print $15; if(!and($6,0x10)){printf(">");} else {printf("<"); for(i=1;i<51;++i) printf(" "); } print $11; if(!and($6,0x10)) {printf(">"); os=s1[$5];}else {printf("<"); for(i=1;i<51;++i) printf(" "); os=dnarev(s1[$5])} print os; print("\n");}NR==FNR{s1[$1]=$2}' $TMPFDR/fa/probe2originalseq.txt - >>validation/${PLATFORM}_${REFCODE}.txt
  ## ch, rs, and nv
  awk '$1!="*" && $10=="50M" && !($5~/^cg/){$2=$2-50;$3=$3+50;print $0;}' $TMPFDR/ref_II | wzseqtk.py getfasta -i - -f $BISCUIT_FASTA | awk -f wanding.awk -e 'NR!=FNR{print $5, "II", s1[$5]; print joinr(1,14); print $15; if(!and($6,0x10)){printf("")} else {printf("<"); for(i=1;i<51;++i) printf(" "); } print $11; if(!and($6,0x10)) {os=s1[$5];}else {printf("<"); for(i=1;i<51;++i) printf(" "); os=dnarev(s1[$5])} print os; print("\n");}NR==FNR{s1[$1]=$2}' $TMPFDR/fa/probe2originalseq.txt - >>validation/${PLATFORM}_${REFCODE}.txt
  gzip -f validation/${PLATFORM}_${REFCODE}.txt
  cd $BASEDIR
}

function count_manifest {
  mft=$1
  echo -e "=== Build tsv manifest "${mft}"===="
  echo "Original manifest rows: " $(zcat -f $TMPFDR/fa/standard_input.tsv | wc -l)
  echo "Final manifest rows non-ctl: " $(zcat -f ${mft} | cut -f9 | awk 'NR>1' | grep -v 'ctl' |  wc -l)
  echo "Final manifest rows non-ctl (uniq): " $(zcat -f ${mft} | cut -f9 | awk 'NR>1' | sort | uniq | grep -v 'ctl' |  wc -l)
  echo "Final manifest rows cg: " $(zcat -f ${mft} | cut -f9 | awk 'NR>1' | grep 'cg' |  wc -l)
  echo "Final manifest rows cg (uniq): " $(zcat -f ${mft} | cut -f9 | awk 'NR>1' | grep 'cg' | sort | uniq |  wc -l)
  echo "Final manifest rows ch: " $(zcat -f ${mft} | cut -f9 | awk 'NR>1' | grep 'ch' |  wc -l)
  echo "Final manifest rows ch (uniq): " $(zcat -f ${mft} | cut -f9 | awk 'NR>1' | grep 'ch' | sort | uniq |  wc -l)
  echo "Final manifest rows rs: " $(zcat -f ${mft} | cut -f9 | awk 'NR>1' | grep 'rs' |  wc -l)
  echo "Final manifest rows rs (uniq): " $(zcat -f ${mft} | cut -f9 | awk 'NR>1' | grep 'rs' | sort | uniq |  wc -l)
  echo "Final manifest rows nv: " $(zcat -f ${mft} | cut -f9 | awk 'NR>1' | grep 'nv' |  wc -l)
  echo "Final manifest rows nv (uniq): " $(zcat -f ${mft} | cut -f9 | awk 'NR>1' | grep 'nv' | sort | uniq |  wc -l)
  echo "Final manifest rows ctl: " $(zcat -f ${mft} | cut -f9 | awk 'NR>1' | grep 'ctl' |  wc -l)
  echo "Final manifest rows ctl (uniq): " $(zcat -f ${mft} | cut -f9 | awk 'NR>1' | grep 'ctl' | sort | uniq |  wc -l)

  echo "Check C/O tag (only valid when using the new ID system):"
  zcat ${mft} | awk 'NR>1' | awk '$9~/cg/{split($9,a,"_"); aa=substr(a[2],2,1); if(aa=="C") { if((and($10,0x10) && $18=="f") || (!and($10,0x10) && $18=="r")) {ok="ok";}else{ok="not ok";} } if(aa=="O") { if((!and($10,0x10) && $18=="f") || (and($10,0x10) && $18=="r")) {ok="ok";}else{ok="not ok";} } print $9,aa,$10,$18,ok;}' >$TMPFDR/CHECK_CO_tag.tsv
  cut -f2- $TMPFDR/CHECK_CO_tag.tsv | sort | uniq -c
  echo "See $TMPFDR/CHECK_CO_tag.tsv for detail."
  
  echo "Check 1/2 tag (only valid when using the new ID system):"
  zcat ${mft} | awk 'NR>1&&$9~/cg/{if($5=="NA"){st="II";}else{st="I";} split($9,a,"_"); aa=substr(a[2],3,1); if(aa==1){bb="I";} else {bb="II";} if(st==bb){ok="ok";}else{ok="not ok";} print $9,st,bb,ok;}' >$TMPFDR/CHECK_12_tag.tsv
  cut -f2- $TMPFDR/CHECK_12_tag.tsv | sort | uniq -c
  echo "See $TMPFDR/CHECK_12_tag.tsv for detail."
}

function set_features_hg19 {
 :
}

function set_features_hg38 {
    export FEATURES="
ABCompartment.20220911
Blacklist.20220304
CGI.20220904
ChromHMM.20220303
CTCFbind.20220911
HM.20221013
ImprintingDMR.20220818
MetagenePC.20220911
nFlankCG.20220321
PMD.20220911
REMCChromHMM.20220911
rmsk1.20220307
rmsk2.20220321
Tetranuc2.20220321
TFBSrm.20221005
"
}

function set_features_mm10 {
    export FEATURES="
Blacklist.20220304
CGI.20220904
ChromHMM.20220318
Chromosome.20221129
EnsRegBuild.20220710
HM.20221013
MetagenePC.20220911
nFlankCG.20220321
PMD.20220911
rmsk1.20220321
rmsk2.20220321
Tetranuc2.20220321
Tetranuc4.20220321
TFBSrm.20221005
"
}

function buildFeatureGenome {
  echo "=== Create genome-based annotation files ===="

  mkdir -p $TMPFDR
  mkdir -p features/${PLATFORM}/${REFCODE}/
  rm -f features/${PLATFORM}/${REFCODE}/*
  
  [[ -f tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz ]] || exit 1;
  zcat tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | awk 'NR>1&&$1!="NA"&&$9~/^c[gh]/' | sortbed >$TMPFDR/${PLATFORM}_${REFCODE}.bed
  
  for f in ${FEATURES[@]}; do
    echo Processing feature $f;
    # to capture all file parts
    for ff in ~/references/${REFCODE}/features/${f}*.bed.gz; do
      fname=$(basename $ff .bed.gz)
      bedtools intersect -b $TMPFDR/${PLATFORM}_${REFCODE}.bed -a $ff -sorted -wo | awk 'BEGIN{print "Probe_ID\tKnowledgebase"}{print $14,$4;}' | sort | uniq | gzip -c >features/${PLATFORM}/${REFCODE}/$fname.gz
      echo "  1 file captured" $(zcat features/${PLATFORM}/${REFCODE}/$fname.gz | wc -l) "rows."
    done
  done

  # tar -zcvf features/${PLATFORM}.${REFCODE}.annotations.tar.gz -C features/${PLATFORM}_${REFCODE}/ .
  echo "Finished processing all features."
}

function buildFeatureStudy {

  echo "=== Create study-based annotation files ===="
  mkdir -p $TMPFDR
  mkdir -p features/${PLATFORM}/Studies/
  rm -f features/${PLATFORM}/Studies/*
  
  [[ -f tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz ]] || exit 1;
  for p2 in studies/*; do
    for f in $p2/*.gz; do
      zcat $f | awk 'BEGIN{print "Probe_ID\tKnowledgebase"}FNR==NR && FNR>1{a[$1]=$2;} NR!=FNR {if($9~/_/) {split($9,kk,"_"); k=kk[1];} else {k = $9} if(k in a) {print $9,a[k];}}' - <(zcat tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz) | gzip -c >features/${PLATFORM}/Studies/$(basename $p2)_$(basename $f .gz).gz
    done
  done

  echo "Finished processing all study annotations."
}

function buildFeatureTechnical {
  echo "=== Create technical annotation files ===="
  mkdir -p $TMPFDR
  mkdir -p features/${PLATFORM}/Technical
  
  echo -n Processing probe type
  zcat tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | awk 'NR==1{print $9,"Knowledgebase";}NR>1{print $9,"ProbeType;"substr($9,1,2);}' | gzip -c >features/${PLATFORM}/Technical/ProbeType.gz
  echo " (captured the following probe types)"
  zcat features/${PLATFORM}/Technical/ProbeType.gz | awk 'NR>1' | cut -f2 | uniq -c

  ## turn Infinium chemistry into KYCG knowledgebases
  echo -n Processing Infinium chemistry
  zcat tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | awk 'NR==1{print "Probe_ID\tKnowledgebase";}NR>1{if($7=="NA"){type="II";} else {type="I."$7;} print $9,"InfiniumChemistry;"type;}' | gzip -c >features/${PLATFORM}/Technical/InfiniumChemistry.gz
  echo " (captured the following Infinium chemistry)"
  zcat features/${PLATFORM}/Technical/InfiniumChemistry.gz | awk 'NR>1' | cut -f2 | sort | uniq -c

  ## turn masks into KYCG knowledgebases
  zcat mask/${PLATFORM}.${REFCODE}.mask.tsv.gz | awk 'BEGIN{print "Probe_ID\tKnowledgebase";}NR>1{split($3,a,","); for(i in a) {print $1,"Mask;"a[i];}}' | gzip -c >features/${PLATFORM}/Technical/Mask_${REFCODE}.gz
  echo "Finished processing all technical features."
}

function buildFeatureGene {
  echo "=== Create gene annotation file ===="
  mkdir -p gene_assoc
  zcat $GMGTF | awk '(!/^#/)&&$3=="transcript"{match($9,/gene_name "([^"]*)"/, genename); match($9,/transcript_id "([^"]*)"/, transid); match($9,/transcript_type "([^"]*)"/, transcripttype); print $1,$4,$5,$7,genename[1],transcripttype[1],transid[1]}' | awk -f wanding.awk -e '$4=="+"{print $1,max($2-1500,0),$3,joinr(4,7),$2}$4=="-"{print $1,max($2,0),$3+1500,joinr(4,7),$3}' | sortbed >$TMPFDR/$GMCODE.wtss1500.bed
  
  zcat tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz | awk 'NR>1 && $1!="NA"{if(and($10, 0x10))s="-"; else s="+";print $1,$2,$3,s,$9;}' | sortbed | bedtools intersect -a - -b $TMPFDR/$GMCODE.wtss1500.bed -sorted -loj | awk '$10!="."{if ($9=="+") a=$2+1-$13; else a=$13-$2; b=$6":"$7"-"$8$9; print $1,$2+1,$3,$4,$5,$10,$11,$12,a,b}' | sort -k1,1 -k2,2n -k4,7 | bedtools groupby -i - -g 1-5 -c 6,6,7,8,9 -o distinct,collapse,collapse,collapse,collapse -delim ";" | awk -f wanding.awk -e 'NR==FNR{gene[$5]=joinr(6,10)}NR!=FNR&&FNR>1{if($9 in gene){g=gene[$9];}else{g="NA\tNA\tNA\tNA\tNA"} if(and($10,0x10)) {st="-";}else{st="+";} print $1,$2,$3,st,$9,g;}' - <(zcat tsv_manifest/${PLATFORM}.${REFCODE}.manifest.tsv.gz) | awk -f wanding.awk -e 'BEGIN{print "CpG_chrm\tCpG_beg\tCpG_end\tprobe_strand\tprobeID\tgenesUniq\tgeneNames\ttranscriptTypes\ttranscriptIDs\tdistToTSS"}1' | gzip -c >gene_assoc/${PLATFORM}.${REFCODE}.manifest.${GMCODE}.tsv.gz
  echo "Created: "gene_assoc/${PLATFORM}.${REFCODE}.manifest.${GMCODE}.tsv.gz
}

## Example usage: snp_file_anno tsv_manifest/EPICv2.hg38.manifest.tsv.gz /mnt/isilon/zhoulab/labprojects/20200704_dbSNP_mutation/GRCh38p7_b151/00-common_all_SNPonly.sorted.bed.gz snp_anno/EPICv2.hg38.snp.tsv.gz
function snp_file_anno {
  in_manifest=$1
  snp_file=$2
  out_file=$3
  tmp_file=tmp/vcf_format/tmp.bed
  zcat $in_manifest | awk '$9~/^rs/' | sortbed |
    bedtools intersect -a - -b <(zcat $snp_file) -loj -sorted | awk -f wanding.awk -e '
     {if($10=="0") st="+"; else st="-"; probeID=$9; IorII=$28; REF_BYPOS=$6;REF=$33;ALT=$34;
     if ($32==".") { rsID="NA"; } else { rsID=$32"."$33">"$34; }
     if (IorII=="I") {
        U_LAST = substr($15,length($15),1);
        if (st == "-") { U_LAST = dnarev(U_LAST); }
        CREF = REF;
        if (probeID~/_[TB]O/) {
           if (st == "+" && CREF == "C") CREF="T";
           if (st == "-" && CREF == "G") CREF="A";
        } else {
           if (st == "+" && CREF == "G") CREF="A";
           if (st == "-" && CREF == "C") CREF="T";
        }
        if (U_LAST == CREF && U_LAST != ALT) { U="REF"; } else { U="ALT"; }
     } else { # type II
        if (st=="+") {
           if (REF_BYPOS~/[AGT]/) { U="REF"; REF="AGT"; } else { U="ALT"; REF="C";}
        } else {
           if (REF_BYPOS~/[ACT]/) { U="REF"; REF="ACT"; } else { U="ALT"; REF="G";}
        }
     }
     print $1,$2,$3,st,rsID,IorII,U,REF,ALT,probeID;}' >$tmp_file"_1"

  zcat $in_manifest | awk '$6!="NA"' |
    awk '$9~/^cg/&&$28=="I"{if($10=="0"){$2=$3;} else {$2=$2-1;} $3=$2+1; print $0;}' |
    sortbed | bedtools intersect -a - -b <(zcat $snp_file) -loj -sorted | awk -f wanding.awk -e '
     {if($10=="0") st="+"; else st="-"; probeID=$9; IorII=$28; REF_BYPOS=$6; COL=$8;
     if ($32==".") { rsID="NA"; } else { rsID=$32"."$33">"$34; }
     U = "REF_InfI";
     if (probeID~/_[TB]O/) { # assume color channel is G
        if (st == "+") { REF="G"; ALT="ACT"; } else { REF="C"; ALT="AGT"; }
     } else {
        if (st == "+") { REF="C"; ALT="AGT"; } else { REF="G"; ALT="ACT"; }
     }
     if (COL=="R") { TMP=REF; REF=ALT; ALT=TMP; } # swap if channel is R
     print $1,$2,$3,st,rsID,IorII,U,REF,ALT,probeID;}' >$tmp_file"_2"

  cat $tmp_file"_1" $tmp_file"_2" | sortbed |
    awk 'BEGIN{print "chrm\tbeg\tend\tstrand\trs\tdesignType\tU\tREF\tALT\tProbe_ID";}1' |
    gzip -c >$out_file
}
