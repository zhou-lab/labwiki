# Table of Content
<!--ts-->
   * [Table of Content](#table-of-content)
   * [JOB SUBMISSION AND MONITORING](#job-submission-and-monitoring)
   * [PROJECT ORGANIZATION](#project-organization)
   * [SYMLINKS](#symlinks)
   * [TWO SHARED LAB FOLDERS](#two-shared-lab-folders)
   * [LEGACY DATA FOLDERS](#legacy-data-folders)
   * [HOW TO MOUNT NETWORK DISK?](#how-to-mount-network-disk)
   * [ENVIRONMENTAL VARIABLES](#environmental-variables)
   * [REFERENCE GENOMES FOLDER](#reference-genomes-folder)
   * [USEFUL COMMAND LINE TOOLS](#useful-command-line-tools)

<!-- Added by: zhouw3, at: Tue Jan 28 12:21:23 EST 2020 -->

<!--te-->

# JOB SUBMISSION AND MONITORING

Clone this for job submission tools

[https://github.com/zhou-lab/labpipelines](https://github.com/zhou-lab/labpipelines)

Make pipelines available for execution

```
export PATH=~/repo/labpipelines/pipelines:$PATH
export WZSEQ_ENTRY=~/repo/labpipelines/entry/wzseq.sh
ln -s `rf ~/repo/labpipelines/pbsgen/pbsgen_respublica.py` ~/bin/pbsgen
```

Alias in `.bashrc` for quick job submission, deletion, monitoring

```
alias qsubi='qlogin -q interactive.q' # interactive job
alias qstatall='qstat -u "*" | less'  # check all job status
alias qstatallrun='qstat -u "*" -s r | less' # check all user jobs
alias qhost='qhost | less' # check queue status
alias qwatch="watch qstat" # keep monitoring jobs
function qdelall {
  qstat | grep 'zhouw3' | awk -F " " '{print $1}' | xargs -I {} qdel {}
}
```

Execute one job

```
pbsgen one "samtools index a_bam.bam" -dest <path_for_script> -submit
```

Submit multiple jobs with script

```
find folder/ -type f -name '*.pbs' | sort | xargs -I {} qsub {}
```

# PROJECT ORGANIZATION
Your project workspace should ideally be sitting at `~/zhou_lab/projects/`.
It'd be better you follow the nomenclature starting with a date when creating your project folder, like `20200102_SPLiTseq_mouse_brain` and `20200106_human_WGBS`.

# SYMLINKS
Symlinks are great ways to keep your path simple and clean. The real path can be seen with `readlink -f`. Here are some common symlinks:

- Genome sequence and annotations: `~/references -> /mnt/isilon/zhou_lab/projects/20191221_references`
- shared lab storage (can be mounted as network disk): `~/zhoulab -> /mnt/isilon/zhoulab`
- shared lab storage (cannot be mounted as network disk): `~/zhou_lab -> /mnt/isilon/zhou_lab`
- personal scratch space (faster in IO, but cannot be mounted as network disk): `~/scr1_zhouw3 -> /scr1/users/zhouw3`
- tools (you should create yours): `~/tools -> ~/zhoulab/HFS10T/2019_08_28_HPC_Laird_secondary/2018_05_02_Wanding_tools`

# TWO SHARED LAB FOLDERS
There are two shared lab folders `/mnt/isilon/zhoulab/` and `/mnt/isilon/zhou_lab`. 
Sorry for the confusing nomenclature but `zhoulab` can be mounted as a network disk on your local computer which means you don't need to sync files back and forth. You can use exactly the same path on HPC and on your local computer by creating a symlink. For example, one my Mac, I have
`ln -s /Volumes/zhoulab/ /mnt/isilon/zhoulab`

But because of that functionality, `zhoulab` has NO write protection, meaning that important data can get deleted at one mistake! I am now syncing the important data to `zhou_lab` which raw data will be kept read-only, just to add a layer of safety.

# LEGACY DATA FOLDERS
There are three of them `~/zhoulab/HFS10T/`, `~/zhoulab/HFS8T/` and `~/zhoulab/HFS3T/`. Please make sure you don't write into them. I will also try make them read-only.

# HOW TO MOUNT NETWORK DISK?
If you use a mac, go to Finder > Go > Connect to Server, then put `smb://ressmb03.research.chop.edu/zhoulab`
Your drive will be at `/Volumes/zhoulab`.

# ENVIRONMENTAL VARIABLES
These are the ones I use (you can consider putting them to your `~/.bashrc`, obviously with replacement of your user names)

```
alias rm='rm -i'
alias lc="wc -l"
alias ll="ls -l"
alias parallel="parallel --gnu --progress"
alias scp='rsync -Pravdtze ssh'
alias awk='awk -F"\t" -v OFS="\t"'
alias les="less -S"
alias rdf="readlink -f"

export PATH=~/bin:~/local/bin:$PATH
```

# REFERENCE GENOMES FOLDER

References genome is shared among users. Let's all agree to use the following link for now.

`~/references -> /mnt/isilon/zhou_lab/projects/20191221_references`

All genome assembly is organized by their name (UCSC id if available, Ensembl id if not). 
Underneath each folder like `~/references/hg38/` you will find annotation which contains the annotation of that genome including cpg island, etc. 
Index for each software will be contained in its own folder like `~/references/hg38/biscuit`.

```
.
├── annotation
│   ├── cytoband
│   │   └── cytoBand.txt.gz
│   ├── rmsk
│   │   ├── rmsk.comp.txt
│   │   ├── rmsk.num_cpg.txt
│   │   ├── rmsk.txt.bed
│   │   ├── rmsk.txt.gz
│   │   ├── rmsk2.txt.bed
│   │   └── rmsk_hg38.gtf
│   └── transcripts
│       ├── CCDS.20180614.release22.txt
│       ├── gencode.v28.annotation.gff3.gz
│       ├── gencode.v28.annotation.gtf
│       ├── gencode.v28.annotation.gtf.gz
│       ├── gencode.v28.annotation.gtf.havana_clean.bed
│       ├── gencode.v28.annotation.gtf.transcript.bed
│       ├── gencode.v28.annotation.gtf.tss.bed
│       ├── gencode.v28.annotation.gtf.tss.lincRNA.bed
│       └── gencode.v28.annotation.gtf.tss.protein_coding.bed
├── biscuit
│   ├── hg38.fa.bis.amb
│   ├── hg38.fa.bis.ann
│   ├── hg38.fa.bis.pac
│   ├── hg38.fa.dau.bwt
│   ├── hg38.fa.dau.sa
│   ├── hg38.fa.par.bwt
│   └── hg38.fa.par.sa
├── composition
│   └── hg38.fa.comp
├── hg38.fa
├── hg38.fa.fai
└── liftOver
    └── hg19ToHg38.over.chain.gz
```

# USEFUL TOOLS
- [FZF](https://github.com/junegunn/fzf) - fuzzy search
- [z.sh](https://github.com/rupa/z/blob/master/z.sh) - jump around based on history
- [gh-md-toc](https://github.com/ekalinin/github-markdown-toc) - create TOC for markdown files
