Table of Contents
=================

<!--ts-->
   * [Table of Contents](#table-of-contents)
      * [Unix Tutorial](#unix-tutorial)
      * [Environment Variables](#environment-variables)
      * [Symlinks](#symlinks)
      * [Job Submission and Monitoring](#job-submission-and-monitoring)
      * [Script-less Submission](#script-less-submission)
      * [pbsgen-style submission](#pbsgen-style-submission)
      * [Keep Job Running After Disconnection](#keep-job-running-after-disconnection)
      * [Submission Script from Scratch](#submission-script-from-scratch)
      * [Project Organization](#project-organization)
      * [Project Documentation](#project-documentation)
      * [Two Shared Lab Folders](#two-shared-lab-folders)
      * [Legacy Data Folders](#legacy-data-folders)
      * [How to Mount Network Drives?](#how-to-mount-network-drives)
      * [Transfer files from/to HPC](#transfer-files-fromto-hpc)
      * [Reference Genome Folder](#reference-genome-folder)
      * [Useful Tools](#useful-tools)

<!-- Added by: zhouw3, at: Sun Apr 26 09:59:15 EDT 2020 -->

<!--te-->

## Use VDI

[https://wiki.chop.edu/pages/viewpage.action?pageId=238785326](https://wiki.chop.edu/pages/viewpage.action?pageId=238785326)


## Unix Tutorial

[http://www.ee.surrey.ac.uk/Teaching/Unix/](http://www.ee.surrey.ac.uk/Teaching/Unix/)

Specs of Respublica, how many nodes, cpus.

[https://wiki.chop.edu/pages/viewpage.action?spaceKey=RISUD&title=Basic+Cluster+Information](https://wiki.chop.edu/pages/viewpage.action?spaceKey=RISUD&title=Basic+Cluster+Information)

## Environment Variables

Append
```
source /mnt/isilon/zhoulab/labtools/bashrc/chop/bashrc_hpc_zhoulab
```

to `~/.bashrc`.

## Symlinks
Symlinks are great ways to keep your path simple and clean. The real path can be seen with `readlink -f`. Here are some common symlinks:

- Genome sequence and annotations:

```
ln -s /mnt/isilon/zhou_lab/projects/20191221_references ~/references
```
- Lab softwares:

```
ln -s /mnt/isilon/zhoulab/labsoftware ~/software
```

- Shared lab storage (can be mounted as network disk):

```
ln -s /mnt/isilon/zhoulab ~/zhoulab
```
- Shared lab storage (cannot be mounted as network disk):

```
ln -s /mnt/isilon/zhou_lab ~/zhou_lab
```
- Personal scratch space (faster in IO, but cannot be mounted as network disk, replace CHOPID by your ID): 

```
ln -s /scr1/users/<CHOPID> ~/scr1_<CHOPID>
```

## Job Submission and Monitoring

We have the follow repo cloned to `/mnt/isilon/zhoulab/labpipelines` for job submission tools

[https://github.com/zhou-lab/labpipelines](https://github.com/zhou-lab/labpipelines)

Alias are defined in `/mnt/isilon/zhoulab/labtools/bashrc/chop/bashrc_hpc_zhoulab` for quick job submission, deletion, monitoring

```
alias qsubi='qlogin -q interactive.q' # interactive job
alias qstatall='qstat -u "*" | less'  # check all job status
alias qstatallrun='qstat -u "*" -s r | less' # check all user jobs
alias qhost='qhost | less' # check queue status
alias qwatch="watch qstat" # keep monitoring jobs
alias qsub1="qsub -pe smp 1 -l m_mem_free=5G -l h_vmem=5G"
alias qsub4="qsub -pe smp 4 -l m_mem_free=5G -l h_vmem=5G"
alias qsub12="qsub -pe smp 12 -l m_mem_free=5G -l h_vmem=5G"
alias qsub24="qsub -pe smp 24 -l m_mem_free=5G -l h_vmem=5G"
export HPCUSERNAME=zhouw3
qdelall
```

qstat on running jobs only

```
qstatz
qwatchz
```

Execute one job (more about pbsgen below)

```
pbsgen "<your command>" -submit
```

Submit multiple jobs with script

```
find folder/ -type f -name '*.pbs' | sort | xargs -I {} qsub {}
```

Find out all jobs running and run

```
qacct -j -o zhouw3
qacct -j -o zhouw3 -d 1 # just since yesterday
```


## Script-less Submission

Notice in the examples below `qsub` can be replaced with `qsub1`, `qsub4`, `qsub12` and `qsub24`.

Here string (importantly no space before and after `EOF`, see [here](https://linuxize.com/post/bash-heredoc/) for some tutorial about heredoc).
```
qsub12 <<'EOF'
<your command>
EOF
```

Pipe in
```
cat <<'EOF' | qsub
<your command>
EOF
```
You can also pipe into both a file and qsub (so that you keep a record)
```
cat <<'EOF' | tee <your file name> | qsub
<your command>
EOF
```

you can also use `echo` if it's just one line.
```
echo <your command> | qsub
```

For R-command, you can do
```
qsub <<'EOF'
Rscript <<'EOF2'
<your R command>
EOF2
EOF
```

## pbsgen-style submission

pbsgen gives more control over the pbs file generated, an example:

```
cat <<'EOF' | pbsgen -submit
<your bash code>
EOF
```

An example for R

```
cat <<'EOF' | pbsgen -submit
Rscript - <<'EOF2'
<your R code>
EOF2
EOF
```

Common options:

- `-submit`: submit the generated script
- `-ppn 4`: number of cores (in this case 4 cores, default to 1)
- `-name`: name of the pbs file and job name. if not an absolute path, use current working directory. if absolute path, the basename will be used for job name.
- `-pbsdir`: name of the default pbs folder. if `$PBSDIR` is not set, use current working directory as default, otherwise use `$PBSDIR` as default.

Example 1: You just don't care where the script file is

```
pbsgen "echo Hello world"
```

This creates job `j<i>_$NAMEROOT` at`$PBSDIR/j<i>_$NAMEROOT.pbs`. `<i>` is auto-incremented.

Example 2: You want specify job name but don't care where.
```
pbsgen "echo Hello World" -name Pearland
```

This creates job `Pearland` at `$PBSDIR/Pearland.pbs`

Example 3: You want to specify both script path and job name.
```
pbsgen "echo Hello World" -name ~/test/Pearland
```

This creates job `Pearland` at `~/test/Pearland.pbs`

Example 4: You want script folder (like current dir) but don't care name
```
pbsgen "echo Hello World" -pbsdir .
```

This creates job `j<i>_$NAMEROOT` at`./j<i>_$NAMEROOT.pbs`. `<i>` is auto-incremented.

Environment variables control the defaults

```
# change this if you want to just change the folder where pbs file is auto-generated
export PBSDIR=/mnt/isilon/zhoulab/tmp/pbs
# change this if you don't like the job name
export NAMEROOT=LabJob
```

Customize these in your `~/.bashrc` file after loading the zhoulab file.

## Keep Job Running After Disconnection

Use `screen` and run everything inside.

- `F2` new panel
- `F11`/`F12` switch left and right
- `screen -r(sr)` reattach
- `Ctrl-a d` detach
- `Ctrl-a K(cap K)` kill the current window
- `Ctrl-a [` copy mode and use `shift+up/down` to scroll up and down
- `Ctrl-a A` set current window title

For more see [video tutorial](https://www.youtube.com/watch?v=HomIzLB-HBc)

## Submission Script from Scratch

Most cases you can use `pbsgen` and `qsubi`, `qsub1-24` for auto-generated submission script (see above). But more details of the submission script can be found at 
[https://wiki.chop.edu/display/RISUD/Grid+Engine](https://wiki.chop.edu/display/RISUD/Grid+Engine). The job specs can be placed either on the command line at the head of the script. There is no time limit.

```
qsub -l h_vmem=4G -l m_mem_free=4G -pe smp 2 script.sh
```

**This is an important thing about SMP, whatever -l <resource> request you make, is multiplied by your SMP count. So if I want 4 cores but 32GB memory, I need to submit with -l h_vmem=8G -l m_mem_free=8G -pe smp 4.** With `pbsgen` it is specified with `-memG 8` option.

## Project Organization
Your project workspace should ideally be sitting at `~/zhou_lab/projects/`.
It'd be better you follow the nomenclature starting with a date when creating your project folder, like `20200102_SPLiTseq_mouse_brain` and `20200106_human_WGBS`.

## Project Documentation
Please document every command needed and working directory for analysis. Create your git repository in `~/zhoulab/labprojects`. See my in `zhouw3` for some examples.

```
cd `/zhoulab/labprojects
mkdir <CHOPID>
cd <CHOPID>
git init
```

Git tutorial: [https://www.youtube.com/watch?v=HVsySz-h9r4](https://www.youtube.com/watch?v=HVsySz-h9r4)

## Two Shared Lab Folders
There are two shared lab folders `/mnt/isilon/zhoulab/` and `/mnt/isilon/zhou_lab`. 
Sorry for the confusing nomenclature but `zhoulab` can be mounted as a network disk on your local computer which means you don't need to sync files back and forth. You can use exactly the same path on HPC and on your local computer by creating a symlink. For example, one my Mac, I have
`ln -s /Volumes/zhoulab/ /mnt/isilon/zhoulab`

But because of that functionality, `zhoulab` has NO write protection, meaning that important data can get deleted at one mistake! I am now syncing the important data to `zhou_lab` which raw data will be kept read-only, just to add a layer of safety.

## Legacy Data Folders
There are three of them `~/zhou_lab/HFS10T/`, `~/zhou_lab/HFS8T/` and `~/zhou_lab/HFS3T/`. Please make sure you don't write into them. I will also try make them read-only.

## How to Mount Network Drives?
If you use a mac, go to Finder > Go > Connect to Server, then put `smb://ressmb03.research.chop.edu/zhoulab`
Your drive will be at `/Volumes/zhoulab`. I usually also do 

```
sudo mkdir -p /mnt/isilon/
ln -sf /Volumes/zhoulab/ /mnt/isilon/zhoulab
```
so that you can use the same path on HPC and local machine.

## Transfer files from/to HPC
transfer files from HPC to local:

```
scp username@respublica.research.chop.edu:/path/on/hpc/ ~/path/on/local
```
transfer files from local to HPC:

```
scp ~/path/on/local username@respublica.research.chop.edu:/path/on/hpc/ 
```
transfer directory from HPC to local:
```
scp -r ~/path/on/local username@respublica.research.chop.edu:/path/on/hpc/ 
```

## use modules
```
module avail
modulefiles_list
module load gcc/6.4.0
```

## Reference Genome Folder

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

## Useful Tools
- [FZF](https://github.com/junegunn/fzf) - fuzzy search
- [z.sh](https://github.com/rupa/z/blob/master/z.sh) - jump around based on history
- [gh-md-toc](https://github.com/ekalinin/github-markdown-toc) - create TOC for markdown files
- [Renv](https://github.com/viking/Renv) - multiple R installation and switches
