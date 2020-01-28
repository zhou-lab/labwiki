## SYMLINKS
Symlinks are great ways to keep your path simple and clean. The real path can be seen with readlink -f
genome sequence and annotations:
`~/references -> /mnt/isilon/zhou_lab/projects/20191221_references`

Lab Storage:
- shared (can be mounted as network disk):
`~/zhoulab -> /mnt/isilon/zhoulab`
- shared (cannot be mounted as network disk):
`~/zhou_lab -> /mnt/isilon/zhou_lab`
- personal scratch space (faster in IO, but cannot be mounted as network disk):
`~/scr1_zhouw3 -> /scr1/users/zhouw3`

# TWO SHARED LAB FOLDERS
There are two shared lab folders `/mnt/isilon/zhoulab/` and `/mnt/isilon/zhou_lab`. 
Sorry for the confusing nomenclature but zhoulab can be mounted as a network disk on your local computer which means you don't need to sync files back and forth. You can use exactly the same path on HPC and on your local computer by creating a symlink. For example, one my computer I have
`ln -s /Volumes/zhoulab/ /mnt/isilon/zhoulab`
But because of that functionality `zhoulab` has NO write protection, meaning that important data can get deleted at one mistake! I am now syncing the important data to zhou_lab which raw data will be kept read-only, just to add a layer of safety.

# HOW TO MOUNT NETWORK DISK?
If you use a mac, go to Finder > Go > Connect to Server, then put `smb://ressmb03.research.chop.edu/zhoulab`
Your drive will be at `/Volumes/zhoulab`.

# PROJECT ORGANIZATION
Your project workspace should ideally be sitting at `~/zhou_lab/projects/`
It's better you follow the nomenclature starting with a date when creating your project folder, like `20200102_SPLiTseq_mouse_brain`

# ENVIRONMENTAL VARIABLES
These are the ones I use (you can consider putting them to your ~/.bashrc, obviously with adaptation)
```
alias rm='rm -i'
alias lc="wc -l"
alias ll="ls -l"
alias parallel="parallel --gnu --progress"
export PATH=~/bin:~/local/bin:$PATH
alias scp='rsync -Pravdtze ssh'
alias awk='awk -F"\t" -v OFS="\t"'
alias les="less -S"
alias qsuball="find pbs/ -type f -name '*.pbs' | sort | xargs -I {} qsub {}"
alias rdf="readlink -f"
alias qsubi='qlogin -q interactive.q'
alias qstatall='qstat -u "*" | less'
alias qstatallrun='qstat -u "*" -s r | less'
alias qhost='qhost | less'
function qdelall {
  qstat | grep 'zhouw3' | awk -F " " '{print $1}' | xargs -I {} qdel {}
}
```

# JOB SUBMISSION AND MONITORING
See PipelineCollection github repo

# REFERENCE GENOMES FOLDER
`~/references -> /mnt/isilon/zhou_lab/projects/20191221_references`

All genome assembly is organized by their name (UCSC id if available, Ensembl id if not). 
Underneath each folder like `~/references/hg38/` you will find annotation which contains the annotation of that genome including cpg island, etc. 
Index for each software will be contained in its own folder like `~/references/hg38/biscuit`.

# USEFUL COMMAND LINE TOOLS
- [FZF](https://github.com/junegunn/fzf) - fuzzy search
- [z.sh](https://github.com/rupa/z/blob/master/z.sh) - jump around based on history

