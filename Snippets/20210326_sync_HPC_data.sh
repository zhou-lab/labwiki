# USAGE

## infer direction by folder names

## remote > local
# sd /mnt/isilon/a/b/c
# sd /home/zhouw3/a/b/c
# sd hpc:~/a/b/c

## local > remote
# sd /Users/zhouw3/a/b/c
# sd ~/a/b/c

### INSTALLATION
### insert the following to your .zshrc/.bashrc
### ----------------
# LOCAL_HOME="/Users/zhouw3"
# REMOTE_HOME="/home/zhouw3"
# HPC_NAME="hpc5"
# source <(curl -s https://raw.githubusercontent.com/zhou-lab/labwiki/master/Snippets/20210326_sync_HPC_data.sh)
### ----------------

function sd() {
  [[ -z "$LOCAL_HOME" ]] && LOCAL_HOME="/Users/zhouw3"
  [[ -z "$REMOTE_HOME" ]] && REMOTE_HOME="/home/zhouw3"
  [[ -z "$REMOTE_HOME2" ]] && REMOTE_HOME2="/mnt/isilon"
  [[ -z "$HPC_NAME" ]] && HPC_NAME="hpc"

  from=$1
  if [[ $from =~ ^$LOCAL_HOME ]]; then # from local to remote
    to=${from/$LOCAL_HOME/$HPC_NAME":"$REMOTE_HOME}
  elif [[ $from =~ ^$REMOTE_HOME2 ]]; then # from remote to local
    from=$HPC_NAME":"$from
    to={$from/$REMOTE_HOME2/$LOCAL_HOME}
  elif [[ $from =~ ^$REMOTE_HOME ]]; then # from remote to local
    from=$HPC_NAME":"$from
    to=${from/$REMOTE_HOME/$LOCAL_HOME}
  elif [[ $from =~ ^$HPC_NAME ]]; then # from remote to local
    from=${from/\~/$REMOTE_HOME}
    to=${from/$HPC_NAME":"/}
    to=${to/#\~/$LOCAL_HOME}
    to=${to/$REMOTE_HOME/$LOCAL_HOME}
  elif [[ $from =~ ^~ ]]; then # from local to remote
    from=${from/\~/$LOCAL_HOME}
    to=${to/#\~/$REMOTE_HOME}
    to=$HPC_NAME":"$from
  else
    return 1;
  fi

  echo "From: "$from
  echo "To:   "$to

  if [[ $from =~ ^$LOCAL_HOME ]]; then
    ssh $HPC_NAME mkdir -p $(dirname ${to/$hpc_name":"/})
    scp $from $to
  elif [[ $from =~ ^$HPC_NAME ]]; then
    mkdir -p $(dirname $to)
    scp $from $to
  fi
}
