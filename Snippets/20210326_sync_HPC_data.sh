# usage:
# sd hpc5:~/a/b/c # sync from remote to local
# sd ~/a/b/c # sync from local to remote

function sd() {
  local_home="/Users/zhouw3"
  remote_home="/home/zhouw3"
  hpc_name="hpc5"

  from=$1
  if [[ $from =~ ^$local_home ]]; then # from local to remote
    to=${from/$local_home/$hpc_name":"$remote_home}
  elif [[ $from =~ ^$remote_home ]]; then # from remote to local
    from=$hpc_name":"$from
    to=${from/$remote_home/$local_home}
  elif [[ $from =~ ^$hpc_name ]]; then # from remote to local
    from=${from/\~/$remote_home}
    to=${from/$hpc_name":"/}
    to=${to/#\~/$local_home}
    to=${to/$remote_home/$local_home}
  elif [[ $from =~ ^~ ]]; then # from local to remote
    from=${from/\~/$local_home}
    to=${to/#\~/$remote_home}
    to=$hpc_name":"$from
  else
    exit 1
  fi

  echo "From: "$from
  echo "To:   "$to

  if [[ $from =~ ^$local_home ]]; then
    ssh $hpc_name mkdir -p $(dirname ${to/$hpc_name":"/})
    scp $from $to
  elif [[ $from =~ ^$hpc_name ]]; then
    mkdir -p $(dirname $to)
    scp $from $to
  fi
}
