### Apply for Digital Identity

- Go to HR service portal [https://chop.service-now.com/hrportal](https://chop.service-now.com/hrportal)
- Click "create a position"
- Fill out the job description, project description
- Non-traditional Personnel (NTP) if applicable
- Export control if applicable (OIVS, contact Adrienne Gigantino).

### Set up Remote Access

- Apply on "Service Now"
[https://chop.service-now.com/esp](https://chop.service-now.com/esp)
- search "remote access"
- fill out the form and submit

ask CHOP help desk (215-590-4357 or 4-HELP from campus) if you have question

### Be Added to HPC Storage Share

- Go to 
Cirrus
[https://ris.research.chop.edu/systems/cirrus](https://ris.research.chop.edu/systems/cirrus)
- Click
"Request Fileshare Access"
- Add yourself or people
- Choose "ReadWrite"
- add both zhou_lab and zhoulab

### HPC Connections

There are two entry nodes. You can put

```
Host hpc2
     HostName respublica-an01.research.chop.edu
     User <your_user_id>

Host hpc3
     HostName respublica-an02.research.chop.edu
     User <your_user_id>
```

in `~/.ssh/config`.

```
ssh <your_user_id>@respublica.research.chop.edu
```

For ssh key entry, follow instruction at [https://www.ssh.com/ssh/keygen](https://www.ssh.com/ssh/keygen)

```
ssh-keygen -t rsa
```

copy `~/.ssh/id_rsa.pub` to the remote's `~/.ssh/authorized_keys/`


### HPC Set up Environment

- `.bashrc` and `~/bin/`, setup the symlinks, see [HPC wiki](https://github.com/zhou-lab/labwiki/blob/master/HPC.md)
- Zhou lab pipelines, see [https://github.com/zhou-lab/labpipelines](https://github.com/zhou-lab/labpipelines)