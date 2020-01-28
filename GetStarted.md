### Apply for digital identity
- Non-traditional Personnel (NTP) if applicable
- Export control if applicable (OIVS, Adrienne Gigantino).
- Jim did this for me last time.

### Get remote access set up

[service now](https://chop.service-now.com/esp)

search "remote access"

fill out the form and submit

ask CHOP help desk (215-590-4357 or 4-HELP from campus) if you have question

### Be added to HPC storage share

do this on cirrhus > Request Fileshare Access

add both zhou_lab and zhoulab

### Connections

There are two entry nodes. You can put

```
Host hpc2
     HostName respublica-an01.research.chop.edu
     User zhouw3

Host hpc3
     HostName respublica-an02.research.chop.edu
     User zhouw3
```

in `~/.ssh/config`.

### Set up your home environment

- `.bashrc` and `~/bin/`, setup the symlinks, see [HPC wiki](https://github.com/zhou-lab/labwiki/blob/master/HPC.md)
- `labpipelines`, see [labpipelines](https://github.com/zhou-lab/labpipelines)