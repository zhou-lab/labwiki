### Apply for Digital Identity (skip if you already have digital identity)

- Go to HR service portal [https://chop.service-now.com/hrportal](https://chop.service-now.com/hrportal)
- Click "create a position"
- Fill out the job description, project description
- Non-traditional Personnel (NTP) if applicable
- Export control if applicable (OIVS, contact Adrienne Gigantino).

### Set up Remote Access

- Apply on "Service Now"
[https://chop.service-now.com/esp](https://chop.service-now.com/esp)
- search "remote access"
- fill out the form and submit (take about a week to get approved)
     - Address: 3501 Civic Center Blvd, Philadelphia 19104
     - Building: CTRB
     - Department: Pathology
     - Room/Cube: 9300
     - Floor: 9
     - Business justification: Access HPC from remote location
     - Internet Provider: Other (if not known)
- download 'Entrust IdentityGuard Mobile' on your phone, and set up following [https://ishelp.chop.edu/](https://ishelp.chop.edu/)
- install 'Cisco AnyConnect' on your PC/Mac [https://www.cisco.com/c/en/us/support/security/anyconnect-secure-mobility-client/tsd-products-support-series-home.html](https://www.cisco.com/c/en/us/support/security/anyconnect-secure-mobility-client/tsd-products-support-series-home.html) or [https://software.cisco.com/download/home/286281283/type/282364313/release/4.8.02042?i=!pp](https://software.cisco.com/download/home/286281283/type/282364313/release/4.8.02042?i=!pp). Drexel student can also get a copy from the University.
- Put 'remote.chop.edu' Cisco AnyConnect, when prompted, put your CHOP ID on the first line, what you get from 'Entrust' for Passcode on the second line, and your CHOP password for the 2nd passcode on the third.

ask CHOP help desk (215-590-4357 or 4-HELP from campus) if you have question

### Be Added to HPC Storage Share

- Go to 
Cirrus
[https://ris.research.chop.edu/systems/cirrus](https://ris.research.chop.edu/systems/cirrus)
- Click "Respublica Access Request"
- put cost center (22105) and activity (7229660000) submit
- Click "Request Fileshare Access"
- Add yourself or people
- Choose "ReadWrite"
- add both zhou_lab and zhoulab

### HPC Connections

```
ssh <CHOPID>@respublica.research.chop.edu
```


`respublica.research.chop.edu` randomly assigns you one of the two login nodes (`respublica-an01.research.chop.edu` and `respublica-an02.research.chop.edu`)


### .ssh configuration

You can put

```
Host hpc2
     HostName respublica-an01.research.chop.edu
     User <CHOPID>

Host hpc3
     HostName respublica-an02.research.chop.edu
     User <CHOPID>
```

in `~/.ssh/config` and use `ssh hpc2` and `ssh hpc3` to connect.

### .ssh password-less entry

For ssh key entry, follow instruction at [https://www.ssh.com/ssh/keygen](https://www.ssh.com/ssh/keygen)

```
ssh-keygen -t rsa
```

Append the content of `~/.ssh/id_rsa.pub` to the remote's `~/.ssh/authorized_keys`

### Apply for external web access

- Go to cirrus portal (https://cirrus.research.chop.edu)
- Add CostCenter Activity for VM Provisioning (if the group is not added yet)
- "Request Access to CostCenter Activity for VM Provisioning" in cirrus
- "Request New Virtual Machine" in cirrus, you won't see it until costcenter access is approved.

### HPC Set up Environment

Put
```
. /mnt/isilon/zhoulab/labtools/bashrc/chop/bashrc_hpc_zhoulab
```

to `~/.bashrc` (customize if you know how and feel like so).

Setup the symlinks, see [HPC wiki](https://github.com/zhou-lab/labwiki/blob/master/HPC.md)

### Web Citrix Workspace Remote Access (an inconvenient alternative)

- Go to [https://connect.chop.edu](connect.chop.edu)
- Put your user name, password and Entrust token (see above for how to get that)
- click on Desktop and you should be able to access CHOP server using terminal.
