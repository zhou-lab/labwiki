### Apply for Digital Identity (skip if you already have digital identity)

- Go to HR service portal [https://chop.service-now.com/hrportal](https://chop.service-now.com/hrportal)
- Click "create a position"
- Fill out the job description, project description
- Non-traditional Personnel (NTP) if applicable
- Export control if applicable (OIVS, contact Adrienne Gigantino).

### Apply for a PennKey

PennKey lets you download paper

https://www.research.chop.edu/services/chop-pennkey-administration

Create a bookmark for

```javascript:void(location.href="https://proxy.library.upenn.edu/login?url="+location.href)```

Click that bookmark everytime you want to access full text of a paper.

### Set up Remote Access (VPN)

- Download 'Microsoft Authenticator' app on your phone
- When you open your CHOP email for the first time it will instruct you to set up Multi-Factor Authentication
- Install 'Cisco AnyConnect' on your PC/Mac [https://www.cisco.com/c/en/us/support/security/anyconnect-secure-mobility-client/tsd-products-support-series-home.html](https://www.cisco.com/c/en/us/support/security/anyconnect-secure-mobility-client/tsd-products-support-series-home.html) or [https://software.cisco.com/download/home/286281283/type/282364313/release/4.8.02042?i=!pp](https://software.cisco.com/download/home/286281283/type/282364313/release/4.8.02042?i=!pp). Drexel student can also get a copy from the University.
- Put 'remote.chop.edu' Cisco AnyConnect, when prompted, input your CHOP ID (Email) and password. You will receive a notification on your phone to approve. Then it will connect.

ask CHOP help desk (215-590-4357 or 4-HELP from campus) if you have question

### Be Added to HPC Storage Share

- Go to 
Cirrus
[https://www.research.chop.edu/cirrus](https://www.research.chop.edu/cirrus)
- Click "Respublica Access Request"
- put your boss's activity number (ask them) and submit
- Click "Request Fileshare Access"
- Add yourself or people
- Choose "ReadWrite"
- add both zhou_lab and zhoulab (please use your fileshare if not part of Zhou lab)

### HPC connections (VDI, recommended)

[https://wiki.chop.edu/pages/viewpage.action?pageId=238785326](https://wiki.chop.edu/pages/viewpage.action?pageId=238785326)

You will receive an IP to your CHOP email address.
You should use the IP and keep the instance live for as long as possible (only disconnect).

### HPC Connections (Old way, retired)

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

### Apply for external web host

- Go to cirrus portal (https://cirrus.research.chop.edu)
- Add CostCenter Activity for VM Provisioning (if the group is not added yet) (https://www.research.chop.edu/sites/default/files/web/sites/default/files/pdfs/cirrus/Request_Access_to_CostCenter_Activity_for_VM_Provisioning.pdf)
- "Request Access to CostCenter Activity for VM Provisioning" in cirrus. 
- "Request New Virtual Machine" in cirrus, you won't see it until costcenter access is approved.(https://www.research.chop.edu/sites/default/files/web/sites/default/files/pdfs/cirrus/Request_New_Virtual_Machine.pdf)
     - Type 'Linux'
     - Find 'resInzhoupub01' in the drop-down. (At the bottom)
     - Fill in a comment like: VM Access for _Your Name_
     - Go to VM Details tab and select Admin/Sudo from the drop-down
     - Submit

ssh instruction (only work for zhoulab)

`ssh reslnzhoupub01`

HTML folder

`/data/www/html/`

Public domain name

`http://zhouserver.research.chop.edu/`


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

### General Information if needed For IT Requests/Forms

- Address: 3501 Civic Center Blvd, Philadelphia 19104
- Building: CTRB
- Department: Pathology
- Room/Cube: 9300
- Floor: 9
- Internet Provider: Other (if not known)

### Coupa Requisition (lab purchase)

The links below will assist you in updating your default settings, and provide additional Coupa job aids.

https://at.chop.edu/supply-chain/coupa

Coupa Default Settings: https://chop.service-now.com/esp?id=kb_article&sys_id=199866f0dbe6c41c3254c3d23996190b

Coupa Job Aids:  https://chop.service-now.com/esp?id=kb_article&sys_id=4001c699dbe744d0682b9b3c8a961990

Additional Coupa Training: You may login to mycareer@CHOP, click ‘Learning’, and search for Coupa Requester. Here you will be able to complete the online Coupa training modules.
