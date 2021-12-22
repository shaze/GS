# Streptococcus surveillance pipelines


This is a work in progress of a Nextflow implementation of `https://github.com/BenJamesMetcalf/GBS_Scripts_Reference` and related piplines. The bioinformatics has not changed only the wrapping. Any credit due to Ben and Sopio -- errors are mine.


# 1. Installing


## 1.1 Essential pre-requirements

To run this workflow you must have Java 8 or later and Nextflow installed. An example of doing this is

```
sudo apt install openjdk-11-jre-headless
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
sudo chmod a+rx /usr/local/bin/nextflow
```

## 1.2 Bioinformatics software requirements

There is an extensive set of software requirements for this. You can either install manually or we suggest using our containerised image with all the necessary requirements. In this case you need only install either `docker` **or** `singularity` and use the `-profile docker` or `-profile singularity` option. If you can't do this or want to install your own versions you can use the `dockers/Dockerfile` file as a guide of what to install.

To install Docker or Singularity requires root privileges. Generally, Docker is not available or recommended in shared computing environments like HPC clusters but Singularity is widely available: you do not have to be root in order to _run_ Singularity.

### 1.2.1 If you choose Singularity

You need to run Singularity 3 for this workflow.
* A package is not available on Ubuntu and you have to install yourself. Instructions can be found here: https://sylabs.io/guides/3.0/user-guide/installation.html on how to do this (look for `Install the Debian/Ubuntu package using apt`). The instructions are clear and the steps are not onerous but if you have a single user machine and are root it may be easier to use Docker.
* On RHEL versions, you can find it in the OSG repo. You need to add this repo to `/etc/yum.repos.d`

### 1.2.2 If you choose Docker

Remember Docker is not meant for shared computer systems because of security constraints. It's a good option if you are the only user of your computer.

*  on Ubuntu install `docker.io`.
*  On RHEL, the package is `docker`.
*  On MacOS, you can install from here: https://docs.docker.com/desktop/mac/install/

On Linux platforms, you also need to add your user to the _docker_ group, e.g.
`sudo usermod -a -G docker scott` so that you can run Docker as a normal user. This isn'tneed on the Mac. 



## 1.3 Installing the workflow itself

This takes one line

`nextflow pull shaze/GS`

Note that if the workflow is updated when you run the workflow you will get a message like "NOTE: Your local project version looks outdated - a different revision is available in the remote repository". If you want to upgrade to the latest version do another `nextflow pull shaze/GS` 


## 1.4 Getting Reference Databases

The workflows require some reference databases (about 400MB in total). There are three databases (each placed in a separate directory) 

The easiest way of getting them is to say

`nexflow run shaze/GS/fetchDB.nf --target XXX`

Replace XXX with the name of the directory which you want to store the databases. This directory can exist already but if it doesn't exist will be created.  For example, if I say `nexflow run shaze/GS/fetchDB.nf --target /data/strep` then (whether or not that directory previously exists), then the three databases will be created and we will not have three directories/databases, viz.,  `/data/strep/GAS_ReferenceDB`, `/data/strep/GBS_ReferenceDB` and `/data/strep/SPN_Reference_DB`. You can use this for the workflow parameters as described below.


# 2. Instructions for running


## 2.1 Workflow parameters
You need the following three parameters (defaults are given in square brackets -- if you omit the parameter the default is used)

* `batch_dir` : the name of the directory containing the read pair files. There must be
  exactly two files per sample, named appropriately [`/dataC/CRDM/testingreads_gbs/191206_M02143`]
* `out_dir`: the name of the output directory where data should go ["output"]
* `strepA_DB` _or_ `strepB_DB` _or_ `SPN_DB`: the name of the database directory [`/dataC/CRDM/G[AB]_Reference_DB`]



*Important assumption* This workflow assumes that all input fastq files end with `L001_R1_001` or `L001_R2_001`, followed by a suffix (this assumption is required by srst2). The default suffix is `fastq.gz`

See advanced parameters below.


## 2.2 Additional Nextflow parameters

In addition there are two _Nextflow_ parameters that you can use (especially the first). Not that for the above parameters you used two dashes `--` while for the ones below you use only one dash `-`: this is very important.

* `-profile slurm`: this causes each job to be submitted to the job scheduler and improves your parallelism. This is highly desirable.  The workflow provides a high degree of parallelism -- using `max_forks` of 37 with 37 input pairs, the workflow took  20 minutes of elapsed time (25 CPU hours running running several hundred jobs).  If you do *not* use this option, Nextflow executes this job on the computer that you happen to be running on.  Nextflow detects the number of actual cores you have and parallelises sensibly as much as possible (a minumum number of four cores is required).

* `-profile docker`: Run the workflow using Docker containers (see above). Docker must be installed on your system.

* `-profile sigularity`: Run the workflow using Singularity containers (see above). Singularity must be installed on your system.

* NB: _if_ you want to combine profiles, e.g. `slurm` and `docker`, say `-profile slurm,docker`.  

* `-resume`: if something crashed in a run for a reason other than an error in the workflow (say the computer's power went down) then you can use this to pick up execution from the point that execution failed (obviously if there's a bug in the workflow or a problem with the data the workflow will just crash again)



# 3. Running the workflows


## 3.1. `GAS_Scripts_Reference` workflow


```

nextflow run shaze/GS/strepA.nf  --batch_dir NAMEOFINPUTDIR --out_dir NAMEOFOUTPUTDIR --strepA_DB DBdirectory --max_forks  NUMBER
                        

```


A typical run might be

```

nextflow run shaze/GS/strepA.nf  --batch_dir august_data --out_dir /dataC/archive/august -profile slurm
                        

```

This takes that data from `august_data` and puts the results in `/dataC/archive/august` using the default parameters for the database and `max_forks`. It also says that jobs should be submitted to the cluster by the SLURM scheduler which allows all the computers in the cluster to be used.

You can run this from any directory. *If you are running on the Wits cluster please do not run from your home directory -- the output and intermediate files are very large and put pressure on our backup. Use your given project directory*.

This run above assumes that the database directory is in the default place. On the Wits cluster this is `/dataC/CRDM/GAS_Reference_DB`. If it is somehwere else use the `--strepA_DB` argument -- for example `--strepA_DB /home/scott/DBS/GAS_A`.

## 3.2. Running the `GBS_Scripts_Reference` workflow



```

nextflow run shaze/GS/strepB.nf  --batch_dir NAMEOFINPUTDIR --out_dir NAMEOFOUTPUTDIR --strepB_DB DBdirectory --max_forks  NUMBER
                        

```

A typical run might be

```

nextflow run shaze/GS/strepB.nf  --batch_dir august_data --out_dir /dataC/archive/august -profile slurm
                        

```

The same remarks as for `strepA` apply for the various parameters except that if the reference DB is not in the default place (`/dataC/CRDM/GBS_Reference_DB`) then use the `--strepB_DB` option. 


## 3.3. Running the `SPN_Scripts_Reference` workflow



```

nextflow run shaze/GS/spn.nf  --batch_dir NAMEOFINPUTDIR --out_dir NAMEOFOUTPUTDIR --SPN_DB DBdirectory --max_forks  NUMBER
                        

```

A typical run might be

```

nextflow run shaze/GS/strepB.nf  --batch_dir august_data --out_dir /dataC/archive/august -profile slurm
                        

```

The same remarks as for `strepA` apply for the various parameters except that if the reference DB is not in the default place (`/dataC/CRDM/GBS_Reference_DB`) then use the `--strepB_DB` option. 


# 5. Output

The output directory contains
* an Excel spreadsheet summarising the results
* MultiQC reports
* a directory for each sample contains a QC report and the Velvet contigs
* a directory `new_mlst` with new MLST data


# 6. Cleaning up

This script has lots of intermediate files and output files. As an example, a test data aset of 7GB led to 12GB of output data and another 43GB of intermediate files.  Nextflow uses the `work` directory to store all the intermediate files and this needs to be cleaned out but some care needs to be taken.
* delete the work directory `/bin/rm -rf work`  (there are fancier ways of doing this but this works).




# 7.  Advanced paramters

* `max_forks`: a parameter limiting the parallelism. Nextflow understands your environment (whether you are running on a stand-alone computer or on a cluster) and handles the parallelism appropriately. If you have an 8-core machine it can run up to 8 processes at the same time for you; if you have a cluster with a scheduler it will interact with the scheduler to run many jobs in parallel.  However, there are times when you need to tell Nextflow to _restrict_ the amount of parallelism. In this pipeline's case, depending on your disk system, disk performance may be a reason to do so. You may have 1000 CPUs in the cluster but if if they each try to read a file on the same disk, the disk will be very, very slow and you will become very unpopular. 
With `max_forks`, only this number of samples are allowed to be in the first phase of the workflow at the same time and acts as a throttle (there will in general be more than this number of processes happening in parallel as this throttles only the number of samples being processed in the first phase not the total work being done. The default is 10 --  the reason for this throttle is not so much to limit the number of jobs running (since you can rely on the scheduler to do this sensibly) but that although these files are not huge if a lot of work is being done in parallel the I/O performance can suffer. Experiment so that you can get things done quickly without making everyone else hate you.


* `suffix`: default is `fastq.gz`. The workflow assumes that all the fastq files end with this suffix. If you have a different suffix, then use this parameter, for example `--suffix fq.gz`. All input files must have the same suffix.

* `max_velvet_cpus`: how many CPUs should be used for Velvet. While Velvet has some parts that are parallelisable, a simple application of Amdahl's law shows that it is not worth allocating many threads to individual runs of Velvet. Since we are running Velvet on many different data sets, it makes more sense to run more of these independent jobs in parallel.

# 8. Creating a config file

You can create a config file like this (say `demo.config`) and then run your code 

```
params {
      	   out_dir = "output"
	   batch_dir="/dataC/CRDM/testingreads_gbs/191206_M02143"
	   allDB_dir="/dataC/CRDM/GBS_Scripts_Reference/GBS_Reference_DB"
	   max_forks=10		   
}
```

with the the line `nextflow run -c demo.config ... ` 



