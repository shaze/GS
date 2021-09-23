# strepB.nf

This is a work in progress of a Nextflow implementation of `https://github.com/BenJamesMetcalf/GBS_Scripts_Reference`. The bioinformatics has not changed only the wrapping. Any credit due to Ben and Sopio -- errors are mine.


# Installing

Run 

`nextflow pull shaze/GS`

Note that if the workflow is updated when you run the workflow you will get a message like "NOTE: Your local project version looks outdated - a different revision is available in the remote repository". If you want to upgrade to the latest version do another `nextflow pull shaze/GS` 

# Instructions for running


## Workflow parameters
You need the following three parameters (defaults are given in square brackets -- if you omit the parameter the default is used)

* `batch_dir` : the name of the directory containing the read pair files. There must be
  exactly two files per sample, named appropriately [`/dataC/CRDM/testingreads_gbs/191206_M02143`]
* `out_dir`: the name of the output directory where data should go ["output"]
* `allDB_dir`: the name of the database directory [`/dataC/CRDM/GBS_Scripts_Reference/GBS_Reference_DB`]
* `max_forks`: a parameter limiting the parallelism. Essentially only this number of samples are allowed to be in the first phase of the workflow at the same time and acts as a throttle (there will in general be more than this number of processes happening in parallel as this throttles only the number of samples being processed in the first phase not the total work being done. The default is 10 --  the reason for this throttle is not so much to limit the number of jobs running (since you can rely on the scheduler to do this sensibly) but that although these files are not huge if a lot of work is being done in parallel the I/O performance can suffer. Experiment so that you can get things done quickly without making everyone else hate you.




## Additional Nextflow parameters

In addition there are two _Nextflow_ parameters that you can use (especially the first). Not that for the above parameters you used two dashes `--` while for the ones below you use only one dash `-`: this is very important.

* `-profile slurm`: this causes each job to be submitted to the job scheduler and improves your parallelism. This is highly desirable.  The workflow provides a high degree of parallelism -- using `max_forks` of 37 with 37 input pairs, the workflow took  23 minutes of elapsed time (32 CPU hours running)

* `-profile resume`: if something crashed in a run for a reason other than an error in the workflow (say the computer's power went down) then you can use this to pick up execution from the point that execution failed (obviously if there's a bug in the workflow or a problem with the data the workflow will just crash again)


# Running the `GBS_Scripts_Reference` workflow



```

nextflow run shaze/GS/strepB.nf  --batch_dir NAMEOFINPUTDIR --out_dir NAMEOFOUTPUTDIR --allDB_dir DBdirectory --max_forks  NUMBER
                        

```

A typical run might be

```

nextflow run shaze/GS/strepB.nf  --batch_dir august_data --out_dir /dataC/archive/august
                        

```

This takes that data from `august_data` and puts the results in `/dataC/archive/august` using the default parameters for the database and `max_forks`


## Output

The output directory contains
* an Excel spreadsheet summarising the results
* MultiQC reports
* a directory for each sample contains a QC report and the Velvet contigs
* a directory `new_mlst` with new MLST data



## Creating a config file

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



