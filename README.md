#strep.nf

This is a work in progress of a Nextflow implementation of https://github.com/BenJamesMetcalf/GBS_Scripts_Reference


# Instructions

You need the following three parameters (defaults are given in square brackets -- if you omit the parameter the default is used)

* `batch_dir` : the name of the directory containing the read pair files. There must be
  exactly two files per sample, named appropriately ["/dataC/CRDM/testingreads_gbs/191206_M02143"]
* `out_dir`: the name of the output directory where data should go ["output"]
* `allDB_dir`: the name of the database directory
* `max_forks`: a parameter limiting the parallelism. Essentially only this number of samples are allowed to be in the first phase of the workflow at the same time and acts as a throttle (there will in general be more than this number of processes happening in parallel as this throttles only the number of samples being processed in the first phase not the total work being done. The default is 10 --  the reason for this throttle is not so much to limit the number of jobs running (since you can rely on the scheduler to do this sensibly) but that although these files are not huge if a lot of work is being done in parallel the I/O perforamnce can suffer.


```

nextflow run strepB.nf  --batch_dir NAMEOFINPUTDIR --out_dir NAMEOFOUTPUTDIR --allDB_dir DBdirectory --max_forks  NUMBER
                        

```

You do not have to be in the same directory as the the Nextflow script (and in fact it's bad practice to do so) to run. Just use the 


In addition there are two _Nextflow_ parameters that you can use (especially the first). Not that for the above parameters you used two dashes `--` while for the ones below you use only one dash `-`: this is very important.

* `-profile slurm`: this causes each job to be submitted to the job scheduler and improves your parallelism. This is highly desirable.

* `-profile resume`: if something crashed in a run for a reason other than an error in the workflow (say the computer's power went down) then you can use this to pick up execution from the point that execution failed (obviously if there's a bug in the workflow or a problem with the data the workflow will just crash again)


# Creating a config file

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



