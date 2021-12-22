
params.src="http://www.bioinf.wits.ac.za/data"


params.target="./"


gas=file(params.src+"/GAS_Reference_DB.zip")
gbs=file(params.src+"/GBS_Reference_DB.zip")
spn=file(params.src+"/SPN_Reference_DB.zip")



process fetch {
   input:
      file(gas)
      file(gbs)
      file(spn)
   output:
      file("*DB") into results
   publishDir params.target, mode:'move'
   script:
     """
     unzip GAS_Reference_DB.zip
     unzip GBS_Reference_DB.zip
     unzip SPN_Reference_DB.zip
     /bin/rm -f `readlink $gas`
     /bin/rm -f `readlink $gbs`
     /bin/rm -f `readlink $spn`          
     """
}
