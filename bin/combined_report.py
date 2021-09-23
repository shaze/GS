#!/usr/bin/env python3

import re
import pandas as pd
import glob
import sys
import os.path
from typing import List

batch_dir = os.path.realpath(sys.argv[1])
out_dir   = os.path.realpath(sys.argv[2])

table_heads="Sample,Serotype,ST,adhP,pheS,atr,glnA,sdhA,glcK,tkt,PBP_1A,PBP_2B,PBP_2X,WGS_ZOX_SIGN,WGS_ZOX,WGS_ZOX_SIR,WGS_FOX_SIGN,WGS_FOX,WGS_FOX_SIR,WGS_TAX_SIGN,WGS_TAX,WGS_TAX_SIR,WGS_CFT_SIGN,WGS_CFT,WGS_CFT_SIR,WGS_CPT_SIGN,WGS_CPT,WGS_CPT_SIR,WGS_CZL_SIGN,WGS_CZL,WGS_CZL_SIR,WGS_AMP_SIGN,WGS_AMP,WGS_AMP_SIR,WGS_PEN_SIGN,WGS_PEN,WGS_PEN_SIR,WGS_MER_SIGN,WGS_MER,WGS_MER_SIR,TET,WGS_TET_SIGN,WGS_TET,WGS_TET_SIR,EC,WGS_ERY_SIGN,WGS_ERY,WGS_ERY_SIR,WGS_CLI_SIGN,WGS_CLI,WGS_CLI_SIR,WGS_LZO_SIGN,WGS_LZO,WGS_LZO_SIR,WGS_SYN_SIGN,WGS_SYN,WGS_SYN_SIR,WGS_ERYCLI,FQ,WGS_CIP_SIGN,WGS_CIP,WGS_CIP_SIR,WGS_LFX_SIGN,WGS_LFX,WGS_LFX_SIR,Other,WGS_DAP_SIGN,WGS_DAP,WGS_DAP_SIR,WGS_VAN_SIGN,WGS_VAN,WGS_VAN_SIR,WGS_RIF_SIGN,WGS_RIF,WGS_RIF_SIR,WGS_CHL_SIGN,WGS_CHL,WGS_CHL_SIR,WGS_SXT_SIGN,WGS_SXT,WGS_SXT_SIR,ALPH,SRR,Pili,HVGA,Contig_num,N50,Longest_contig,Total_bases,ReadPair_1,Contig_path".split(",")

bin_heads="Sample,MLST,Serotype,PBP1A,PBP2B,PBP2X,HVGA,PI1,PI2A1,PI2A2,PI2B,SRR1,SRR2,ALP1REF,ALP23REF,ALPHAREF,RIBREF".split(",")



tables = glob.glob("*_TABLE_*")
table_data=[f.readline().split() for f in map(open, tables)]
for i, t in enumerate(table_data):
    curr_line:List[str] = table_data[i]
    sample=curr_line[0]
    curr_line.append("%s/%s/velvet_output/contigs.fa"%(out_dir,sample))
    curr_line[-2]=glob.glob(os.path.join(batch_dir,sample)+"*R1*gz")[0]
    
table = pd.DataFrame(table_data,columns=table_heads)


bins = glob.glob("*_BIN_*")
bin_data=[f.readline().split(",") for f in map(open, bins)]
bin=pd.DataFrame(bin_data,columns=bin_heads)

mlst=glob.glob("*new_mlst*")
if len(mlst)==0: mlst=["---None---"]
mdf =pd.DataFrame(mlst)


plist=[]
pbp = glob.glob("*_newPBP_allele_info.txt")
for p in pbp:
    data=open(p).readline().strip().split()
    plist.append([data[0],data[-1]])
if len(plist)==0: plist=[["---None---", "---None---"]]
pdf = pd.DataFrame(plist,columns=["Sequence","Allele"])

outf:str=sys.argv[3]
with pd.ExcelWriter(outf,engine_kwargs={'strings_to_numbers': True}) as writer: 
    table.to_excel(writer, "TABLE")
    bin.to_excel(writer, "BIN")
    pdf.to_excel(writer,"PBP")
    mdf.to_excel(writer,"MLST")
