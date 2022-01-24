#!/usr/bin/env python3


#  (C) University of the Witwatersrand, Johannesburg 2021
#  Scott Hazelhurst
#  MIT License as specified in https://github.com/shaze/GS/blob/master/LICENSE.md

import re
import pandas as pd
import glob
import sys
import os.path
from typing import List

report    = sys.argv[1]
suffix    = sys.argv[2]
batch_dir = sys.argv[3]
out_dir   = sys.argv[4]
outf:str=sys.argv[5]
csvf=outf.replace("xlsx","tsv")
low_cov   = sys.argv[6].split(",")

if report not in ["SPN","A","B"]:
    sys.exit("Unknown report <%s>"%report)


all_heads={}

all_heads['A']="Sample,emm_Type,ST,gki,gtr,murI,mutS,recP,xpt,yqiL,T_Type,Group_A,EMM_Family,Other_Surface_Proteins,Capsule,SDA1,SLAA,SIC,ROCA,PNGA3,NADase_D330G,Exotoxins,PBP_ID,WGS_ZOX_SIGN,WGS_ZOX,WGS_ZOX_SIR,WGS_FOX_SIGN,WGS_FOX,WGS_FOX_SIR,WGS_TAX_SIGN,WGS_TAX,WGS_TAX_SIR,WGS_CFT_SIGN,WGS_CFT,WGS_CFT_SIR,WGS_CPT_SIGN,WGS_CPT,WGS_CPT_SIR,WGS_CZL_SIGN,WGS_CZL,WGS_CZL_SIR,WGS_AMP_SIGN,WGS_AMP,WGS_AMP_SIR,WGS_PEN_SIGN,WGS_PEN,WGS_PEN_SIR,WGS_MER_SIGN,WGS_MER,WGS_MER_SIR,ER_CL,WGS_ERY_SIGN,WGS_ERY,WGS_ERY_SIR,WGS_CLI_SIGN,WGS_CLI,WGS_CLI_SIR,WGS_LZO_SIGN,WGS_LZO,WGS_LZO_SIR,WGS_SYN_SIGN,WGS_SYN,WGS_SYN_SIR,WGS_ERY/CLI,TET,WGS_TET_SIGN,WGS_TET,WGS_TET_SIR,GYRA_PARC,WGS_CIP_SIGN,WGS_CIP,WGS_CIP_SIR,WGS_LFX_SIGN,WGS_LFX,WGS_LFX_SIR,OTHER,WGS_DAP_SIGN,WGS_DAP,WGS_DAP_SIR,WGS_VAN_SIGN,WGS_VAN,WGS_VAN_SIR,WGS_RIF_SIGN,WGS_RIF,WGS_RIF_SIR,WGS_CHL_SIGN,WGS_CHL,WGS_CHL_SIR,WGS_SXT_SIGN,WGS_SXT,WGS_SXT_SIR,Contig_Num,N50,Longest_Contig,Total_Bases,ReadPair_1,Contig_Pathn".split(",")
all_heads['B']="Sample,Serotype,ST,adhP,pheS,atr,glnA,sdhA,glcK,tkt,PBP_1A,PBP_2B,PBP_2X,WGS_ZOX_SIGN,WGS_ZOX,WGS_ZOX_SIR,WGS_FOX_SIGN,WGS_FOX,WGS_FOX_SIR,WGS_TAX_SIGN,WGS_TAX,WGS_TAX_SIR,WGS_CFT_SIGN,WGS_CFT,WGS_CFT_SIR,WGS_CPT_SIGN,WGS_CPT,WGS_CPT_SIR,WGS_CZL_SIGN,WGS_CZL,WGS_CZL_SIR,WGS_AMP_SIGN,WGS_AMP,WGS_AMP_SIR,WGS_PEN_SIGN,WGS_PEN,WGS_PEN_SIR,WGS_MER_SIGN,WGS_MER,WGS_MER_SIR,TET,WGS_TET_SIGN,WGS_TET,WGS_TET_SIR,EC,WGS_ERY_SIGN,WGS_ERY,WGS_ERY_SIR,WGS_CLI_SIGN,WGS_CLI,WGS_CLI_SIR,WGS_LZO_SIGN,WGS_LZO,WGS_LZO_SIR,WGS_SYN_SIGN,WGS_SYN,WGS_SYN_SIR,WGS_ERYCLI,FQ,WGS_CIP_SIGN,WGS_CIP,WGS_CIP_SIR,WGS_LFX_SIGN,WGS_LFX,WGS_LFX_SIR,Other,WGS_DAP_SIGN,WGS_DAP,WGS_DAP_SIR,WGS_VAN_SIGN,WGS_VAN,WGS_VAN_SIR,WGS_RIF_SIGN,WGS_RIF,WGS_RIF_SIR,WGS_CHL_SIGN,WGS_CHL,WGS_CHL_SIR,WGS_SXT_SIGN,WGS_SXT,WGS_SXT_SIR,ALPH,SRR,Pili,HVGA,Contig_num,N50,Longest_contig,Total_bases,ReadPair_1,Contig_path".split(",")
all_heads['SPN']="Sample,WGS_Serotype,Pili,ST,aroe,gdh,gki,recP,spi,xpt,ddl,PBP1A,PBP2B,PBP2X,WGS_PEN_SIGN,WGS_PEN,WGS_PEN_SIR_Meningitis,WGS_PEN_SIR_Nonmeningitis,WGS_AMO_SIGN,WGS_AMO,WGS_AMO_SIR,WGS_MER_SIGN,WGS_MER,WGS_MER_SIR,WGS_TAX_SIGN,WGS_TAX,WGS_TAX_SIR_Meningitis,WGS_TAX_SIR_Nonmeningitis,WGS_CFT_SIGN,WGS_CFT,WGS_CFT_SIR_Meningitis,WGS_CFT_SIR_Nonmeningitis,WGS_CFX_SIGN,WGS_CFX,WGS_CFX_SIR,WGS_AMP_SIGN,WGS_AMP,WGS_AMP_SIR,WGS_CPT_SIGN,WGS_CPT,WGS_CPT_SIR,WGS_ZOX_SIGN,WGS_ZOX,WGS_ZOX_SIR,WGS_FOX_SIGN,WGS_FOX,WGS_FOX_SIR,EC,WGS_ERY_SIGN,WGS_ERY,WGS_ERY_SIR,WGS_CLI_SIGN,WGS_CLI,WGS_CLI_SIR,WGS_SYN_SIGN,WGS_SYN,WGS_SYN_SIR,WGS_LZO_SIGN,WGS_LZO,WGS_LZO_SIR,WGS_ERY/CLI,Cot,WGS_COT_SIGN,WGS_COT,WGS_COT_SIR,Tet,WGS_TET_SIGN,WGS_TET,WGS_TET_SIR,WGS_DOX_SIGN,WGS_DOX,WGS_DOX_SIR,FQ,WGS_CIP_SIGN,WGS_CIP,WGS_CIP_SIR,WGS_LFX_SIGN,WGS_LFX,WGS_LFX_SIR,Other,WGS_CHL_SIGN,WGS_CHL,WGS_CHL_SIR,WGS_RIF_SIGN,WGS_RIF,WGS_RIF_SIR,WGS_VAN_SIGN,WGS_VAN,WGS_VAN_SIR,WGS_DAP_SIGN,WGS_DAP,WGS_DAP_SIR,Contig_Num,N50,Longest_Contig,Total_Bases,ReadPair_1,Contig_Path".split(",")

abh="Sample,MLST,Serotype,PBP1A,PBP2B,PBP2X,HVGA,PI1,PI2A1,PI2A2,PI2B,SRR1,SRR2,ALP1REF,ALP23REF,ALPHAREF,RIBREF".split(",")
all_bin_heads={}
all_bin_heads['A']=abh+["-" for i in range(len(abh),49)]
all_bin_heads['B']="Sample,MLST,Serotype,PBP1A,PBP2B,PBP2X,HVGA,PI1,PI2A1,PI2A2,PI2B,SRR1,SRR2,ALP1REF,ALP23REF,ALPHAREF,RIBREF".split(",")
all_bin_heads['SPN']=[]

table_heads=heads=all_heads[report]
    
bin_heads=all_bin_heads[report]

    


print(table_heads)

tables = glob.glob("*_TABLE_*")
table_data=[f.readline().split() for f in map(open, tables)]
for i, t in enumerate(table_data):
    curr_line:List[str] = table_data[i]
    sample=curr_line[0]
    fq = os.path.join(batch_dir,sample)+"_R1_001."+suffix
    if report=="A":
        curr_line.append(fq)
        curr_line.append("%s/%s/velvet_output/contigs.fa"%(out_dir,sample))
    else:    
        curr_line.append("%s/%s/velvet_output/contigs.fa"%(out_dir,sample))
        curr_line[-2]=fq
    # makes assumptions about the input format

table = pd.DataFrame(table_data,columns=table_heads,dtype='string')
table.sort_values(by="Sample",inplace=True)
table.to_csv(csvf,sep="\t")


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


low_f = pd.DataFrame(low_cov,columns=["Sample"])



with pd.ExcelWriter(outf,engine_kwargs={'options': {'strings_to_numbers': True}}) as writer: 
    table.to_excel(writer, "TABLE")
    if bin_heads:
        bins = glob.glob("*_BIN_*")
        bin_data=[f.readline().split(",") for f in map(open, bins)]
        bin=pd.DataFrame(bin_data,columns=bin_heads)
        bin.to_excel(writer, "BIN")
    pdf.to_excel(writer,"PBP")
    mdf.to_excel(writer,"MLST")
    low_f.to_excel(writer,"Bad samples")
