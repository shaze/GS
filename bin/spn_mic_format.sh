#!/bin/bash
fout=$1
fin=$2
awk 'NR<=2'  $fout | awk -F"," '{$1=""; print $0}'  | sed 's/"//g' | sed 's/ //' > temp1.txt
echo "WGS_AMP_SIGN WGS_AMP WGS_AMP_SIR \
WGS_CPT_SIGN WGS_CPT WGS_CPT_SIR \
WGS_ZOX_SIGN WGS_ZOX WGS_ZOX_SIR \
WGS_FOX_SIGN WGS_FOX WGS_FOX_SIR" > temp2.txt

echo "NA NA NA \
NA NA NA \
NA NA NA \
NA NA NA " >> temp2.txt
paste -d ' ' temp1.txt temp2.txt> $fin
