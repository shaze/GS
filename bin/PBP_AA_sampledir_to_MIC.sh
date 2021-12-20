#!/bin/bash -l


mkdir samples/

sed -e '1,/1A-S2 Query/d' EXTRACT_1A-S2_target.fasta > temp1.fna
transeq temp1.fna temp1.faa -frame=1
echo ">Sample1" > Sample_PBP1A_AA.faa
grep -v ">" temp1.faa >> Sample_PBP1A_AA.faa


rm -f temp*
sed -e '1,/2B-S2 Query/d' EXTRACT_2B-S2_target.fasta > temp1.fna
transeq temp1.fna temp1.faa -frame=1
echo ">Sample1" > Sample_PBP2B_AA.faa
grep -v ">" temp1.faa >> Sample_PBP2B_AA.faa

rm -f temp*
sed -e '1,/2X-S2 Query/d' EXTRACT_2X-S2_target.fasta > temp1.fna
transeq temp1.fna temp1.faa -frame=1
echo ">Sample1" > Sample_PBP2X_AA.faa
grep -v ">" temp1.faa >> Sample_PBP2X_AA.faa

rm -f temp*

