 #!/bin/bash


FILENAMES=$(cut -d , -f 1 ../STARInputList.csv)
INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment


for f in $FILENAMES
         do grep -r "concordant pair alignment rate"  ${INPDIR}/AlignSalmonGenomeBowtie/${f} | cut -d / -f 8,9
done > ../Results/summary_bowtie.txt

echo $FILENAMES | cut 


for f in $FILENAMES
         do grep -r "concordant pair alignment rate"  ${INPDIR}/AlignSalmonGenomeBowtie/${f}_10mismatch | cut -d / -f 8,9
done > ../Results/summary_bowtie_10mismatch.txt      

for f in $FILENAMES
         do grep -H "Uniquely mapped reads % |" ${INPDIR}/AlignSalmonGenomeStar/${f}Log.final.out | cut -d / -f 8
 done > ../Results/summary_star.txt

