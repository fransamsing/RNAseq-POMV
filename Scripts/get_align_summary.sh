 #!/bin/bash


FILENAMES=$(cut -d , -f 1 ../STARInputList.csv)
INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment


for f in $FILENAMES
         do grep -r "overall read mapping rate"  ${INPDIR}/AlignSalmonGenomeBowtie/TophatDefaults/${f} | cut -d / -f 9,10
done > ../Results/summary_tophat_defaults.txt


for f in $FILENAMES
         do grep -r "overall read mapping rate"  ${INPDIR}/AlignSalmonPOMVCombined/${f} | cut -d / -f 8,9
done > ../Results/summary_tophat_salmon_POMV.txt

for f in $FILENAMES
         do grep -r "overall read mapping rate"  ${INPDIR}/AlignSalmonISAVCombined/${f} | cut -d / -f 8,9
done > ../Results/summary_tophat_salmon_ISAV.txt     

for f in $FILENAMES
         do grep -Hi "Uniquely mapped reads % |" ${INPDIR}/AlignSalmonGenomeStar/STAR/${f}Log.final.out | cut -d / -f 9
 done > ../Results/summary_star.txt






