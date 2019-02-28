##################################################################
#                          README                                #
##################################################################


The following directory contains all shell scripts necessary for the analysis of an RNAseq experiment, mainly alignment using star and de novo assembly of a viral genome using Trinity and SPAdes. 

Feature counts was used to count reads mapped to the genome, and the raw counts were later using in the Limma-voom pipeline for differential expression analysis. This is contained in the Differential expression directory. 

Directory structure: 

/Scripts 

The scripts directory contains all the scripts to perform the alignment and read counts of raw reads from RNAseq experiment. The pipeline used in this analysis is presented as an image in the /Docs directory





                                                                              ,///////.           ./
          __                    ___////,,,.._        ,//               __,---//////////,        .///
         /o \/        __       /O            ``--.._///             ,-'  ) ) ) ) ) )''////_    /////
<><      \__/\ __    /o \/     \__  \               ,'          _,-' ))`. ) ) ) ) ) ) ) ) )`-.//////
 <><          /o \/  \__/\      __)  \           ___`.         / ()_)))))\ ) ) ) ) ) ) ) ) )////////
   <><     __ \__/\_            \____/__,,,---''''  \\\        \____ )))))) ) ) ) ) ) ) ) ) \\\\\\\\
 <><      /o \/  /o \/                               `\\        `````.)))/ ) ) ) ) ) ) ) ),-`\\\\\\\
<><       \__/\  \__/\                                           ___,')),') ) ) ) )_),,--'    \\\\\\
<><       __     __                  ___////,,,.._        ,//   (_______.\\)_),--'"            `\\\\
>< <><   /o \/  /o \/     __        /O            ``--.._///             -\\\                    `\\\
><  <><  \__/\  \__/\    /o \/      \__  \               ,'                \\\                     `\
  <><                    \__/\       __)  \           ___`.
<><                __                \____/__,,,---''''  \\\
                  /o \/                                   `\\
           __     \__/\    ___////,,,.._        ,//                                      ,///////.           ./
          /o \/           /O            ``--.._///                                __,---//////////,        .///
          \__/\           \__  \               ,'                              ,-'  ) ) ) ) ) )''////_    /////
                           __)  \           ___`.                          _,-' ))`. ) ) ) ) ) ) ) ) )`-.//////
                           \____/__,,,---''''  \\\                        / ()_)))))\ ) ) ) ) ) ) ) ) )////////
                                                `\\                       \____ )))))) ) ) ) ) ) ) ) ) \\\\\\\\
                                                      Krogg                `````.)))/ ) ) ) ) ) ) ) ),-`\\\\\\\
                                                                            ___,')),') ) ) ) )_),,--'    \\\\\\
                                                                           (_______.\\)_),--'"            `\\\\
                                                                                    -\\\                    `\\\
                                                                                      \\\                     `\




