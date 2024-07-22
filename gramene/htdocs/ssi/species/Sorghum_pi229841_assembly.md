### Assembly
The preliminary assembly was constructed using Canu and 28x coverage Pacbio long read data and then aligned to the reference (BTx623) using MUMmer [(Voelker et al, 2022)](https://www.sorghumbase.org/paper/ten-new-high-quality-genome-assemblies-for-diverse-bioenergy-sorghum-genotypes). Contigs were sorted using coordinate data from MUMmer and code developed by the Cooper Lab at UNCC (https://github.com/cponce-uncc/bioinformatics_tools). Contigs under 50k base pairs were discarded, as were contigs that were manually determined to not clearly map to a position after sorting. Final assembly statistics are as follows:

Final Contig Count: 1976

Final Base Pair Count: 598235854

Percent Coverage of Reference:81.950%
