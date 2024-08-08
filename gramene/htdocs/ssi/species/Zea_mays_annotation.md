### Annotation

The gene annotations were produced with the CSHL gene pipeline developed under the [NAM project](https://nam-genomes.org). 
In summary, it is an automated, evidence-based method combining third-party software including Mikado, BRAKER and PASA. 
Gene models were filtered by conservation and Maker Annotation Edit Distance (AED) score, 
and then classified into protein_coding and misc_non-coding sets.

Gene models from version B73_RefGen_v4 have been mapped to the current Zm-B73-REFERENCE-NAM-5.0 and can be retrieved with the
[ID History Converter](http://plants.ensembl.org/Oryza_sativa/Tools/IDMapper)

Evidence-based predictions were directly inferred from assembled transcripts, which were generated using five different genome-guided transcript assembly programs and processed using Mikado to pick the optimal set of transcripts for each locus. To generate assembled transcripts, quality inspected RNA-seq reads were mapped to the genome. In order to pick the final transcripts, Mikado uses assembled transcripts combined with high-confidence splice junctions with the mapped reads as input, predicted ORFs for the assembled transcripts generated and homology results of transcripts to SwissProt (Viridiplantae) sequences.

Ab initio predictions were performed using BRAKER with both evidence-based predicted proteins and mapped RNA-seq reads as input. 

A working set (WS) of models was generated to capture the complete gene space by combining evidence based and non-overlapping BRAKER gene models using BEDtools. Additional structural improvements on the WS models were completed using the software PASA. Transposable element related genes were filtered from the evidence and non-overlapping BRAKER sets using the TEsorter tool, which uses the REXdb database of TEs.

The TE filtered WS models were given AED scores using MAKER-P (v.3.0). Only models with AED < 0.75 passed to the high-confidence set (HCS). The HCS gene models were further classified based on homology to related species, and assigned coding and non-coding biotypes. The HCS gene models were checked for missing start and stop codons. The CDS boundaries of the transcripts were modified based on conserved start codon positions or extended to a start or stop codon whenever possible. All conserved genes in addition to lineage-specific genes that had a complete CDS were marked as protein-coding. The remaining lineage-specific genes were marked as non-coding. HCS gene models were checked and potentially split or merged using the GFF3toolkit. Gene ID assignment was made as per [MaizeGDB nomenclature schema](https://www.maizegdb.org/nomenclature). 

#### Repeats

Repeat features were annotated with the RepeatMasker pipeline with 
RepBase and wessler-bennetzen-2015/TE_12-Feb-2015_15-35 libraries.
