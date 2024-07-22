### Annotation
#### Method
Genome-guided transcript assemblies were made from close to 1 billion bp of 2x151bp paired-end Illumina RNAseq reads using PERTRAN (Shu, unpublished; see Cooper *et al*, 2019). PASA (Haas *et al*, 2003) alignment assemblies were constructed using the PERTRAN output from the Rio RNAseq data along with sequences from known *S. bicolor* expressed sequence tags (ESTs) associated with the current reference genome.

As further described in Phytozome, loci were determined by transcript assembly alignments and/or EXONERATE alignments of proteins from *Arabidopsis thaliana*, soybean, maize, rice, foxtail, *Sorghum bicolor* BTx623, brachy, grape, and Swiss-Prot proteomes to the repeat-soft-masked *Sorghum bicolor* Rio genome using RepeatMasker (RepeatMasker Open-3.0 by AFA Smit, R Hubley & P Green, 1996-2011) with up to 2K BP extension on both ends unless extending into another locus on the same strand. Gene models were predicted by homology-based predictors, FGENESH+ (Salamov and Solovyev, 2000), FGENESH_EST (similar to FGENESH+, EST as splice site and intron input instead of protein/translated ORF), and GenomeScan (Yeh *et al*, 2001), PASA assembly ORFs (in-house homology constrained ORF finder) and from AUGUSTUS via BRAKER1 (Hoff *et al*, 2016). The best scored predictions for each locus were selected using multiple positive factors including EST and protein support, and one negative factor: overlap with repeats. The selected gene predictions were improved by PASA (Haas *et al*, 2003). PASA-improved gene model proteins were subject to protein homology analysis to above mentioned proteomes to obtain Cscore and protein coverage; PASA-improved transcripts were selected based on Cscore, protein coverage, EST coverage, and its CDS overlapping with repeats. Selected gene models were subject to Pfam analysis and gene models whose protein was more than 30% in Pfam TE domains were removed. For additional details, see [*Sorghum bicolor Rio v2.1 (Sorghum Rio) in Phytozome v12.1*](https://phytozome.jgi.doe.gov/pz/portal.html#!info?alias=Org_SbicolorRio_er).

<table>
  <tbody>
  <tr>
    <td>Sorghum line</td>
    <td> </td>
    <td>Rio </td>
  </tr>
  <tr>
    <td>Assembly information</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td>Assembly name	</td>
    <td></td>
    <td>SbicolorRio_v2</td>
  </tr>
    <tr>
      <td>Assembly date</td>
      <td></td>
      <td>n/a</td>
    </tr>
    <tr>
      <td>Assembly accession</td>
      <td></td>
      <td>n/a</td>
    </tr>
    <tr>
      <td>WGS accession</td>
      <td></td>
      <td>JADDXV000000000</td>
    </tr>
    <tr>
      <td>Assembly provider	</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Sequencing description</td>
      <td>Sequencing technologies:	</td>
      <td>Illumina HiSeq 2500</td>
    </tr>
    <tr>
      <td></td>
      <td>Sequencing method</td>
      <td></td>
    </tr>
    <tr>
      <td></td>
      <td>Sequencing method	</td>
      <td></td>
    </tr>
    <tr>
      <td></td>
      <td>Genome coverage:	</td>
      <td>74.928x</td>
    </tr>
    <tr>
      <td>Assembly description	</td>
      <td>Assembly methods:	</td>
      <td>FALCON v. 2.2</td>
    </tr>
    <tr>
      <td></td>
      <td>Construction of pseudomolecules</td>
      <td></td>
    </tr>
    <tr>
      <td>Finishing strategy	</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
    <td>NCBI submission	</td>
    <td></td>
    <td>Submitted (26-OCT-2020)</td>
    </tr>
    <tr>
      <td>Publication:	</td>
      <td></td>
      <td>Cooper <em>et al</em> (2019)</td>
    </tr>
    <tr>
      <td>Assembly statistics	</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Number of contigs	</td>
      <td></td>
      <td>3,830</td>
    </tr>
    <tr>
      <td>Total assembly length (Mb)	</td>
      <td></td>
      <td>729</td>
    </tr>
    <tr>
      <td>Contig N50 (Mb)	</td>
      <td></td>
      <td>0</td>
    </tr>
    <tr>
      <td>Annotations stats</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Total number of genes	</td>
      <td></td>
      <td>35,490</td>
    </tr>
    <tr>
      <td>Total number of transcripts	</td>
      <td></td>
      <td>41,048</td>
    </tr>
    <tr>
      <td>Average gene length</td>
      <td></td>
      <td>3,322</td>
    </tr>
    <tr>
      <td>Exons per transcript</td>
      <td></td>
      <td>5</td>
    </tr>
    </tbody>
</table>

*Source: NCBI, April 2021.*
