### Annotation
Genome annotation was carried out as described by (Deschamps *et al*, 2018). First, genome repeats were masked using RepeatMasker and a curated sorghum specific repeats file from Repbase. The repeat-masked genome was used as input to two categories of gene predictors. De novo gene prediction programs Fgenesh (Solovyev *et al*, 1994), Augustus (Stanke and Waack, 2003), and SNAP (Bromberg and Rost, 2007) were run under default parameters and the training sets used were monocots, maize, and rice, respectively. The EST, cDNA, long-read evidence-based gene structure modelers GMAP and PASA, as well as the protein evidence-based gene structure modeler SPLAN were also run. Long read sequences of BTx623 line of sorghum from NCBI, along with other sorghum EST’s and cDNA were used as the evidence set to PASA. Other non-sorghum Poales EST, cDNA sequences from NCBI, and monocot transcripts from phytozome were used as additional closely related species evidence for gene prediction with GMAP. Uniref100 plant protein sequences were used as an evidence dataset for gene structure prediction using SPLAN. All gene annotation files were run through EvidenceModeler and the output used to polish the gene boundaries in PASA. The final PASA annotation file was combined with tRNA predictions file from tRNA-ScanSE to obtain the final structural annotation file, along with fasta sequences of protein, CDS, cDNA and gene. For additional details, see (Deschamps *et al*, 2018).

<table>
  <tbody>
    <tr>
      <td>Sorghum line	</td>
      <td></td>
      <td>**Tx430**</td>
    </tr>
    <tr>
      <td>Assembly information	</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Assembly name	</td>
      <td></td>
      <td>Corteva_Sorghum_ONT_TX430_1.0
</td>
    </tr>
    <tr>
      <td>Assembly date	</td>
      <td></td>
      <td>n/a
</td>
    </tr>
    <tr>
      <td>Assembly accession	</td>
      <td></td>
      <td>n/a
</td>
    </tr>
    <tr>
      <td>WGS accession	</td>
      <td></td>
      <td>QWKM00000000
</td>
    </tr>
    <tr>
      <td>Assembly provider	</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Sequencing description	</td>
      <td>Sequencing technologies:	</td>
      <td>Oxford Nanopore MiniION
</td>
    </tr>
    <tr>
      <td></td>
      <td>Sequencing method	</td>
      <td></td>
    </tr>
    <tr>
      <td></td>
      <td>Genome coverage:	</td>
      <td>90.0x
</td>
    </tr>
    <tr>
      <td>Assembly description	</td>
      <td>Assembly methods:	</td>
      <td>SMARTdenovo v. June-2017; CANU v. 1.6
</td>
    </tr>
    <tr>
      <td></td>
      <td>Construction of pseudomolecules	</td>
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
      <td>Submitted (28-AUG-2018)

</td>
    </tr>
    <tr>
      <td>Publication:	</td>
      <td></td>
      <td>Deschamps et al (2018)
</td>
    </tr>
    <tr>
      <td>Assembly statistics	</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Number of contigs	</td>
      <td></td>
      <td>723</td>
    </tr>
    <tr>
      <td>Total assembly length (Mb)	</td>
      <td></td>
      <td>666
</td>
    </tr>
    <tr>
      <td>Contig N50 (Mb)	</td>
      <td></td>
      <td>3</td>
    </tr>
    <tr>
      <td>Annotations stats	</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td>Total number of genes	</td>
      <td></td>
      <td>36,937
</td>
    </tr>
    <tr>
      <td>Total number of transcripts	</td>
      <td></td>
      <td>49,928
</td>
    </tr>
    <tr>
      <td>Average gene length	</td>
      <td></td>
      <td>3,252
</td>
    </tr>
    <tr>
      <td>Exons per transcript	</td>
      <td></td>
      <td>4</td>
    </tr>
  </tbody>
</table>

*Source: NCBI, April 2021.*

