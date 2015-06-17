<div class="plain-box"><h2>What's in Version 2</h2>
This release is the second stable release of oge.gramene.org. It includes the most complete and up-to-date assemblies and annotations of the genomes studied by Oryza Genome Evolution (OGE) project, with 12 complete whole genomes from International Oryza Map Alignment Project (I-OMAP) and 4 chromosome 3 short arms from <a href='http://www.omap.org/'>OMAP project</a>.

<ul> 
<li> 
Compared with version 1 release, this release has updated Oryza meridionalis assembly, which is the final version of the assembly that was accepted by NCBI genbank. This version adjusted many gap sizes in order to conform to the NCBI submission guidelines. Accordingly MAKER-P annotation was reran and subsequence annotations were repeated with the new version of gene models. However the protein compara analysis which generates gene trees was not redone, instead the O.meridionalis version2 MAKER-P gene names were mapped to version1 names based on location projection and protein sequence similarity search.
</li>
<li>
RNA-seq expression profiles from three tissues (root, leaf, panicle) were added for 10 whole genomes. They are O.sativa. Japonica, O.rufipogon, O.nivara, O.meridionalis, O.glaberrima, O.barthii, O.glumaepatula, O.punctata, O.brachyantha, L.perrieri.      
</li>
<li>
Nuclear organelle DNA insertion (norgDNA) annotation were provided by Christos in Doreen Ware's lab and displayed as a dna feature track.
</li>
<li>
lastz-net WGA of O.meridionalis against the other 10 whole genome assemblies were repeated with the updated O.meridionalis assemlby.
Below is the summary table of the WGAs
<ul>
<li>Whole Genomes</li>
<style type="text/css">
	table.tableizer-table {
	border: 1px solid #CCC; font-family: Arial, Helvetica, sans-serif
	font-size: 12px;
} 
.tableizer-table td {
	padding: 2px;
	margin: 1px;
	border: 1px solid #CCC;
}
.tableizer-table th {
	background-color: #0B610B; 
	color: #FFF;
	font-weight: bold;
}
</style><table class="tableizer-table">
<tr class="tableizer-firstrow"><th>Species</th><th>OGE version</th><th>Genbank WGS Proj Accession</th><th>Genbank Assembly ID</th><th>Genbank Bioproject</th></tr>
 <tr><td><a href='/Oryza_sativa/Info/Index'><i>O. sativa ssp. japonica</i></a></td><td>v7</td><td>NA</td><td>NA</td><td>NA</td></tr>
 <tr><td><a href='/Oryza_indica/Info/Index'><i>O. sativa ssp. indica</i></a></td><td>v2</td><td>AAAA00000000.2</td><td>GCA_000004655.2 </td><td>PRJNA361</td></tr>
 <tr><td><a href='/Oryza_nivara/Info/Index'><i>O. nivara</i></a></td><td>v1.0</td><td>AWHD00000000.1</td><td>GCA_000576065.1</td><td>PRJNA48107</td></tr>
 <tr><td><a href='/Oryza_rufipogon/Info/Index'><i>O. rufipogon</i></a></td><td>v1.1</td><td>NA</td><td>NA</td><td>PRJEB4137</td></tr>
 <tr><td><a href='/Oryza_barthii/Info/Index'><i>O. barthii</i></a></td><td>v3</td><td>ABRL00000000.2</td><td>GCA_000182155.2</td><td>PRJNA30379</td></tr>
 <tr><td><a href='/Oryza_glaberrima/Info/Index'><i>O. glaberrima</i></a></td><td>v2</td><td>ADWL00000000.1</td><td>NA</td><td>PRJNA13765</td></tr>
 <tr><td><a href='/Oryza_glumaepatula/Info/Index'><i>O. glumaepatula</i></a></td><td>v1.5</td><td>ALNU00000000.2</td><td>GCA_000576495.1</td><td>PRJNA48429</td></tr>
 <tr><td><a href='/Oryza_meridionalis/Info/Index'><i>O. meridionalis</i></a></td><td>v1.3</td><td>ALNW00000000.2</td><td>GCA_000338895.2</td><td>PRJNA48433</td></tr>
 <tr><td><a href='/Oryza_longistaminata/Info/Index'><i>O. longistaminata</i></a></td><td>v117</td><td>NA</td><td>NA</td><td>PRJNA245492</td></tr>
 <tr><td><a href='/Oryza_punctata/Info/Index'><i>O. punctata</i></a></td><td>v1.2</td><td>AVCL00000000.1</td><td>GCA_000573905.1</td><td>PRJNA13770</td></tr>
 <tr><td><a href='/Oryza_brachyantha/Info/Index'><i>O. brachyantha</i></a></td><td>v1.4b</td><td>AGAT00000000.1</td><td>GCA_000231095.2</td><td>PRJNA70533</td></tr>
 <tr><td><a href='/Leersia_perrieri/Info/Index'><i>L. perrieri</i></a></td><td>v1.4</td><td>ALNV00000000.2</td><td>GCA_000325765.3</td><td>PRJNA163065</td></tr>
</table>

<li>Chromosome 3 short arms</li>
<table class="tableizer-table">
<tr class="tableizer-firstrow"><th>Species</th><th>OGE version</th><th>Genbank WGS Proj Accession</th><th>Genbank Assembly ID</th><th>Genbank Bioproject</t
h></tr>
 <tr><td><a href='/Oryza_granulata3s/Info/Index'><i>O. granulata</i></a></td><td>v03</td><td>ALNT00000000.1</td><td>GCA_000325655.2</td><td>PRJNA13788</td></tr>
 <tr><td><a href='/Oryza_officinalis/Info/Index'><i>O. officinalis</i></a></td><td>v03</td><td>JJMQ00000000.1</td><td>GCA_000717455.1</td><td>PRJNA239525</td></
tr>
 <tr><td><a href='/Oryza_minutabb3s/Info/Index'><i>O. minuta</i> (BB)</a></td><td>v03</td><td>JJNN00000000.1</td><td>GCA_000632695.1</td><td>PRJNA239524</td></t
r>
<tr><td><a href='/Oryza_minutacc3s/Info/Index'><i>O. minuta</i> (CC)</a></td><td>v03</td><td>JJNN00000000.1</td><td>GCA_000632695.1</td><td>PRJNA239524</td></tr>
</table>
</ul>
</li>
</ul>
</div>
