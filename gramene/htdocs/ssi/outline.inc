<div class="plain-box" style="float:right; padding:20px"><h2>What's New in Release 7</h2>

<h3>New Functionality</h3>

<ul>
	<li>
	<p>Oryza CLIMtools portal with interactive web-based views of environment x genome associations, RiboSNitch prediction, and correlations between the local environment and a pool of curated phenotypes.</p>
	</li>
	<li>
	<p>Updated RAP-DB gene annotations (Sept. 2022), see for example <a href="http://oryza.gramene.org?idList=Os01g0104200">Os01g0104200</a>. The <a href="https://oryza.gramene.org/?filters=%7B%22status%22:%22init%22,%22operation%22:%22AND%22,%22negate%22:false,%22marked%22:false,%22leftIdx%22:0,%22rightIdx%22:3,%22children%22:%5B%7B%22fq_field%22:%22capabilities%22,%22fq_value%22:%22pubs%22,%22name%22:%22publication%22,%22category%22:%22Curated%22,%22leftIdx%22:1,%22rightIdx%22:2,%22negate%22:false,%22showMenu%22:false,%22marked%22:true%7D%5D,%22showMarked%22:true,%22showMenu%22:false,%22moveCopyMode%22:%22%22,%22searchOffset%22:0,%22rows%22:20%7D&amp;genomes=" rel="nofollow">functionally curated gene annotations from GeneRIF and RAP-DB</a> include a Papers tab listing directly relevant PubMed articles. The associated Plant Ontology (PO) and Trait Ontology (TO) terms are also searchable, for example: <a href="https://oryza.gramene.org/?filters=%7B%22status%22:%22init%22,%22operation%22:%22AND%22,%22negate%22:false,%22marked%22:false,%22leftIdx%22:0,%22rightIdx%22:3,%22children%22:%5B%7B%22fq_field%22:%22PO__ancestors%22,%22fq_value%22:%229005%22,%22name%22:%22root%22,%22category%22:%22Plant%20Ontology:%20anatomy%22,%22leftIdx%22:1,%22rightIdx%22:2,%22negate%22:false,%22showMenu%22:false,%22marked%22:true%7D%5D,%22showMarked%22:true,%22showMenu%22:false,%22moveCopyMode%22:%22%22,%22searchOffset%22:0,%22rows%22:20%7D&amp;genomes=" rel="nofollow">PO: Anatomy: root</a> or <a href="https://oryza.gramene.org/?filters=%7B%22status%22:%22init%22,%22operation%22:%22AND%22,%22negate%22:false,%22marked%22:false,%22leftIdx%22:0,%22rightIdx%22:3,%22children%22:%5B%7B%22fq_field%22:%22TO__ancestors%22,%22fq_value%22:%22276%22,%22name%22:%22drought%20tolerance%22,%22category%22:%22Trait%20ontology%22,%22leftIdx%22:1,%22rightIdx%22:2,%22negate%22:false,%22showMenu%22:false,%22marked%22:true%7D%5D,%22showMarked%22:true,%22showMenu%22:false,%22moveCopyMode%22:%22%22,%22searchOffset%22:0,%22rows%22:20%7D&amp;genomes=" rel="nofollow">TO: drought tolerance</a>.</p>
	</li>
	<li>
	<p>Added new links to the Gene Tree Curation Tool in the Homology tab of Nipponbare genes. Example: <a href="http://curate.gramene.org/admin/curationUI/oryza_v5/?since=4479&amp;gene=Os05g0569300" rel="nofollow">Os05g0569300</a>.</p>
</li>
	<li>
<p>Added new links to the GrameneOryza search interface results, see for example: <a href="https://oryza-ensembl.gramene.org/Oryza_sativa/Gene/GrameneSearch?g=Os05g0113900;r=5:738208-739079;t=Os05t0113900-01">Os05g0113900</a>.</p>
</li>
</ul>


<h3>New & Updated Data</h3>

<p><em><strong>Genes</strong></em></p>

<ul>
<li>Updated gene names for the Magic16 genomes (except Nipponbare) with new prefixes: OsZS97, OsMH63, OsIR64, OsPr106, OsGoSa, OsLaMu, OsLima, OsKYG, OsLiXu, OsN22, OsNaBo, OsARC, OsKeNa, OsCMeo, and OsAzu.
</li>
</ul>

<p><em><strong>Genetic Variation</strong></em></p>

<p>Approximately 19 million SNPs called per each of four <a href="https://github.com/warelab/release-notes/blob/main/oryza/www.nature.com/articles/s41597-020-0438-2">Magic16</a> genomes from the <a href="https://doi.org/10.1186/2047-217x-3-7" rel="nofollow">3000 rice genome project</a> were called by Dr. Rod Wing&#39;s group from resequencing reads using <a href="https://gatk.broadinstitute.org/hc/en-us" rel="nofollow">GATK4</a> software:</p>

<ul>
	<li><a href="https://oryza-ensembl.gramene.org/Oryza_sativa117425"><em>Oryza sativa</em> circum-Basmati var. ARC 10497</a></li>
	<li><a href="https://oryza-ensembl.gramene.org/Oryza_sativazs97"><em>Oryza sativa</em> Xian/Indica-1A var. Zhenshan 97</a></li>
	<li><a href="https://oryza-ensembl.gramene.org/Oryza_sativa125827"><em>Oryza sativa</em> Xian/Indica-3B2 var. Liu Xu</a></li>
	<li><a href="https://oryza-ensembl.gramene.org/Oryza_aus"><em>Oryza sativa</em> circum-Aus1 var. N22</a></li>
</ul>

<p>Genetic variation for <em>O. sativa</em> Japonica Nipponbare IRGSP1 was imported from the Euorpean Variation Archive (EVA). It consists of the following data sets:</p>

<ul>
	<li>25.8 million SNPs from Duitama <em>et al</em> (2015)</li>
	<li>3 million SNPs from BGI (2004)</li>
	<li>1.6 million SNPs from OMAP (2007)</li>
	<li>366K SNPs from the <a href="https://doi.org/10.1186/2047-217x-3-7" rel="nofollow">3K rice genome project</a></li>
	<li>157K SNPs from McNally <em>et al</em> (2009)</li>
	<li>1.3K SNPs from Zhao et al (2010)</li>
</ul>


