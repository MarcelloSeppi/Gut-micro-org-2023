############################################################
0.	General remarks.
############################################################

This is a December 2014 release of archaeal COGs constructed by Eugene Koonin's group at the National Center for Biotechnology Information (NCBI), National Library of Medicine (NLM), National Institutes of Health (NIH).

#-----------------------------------------------------------
0.1.	Citation.

Makarova KS, Wolf YI, Koonin EV.

Archaeal clusters of orthologous genes (arCOGs): An update and application for analysis of shared features between Thermococcales, Methanococcales and Methanobacteriales.

Life 2015 [submitted]

#-----------------------------------------------------------
0.2.	Contacts.

Biology of Archaea, data usage, collaborations:

Eugene V. Koonin <koonin@ncbi.nlm.nih.gov>

Biology of Archaea, details of arCOG annotation:

Kira S. Makarova <makarova@ncbi.nlm.nih.gov>

Software, procedures, technical questions etc.:

Yuri I. Wolf <wolf@ncbi.nlm.nih.gov>

############################################################
1.	Notes
############################################################

#-----------------------------------------------------------
1.1.	December 2014 release.

This release corresponds is described in a paper, submitted to Life journal in January 2015.

#-----------------------------------------------------------
1.2.	File formats.

Please note that we introduce minor modifications to the file formats compared to the previous versions. See section #2 for details.

#-----------------------------------------------------------
1.3.	Unannotated proteins.

In a number of published genomes several important proteins were missed in the annotation. Where identified, missing proteins were given internal numerical IDs and added to the arCOG list; sequences were added to the arCOG FASTA file.

############################################################
2.	Files
############################################################

   34178 ar14.genome.tab
   10610 ar14.genclade.tab
75013928 ar14.fa.gz
35690913 ar14.arCOG.csv
  861463 ar14.arCOGdef.tab
    1117 funclass.tab
  256389 ar14.sig.tab
 1052399 ar14.tm.tab

#-----------------------------------------------------------
2.1.	ar14.genome.tab

Contains list of genomes (168). Tab-delimited, format:

<genome-name> <ncbi-taxid> <status> <ncbi-taxonomy>

* Example:

Methanopyrus_kandleri_AV19_uid57883	190192	OLD	cellular organisms;Archaea;Euryarchaeota;Methanopyri;Methanopyrales;Methanopyraceae;Methanopyrus;Methanopyrus kandleri;Methanopyrus kandleri AV19

* Comments:

The field <genome-name> serves as a subdirectory name on NCBI Genomes FTP site <ftp://ftp.ncbi.nih.gov/genomes/Bacteria> (as of February 2014). The field <status> shows the status of this genome relative to the 2012 arCOG release. "OLD" means present in 2012; "NAME_CHANGE" means present in 2012 under a different name; "NEW" means new for 2014 and available from the NCBI FTP site at February 2014; "STILL_MISSING" means new for 2014 and not available from the NCBI FTP site at February 2014 (downloaded from GenBank). The field <ncbi-taxonomy> is a semicolon-delimited lineage according to the NCBI Taxonomy database (as of February 2014).

#-----------------------------------------------------------
2.2.	ar14.genclade.tab

Contains the breakdown of genomes into operational taxonomic clades and genome weights. Tab-delimited, format:

<genome-name> <clade> <weight>

* Example:

Pyrococcus_horikoshii_OT3_uid57753	E.Thermococci	3.48E-01

* Comments:

Clades and genome weights are described in the paper (see section #0.1).

#-----------------------------------------------------------
2.3.	ar14.fa.gz

* Contains:

Archaeal protein sequences (394110) in FASTA format (gzipped)

* Example:

>gi|118430838|ref|NP_146899.2| putative mercury ion binding protein
[Aeropyrum pernix K1]
MIIFKRHSQAILFSHNKQEKALLGIEGMHCEGCAIAIETALKNVKGIIDTKVNYSRGSAI
VTFDDTLVSINDILEHYIFKVPSNYRAKLVSFIS

* Defline:

The first word of the defline always starts with "gi|<protein-id>"; the rest of the defline is optional.

#-----------------------------------------------------------
2.4.	ar14.arCOG.csv

Contains list of archaeal orthology domains. Comma-delimited, format:

<domain-id>, <genome-name>, <protein-id>, <protein-length>, <domain-start>, <domain-end>, <arCOG-id>, <membership-class>, [deprecated fields]

* Example:

529078083,Halorhabdus_tiamatea_SARL4B_uid214082,529078083,232,1,232,arCOG03458,0,,

* Comments:

In this version the fields <domain-id> and <protein-id> are identical and both normally refer to GenBank GIs (see section #1.1 for exceptions). Thus neither <domain-id> nor <protein-id> are necessarily unique in this file (this happens when a protein consists of more than one orthology domains, e.g. 48478501).

The <membership-class> field indicates the nature of the match between the sequence and the arCOG consensus:

0 - the domain matches the arCOG consensus;

1 - the domain is significantly longer than the arCOG consensus;

2 - the domain is significantly shorter than the arCOG consensus;

3 - partial match between the domain and the arCOG consensus.

Deprecated fields should be ignored in this version.

#-----------------------------------------------------------
2.5.	ar14.arCOGdef.tab

Contains list of arCOG annotations. Tab-delimited, format:

<arCOG-id> <functional-class> <gene-name> <arCOG-annotation> <supercluster> <profile-pfam> <profile-cdd> <profile-tigrfam>

* Example:

arCOG00024	C	GlpK	Glycerol kinase	COG00554	pfam00370,pfam02782	cd07786	TIGR01311

* Comments:

Functional classes are described in the file funclass.csv (see section #2.6). Undefined gene names are denoted with "-" character. The <supercluster> field provides hierarchy above the arCOG level; where possible superclusters correspond to 2003 NCBI COGs [doi:10.1186/1471-2105-4-41]. The <profile-pfam>, <profile-cdd> and <profile-tigrfam> fileds are comma-delimited lists of CDD [doi:10.1093/nar/gku1221] profiles with significant (e-value<1e-5) hits to arCOG members.

#-----------------------------------------------------------
2.6.	funclass.tab

Contains list of functional classes and major categories. Tab-delimited, format:

<class-id> <category> <description>

* Example:

3	3	METABOLISM
C	3	Energy production and conversion

* Comments:

N/A

#-----------------------------------------------------------
2.7.	ar14.sig.tab

Contains list of proteins with predicted signal peptides. Tab-delimited, format:

<protein-id> Y

* Example:

330508413       Y

* Comments:

Signal peptides were predicted using the SignalP v. 4.1c program [doi:10.1038/nmeth.1701]; the union of the three predictions (gram-negative, gram-positive and eukaryotic models) was used.

#-----------------------------------------------------------
2.8.	ar14.tm.tab

Contains list of proteins with predicted transmembrane segments. Tab-delimited, format:

<protein-id> <num-segments>

* Example:

337283700	3

* Comments:

Transmembrane segments were predicted using the TMMHMM v. 2.0c program with default parameters [doi:10.1006/jmbi.2000.4315].
