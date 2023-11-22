############################################################
#
############################################################

#-----------------------------------------------------------
# cog-20.cog.csv
#-----------------------------------------------------------

Comma-delimited plain text file assigning proteins to COGs

Columns:

1.	Gene ID (GenBank or ad hoc)

2.	NCBI Assembly ID

3.	Protein ID (GenBank if conforms to [A-Za-z0-9_]+\.[0-9]+ regex; ad hoc otherwise)

4.	Protein length

5.	COG footprint coordinates on the protein. "201-400" means "from position 201 to position 400"; "1-100=201-300" indicates a segmented footprint, 1-100 AND 201-300

6.	Length of the COG footprint on the proteins

7.	COG ID

8.	reserved

9.	COG membership class (0: footprint covers most of the protein and most of the COG profile; 1: footprint covers most of the COG profile and part of the protein; 2: footprint covers most of the protein and part of the COG profile; 3: partial match on both protein and COG profile)

10.	PSI-BLAST bit score for the match between the protein and COG profile

11.	PSI-BLAST e-value for the match between the protein and COG profile

12.	COG profile length

13.	Protein footprint coordinates on the COG profile

Each line corresponds to one instance of a COG in a protein-coding gene. Multidomain proteins are represented by multiple lines. Protein IDs can be shared between multiple genes. A combination of a gene ID and assembly ID uniquely identifies a location in the genome. The order of entries if arbitrary.

#-----------------------------------------------------------
# fun-20.tab
#-----------------------------------------------------------

Tab-delimited plain text file with descriptions of COG functional categories

Columns:

1.	Functional category ID (one letter)

2.	Hexadecimal RGB color associated with the functional category

3.	Functional category description

Each line corresponds to one functional category. The order of the categories is meaningful (reflects a hierarchy of functions; determines the order of display)

#-----------------------------------------------------------
# cog-20.def.tab
#-----------------------------------------------------------

Tab-delimited plain text file with COG descriptions

Columns:

1.	COG ID

2.	COG functional category (could include multiple letters in the order of importance)

3.	COG name

4.	Gene associated with the COG (optional)

5.	Functional pathway associated with the COG (optional)

6.	PubMed ID, associated with the COG (multiple entries are semicolon-separated; optional)

7.	PDB ID of the structure associated with the COG (multiple entries are semicolon-separated; optional)

Each line corresponds to one COG. The order of the COGs is arbitrary (displayed in the lexicographic order)

#-----------------------------------------------------------
# cog-20.org.csv
#-----------------------------------------------------------

Comma-delimited plain text describing genome assemblies

Columns:

1.	NCBI Assembly ID

2.	Organism (genome) name

3.	NCBI Tax ID of the assembly

4.	Taxonomic category used in COGs


Each line corresponds to one genome assembly. The order of the assemblies is meaningful (roughly corresponds to the taxonomic order; determines the order of display)

#-----------------------------------------------------------
# cog-20.tax.csv
#-----------------------------------------------------------

Comma-delimited plain text describing taxonomic categories

Columns:

1.	Taxonomic category in COGs

2.	Parent taxonomic category (self, if top of the hierarchy)

3.	NCBI Tax ID of the assembly

Each line corresponds to one genome taxonomic category. The order of the taxonomic category is meaningful (determines the order of display)

#-----------------------------------------------------------
# cog-20.fa.gz
#-----------------------------------------------------------

Gzipped FASTA file, containing sequences from all COGs
