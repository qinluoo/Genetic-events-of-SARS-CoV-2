# Genetic-events-of-SARS-CoV-2
A method to identify the genetic events in SARS-CoV-2 genome.
# Introduction
In our previous study “Co-mutation modules capture the evolution and transmission patterns of SARS-CoV-2.” (https://doi.org/10.1093/bib/bbab222), the SARS-CoV-2 population was classified into different groups based on the co-mutation modules instead of the phylogenetic tree. Each group corresponds to a set of specific co-mutations that captured the vital evolutionary information of SARS-CoV-2 and the evolutionary relationship between variants accurately. Further, we used specific co-mutations in different variants and the association of mutation and recombination to identify genetic events in SARS-CoV-2 genomes.
# Software Requirment
We require a Windows or linux system with the Python software.Python (3.7.3 64-bit) were used in our experiments.
# Usage of "cmm-recom.py"
"cmm-recom.py" was used to label the co-mutations in the SARS-CoV-2 genomes and to determine the genetic events (i.e. mutation or recombination).

Open the cmd in windows and import:

python cmm-recom.py sequence.fasta reference.fasta co-mutations group_mut recombination.xls mutation.xls group.xls G0.xls

Input Arguments:

sequence.fasta: the SARS-CoV-2 genomes after alignment with reference sequence (EPI_ISL_402125).

reference.fasta: EPI_ISL_402125 is used as the reference sequence.

co-mutations: a list including the identified co-mutations.

group_mut: a list including groups and corresponding co-mutation modules.

Output files:

recombination.xls: sequences identified as recombinants.

mutation.xls: sequences identified as mutants of a group.

group.xls: sequences belonging to groups.

G0.xls: sequences belonging to G0.

# Usage of "world_map.py"
"world map.py" was used to draw the spatial distribution map.

Package required: pyecharts 0.5.10

Open the cmd in windows and import:

python world_map.py country_count.xls

Input Arguments:

country_count.xls: a list including countries and corresponding sequence counts

output: the world heatmap
