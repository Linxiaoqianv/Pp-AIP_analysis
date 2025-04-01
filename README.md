# New Autoinducing Peptides Regulate Antibiotic Production for Sculpting Microbiome

Scripts for data analysis, statistics, and plotting figures in the manuscript.

## Agr_mining
This folder is for genome mining of the agr gene cluster.

- Requirement:
  - antiSMASH 6.1.1 (For installation, please refer to: https://docs.antismash.secondarymetabolites.org/install/ ）

- ./Agr_mining/Agr_mining.sh
  
    The Linux script for genome mining of BGC.

    input: Total genome <*.fna.gz>.

- ./Agr_mining/00.data/antismash_pfam_PF04647_protein.tsv
  
    Final result file.

## Metagenome_analysis
This folder is for root metagenome analysis, including bacterial diversity analysis and figure drawing.

- Requirement:
  - fastp 
  - bowtie2
  - kraken2
  - bracken
  - R

- ./Metagenome_analysis/Metagenome_analysis.sh
  
    Linux script for obtaining taxonomic matrix.

    input: Raw reads of sequencing <*.fq.gz>.

- ./Metagenome_analysis/Pp-AIP_metagenome_Figure.R
  
    R script for bacterial diversity analysis and data visualization

## Phylogenetic_tree_of_isolates
This folder is for isolated strain analysis, related to Figure 4C.

- Requirement:
  - mafft 
  - trimal
  - fasttree
  - itol (online tools: https://itol.embl.de/itol.cgi)

- ./Phylogenetic_tree_of_isolates/Construction_tree.sh
  
    Linux script for construction of 16S rRNA gene phylogenetic tree.

    input: Total 16S rRNA gene sequences <latest2.total.16s.fa>.

    output: <latest2.16s_all_aligned_filtered.tree>.

    Output file is imported into itol for phylogenetic tree visualization. <dataset_color_strip_template_family.txt>, <dataset_external_shapes_template_diameter.txt>, and <dataset_gradient_template_rate.txt> are annotation files.



Reference: New Autoinducing Peptides Regulate Antibiotic Production for Sculpting Microbiome
