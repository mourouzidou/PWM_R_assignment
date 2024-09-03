# Gene Expression Regulation in Gluconeogenesis and Glycolysis Pathways

This repository contains the analysis and code related to the study of gene expression regulation in the Gluconeogenesis (REACTOME ID: R-HSA-70263) and Glycolysis (REACTOME ID: R-HSA-70171) pathways. These pathways are crucial in the metabolism of sugars and understanding their regulation provides insights into the metabolic processes in humans.

## Overview

The goal of this project is to analyze the regulatory elements present in the upstream regions of genes involved in the Gluconeogenesis and Glycolysis pathways

  Gene Data Collection:
        Download gene names from the ENSEMBL database using the REACTOME IDs.
        Retrieve upstream 500bp sequences for these genes.
        Collect upstream sequences for a random set of 1000 genes to serve as a background dataset.

   Hypergeometric Analysis:
        Determine whether the overlap between the Gluconeogenesis and Glycolysis gene lists is greater than expected by chance using a hypergeometric test.

   Motif Analysis:
        Identify overrepresented 8-mer sequences in the upstream regions of genes in each pathway compared to the random background.
        Calculate the p-values for these motifs to assess their significance.

   Position Weight Matrix (PWM) Construction:
        Build PWMs for the most overrepresented motifs in each pathway.
        Scan the sequences from each pathway using the constructed PWMs and analyze the scores to determine regulatory similarities between pathways.

## Methods

* Retrieve gene names and upstream 500bp sequences for the Gluconeogenesis and Glycolysis pathways from the ENSEMBL database.
* Perform hypergeometric test to determine wether there are more shared genes between these two pathways than would be expected ny chance
* Identify 8-mer motifs that are overrepresented in each pathway compared to a random set of 1000 genes.
* Calculate the p-values for these motifs, focusing on those with p-value < 0.001 to identify significant regulatory sequences.
* Build Position Weight Matrices (PWMs) for the most significant motifs and use them to scan the sequences in both pathways
* Compare the regulatory patterns and identify similarities/differences in how these pathways might be regulated.
* Visualize the results using heatmaps and sequence logos to highlight the presence of such motifs across the gene sequences.

The substring with the smallest p-value for gluconeogenesis pathway is AATCTCGC and for glycolysis pathway AAAGACGC.

These logos display the frequency of each nucleotide at each position within the motif, indicating the most conserved bases that likely play a significant role in gene regulation.
<p align="center">
  <img src="https://github.com/user-attachments/assets/3639ecd4-af6e-4fb9-b4eb-559e8aa615e2" alt="Sequence Logos" width="400"/>
</p>
The heatmaps below represent the occurrence of these motifs across the genes in each pathway. Each row corresponds to a gene, and each column represents a position within the sequence. The red regions highlight where the motifs are strongly present. In both heatmaps, the widespread red color indicates a high occurrence of the key motifs across the genes in each pathway, signifying their potential importance in regulating the gene expressions in Glycolysis and Gluconeogenesis.
<p align="center">
  <img src="https://github.com/user-attachments/assets/20dad610-fca6-439b-a5ce-91ffd0acf630" alt="Heatmaps" width="400"/>
</p>




  
