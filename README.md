# ROAR-DB 

This repo contains codes to find similar proteins and their relative nucleotide sequences, the extracted data was used specifically building the ROAR-DB.  

## Finding ROs

Experimentally validated Rieske oxygenase (RO) enzymes were extracted from the literature, and their alpha subunit protein sequences were retrieved. A total of 72 RO sequences were used as queries to search for close homologs using Blastp (e-value of 0). 19207 sequences were then used to construct a Hidden Markov Model (HMM) with HMMER. The PF00355 domain code was used to get Rieske domain-containing proteins from UniProt. The HMM model of ROs (HMMRO) was used to classify these proteins with a threshold (e-value of 1e-10), and the results were stored in the ROAR MongoDB database.

![example_output](/img/ROARMongoDB.png)

## Extracting Operons

**Equation 1:** To calculate the relative importance of the gene group:

The flanking genes around the alpha subunits were determined for the ROAR DB. The surrounding genes were filtered according to 10kb left and 10kb right from the alpha subunit. Annotations were grouped and the gene groups for each RO were scored according to their relative distance and occurrence using Equation 1. These scores were normalized while calculating the heatmaps for each group. 

Score = (1 / ((Σ relative distances + ε) / #genes)) × ((#genes)^2 / sample size) × 100


![example_output](/img/Operon.png)
## Analysing Co-occurences 

Clustered heatmaps for each group were prepared using the libraries of matplotlib and seaborn to analyze co-occurred gene groups.

![example_output](/img/WordCloudR.png)