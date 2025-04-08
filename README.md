# ROAR-DB 

This repo contains codes to find similar proteins and their relative nucleotide sequences, the extracted data was used specifically to build the ROAR-DB.  

## Requirements

**Libraries**
- [Linux]()
- [Python]()
- [PyMongo](), [matplotlib](), [biopython](), [pandas](), [seaborn](), [numpy]() 

**Externals**
- [blastn](https://www.ncbi.nlm.nih.gov/books/NBK569861/)
- [makeblastdb](https://www.ncbi.nlm.nih.gov/books/NBK569861/) 
- [hmmer](http://hmmer.org/) 

**Data**
- [InterPro](https://www.ebi.ac.uk/interpro/)
- [NCBI](https://www.ncbi.nlm.nih.gov/)

## Finding ROs

Experimentally validated Rieske oxygenase (RO) enzymes were extracted from the literature, and their alpha subunit protein sequences were retrieved. A total of 72 RO sequences were used as queries to search for close homologs using Blastp (e-value of 0). 19207 sequences were then used to construct a Hidden Markov Model (HMM) with HMMER. The PF00355 domain code was used to get Rieske domain-containing proteins from UniProt. The HMM model of ROs (HMMRO) was used to classify these proteins with a threshold (e-value of 1e-10), and the results were stored in the ROAR MongoDB database.

![example_output](/img/ROARMongoDB_low.png)

## Extracting Operons


The flanking genes around the alpha subunits were determined for the ROAR DB. The surrounding genes were filtered according to 10kb left and 10kb right from the alpha subunit. Annotations were grouped and the gene groups for each RO were scored according to their relative distance and occurrence using Equation 1. These scores were normalized while calculating the heatmaps for each group. 
 ```
Equation 1: To calculate the relative importance of the gene group:
Score = (1 / ((Σ relative distances + ε) / #genes)) × ((#genes)^2 / sample size) × 100
 ```

![example_output](/img/Operon.png)
## Analysing Co-occurences 

By analyzing co-occurring gene groups, consistently associated gene clusters have been identified for each RO group and each RO enzyme. For this purpose, heatmap and word cloud analyses are performed, and the results can be visualized using Sankey diagrams. In this way, functional predictions can be made for previously unknown or unannotated genes and enzymes, providing data that can support the formulation of hypotheses.

![example_output](/img/WordCloudR.png)