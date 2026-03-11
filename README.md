# bakta_kofam_scan

A program that uses a query_table.tsv for searching up genes in bakta annotated genomes.

## Requirements
`main.py` assumes presence of the query_table.tsv and .tsv summary output from bakta. While bakta result is universal a 
query_table.tsv is in-house made.  
Query table contains genes, their function and probiotic trait.

## Usage
`python main.py --genes query_table.tsv --genomes genomes/`

Directory structure - in genomes directory one should place a directory for a specific genome. In that directory a bakta
summary output should be placed. In that directory an output will be produced.

## Additional info
This repository was created specificaly for a work-related task that goes: "I have a table with probiotic genes found in
bacteria X. Using bakta annotation for our bacterial strains, please find me if and what probiotic genes does these
strains contain"