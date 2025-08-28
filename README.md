# BAKTA_KOFAM_SCAN

This repo contains results and jupyter notebook scripts in order to conduct a scan for genes of interest in the WGS results folder.

## USAGE

If you wish to run the scan by yourself you should place the .ipnyb files in the WGS repository.

e.g.
```
WGS
|-bakta_quick_scan.ipnyb    # HERE
|-kofam_analysis.ipnyb      # HERE
|-kofamscan.sh              # HERE
|-3500
  |-3500_bafta_annotation
    |-3500.faa
    |-etc.
  |-etc.
|-etc.
```

Scripts were written in a way that they work only when placed correctly and only with bakta annotation done already. 

### Kofam scan

Firstly run `kofamscan.sh` in CLI.
1. Open terminal (Windows: open WSL, MacOS: open terminal)
2. Navigate to the WGS directory (copy this command: `cd <path/to/WGS/folder>` replace `<path/to/WGS/folder>` with your path)
3. Then write: `bash kofamscan.sh <path/to/fasta/file>` (eg. `bash kofamscan.sh 3500/3500_bakta_annotation/3500.faa`)

Then open jupyter notebook in order to run `kofam_analysis.ipnyb` (you might need to install jupyter first). Run each code block in jupyter by shift + Enter. Results should be written as .csv file in the bakta_annotation folder.

### Bakta scan

Quikc bakta scan is run by using a jupyter notebook the same way as in kofam scan (no need to run the .sh script beforehand)

## Results explanation

The csv (for kofam) and tsv (for bakta) contains results of the scan. They should ideally present all genes in the sequenced organisms that were annotaded with terms specified in jupyter notebook.