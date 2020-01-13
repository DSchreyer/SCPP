# SCPP - Single-cell RNA-seq processing pipeline
This repository contains scripts that were created during my Master's thesis including the script for the SCPP. SCPP was tested under a Linux environment and is therefore not applicable with macOS or Windows. Only the SCPP script ("SCPP.sh") and the quality control and filtering script "sc_analysis_qc.R" are necessary to start the pipeline. The other scripts can be disregarded. For beginners, the "start_SCPP.sh" script might be interesting, because it contains information about the required and the optional parameters. 
SCPP is also runnable with the "start_SCPP.sh" script. To execute this script, make it executable with `$ chmod +x start_SCPP.sh` and execute it with `$ ./start_SCPP.sh`.

## How to use
1. Clone this repository into a directory on your computer
```bash
$ git clone https://github.com/DSchreyer/Master-Thesis.git
```
2. make SCPP.sh executable
```bash
$ chmod +x SCPP.sh
```
3. start SCPP from the command line
```bash
$ ./SCPP.sh --option1 " " --option2  ...
```



### Options available in this pipeline
| Option  | Default | Function                                    | Branch |
|:---------------------:|:------------------:|:------------------------------------------------------:|:--------------:|
| --data           | -                | Path to directory with scRNA-seq FASTQ files         |                 |
| --output         | "./output"       | Path to directory to store output files in           |                 |
| --threads        | 1                | Number of threads to use                             |                 |
| --qualityControl | "no"             | Perform quality control with FastQC ["no"/"yes"] |                 |
| --fastqc         | -                | Path to executable                                   |                 |
| --trimming       | "no"             | Perform trimming with Trimmomatic ["no"/"yes"]   |                 |
| --trimmomatic    | -                | Path to executable Trimmomatic                       |                 |
| --trimOptions             | -           | Trimming options Trimmomatic is using E.g. "TRAILING:20 MINLEN:75"                         |   |
||||
| --useCellranger           | "yes"       | Use the CellRanger branch ["no/"yes"]                                                      | CellRanger |
| --cellranger              | -           | Path to executable CellRanger                                                              |   |
|--cellrangerTranscriptome| -           | Path to reference data set required for CellRanger Available on 10x Genomics download page |   |
| --CRoptions               | -           | Alternative parameters for CellRanger                                                      |   |
||||
| --useSTARsolo             | "no"        | Use STARsolo branch ["no"/"yes"]                                                           | STARsolo |
| --STARwhitelist           | -           | Path to barcode whitelist for STARsolo                                                     |   |
||||
| --genome                  | -           | Path to reference genome fasta file                                                        | STARsolo UMI-tools|
| --annotation              | -           | Path to reference genome annotation file                                                   |   |
| --star                    | -           | Path to executable STAR                                                                    |   |
| --STARoptions             | -           | Alternative parameters for STAR. See STAR manual                                           |   |
| --index                   | "no"        | Genome index is already generated ["yes"/"no"]                                             |   |
| --indicesDir              | "./indices" | Directory path to store genome indices If --index "yes", specify path with file prefix     |   |
| --read                    | "R2"        | Read containing cDNA sequence ["R1"/"R2"]                                                  |   |
| --barcode                 | "R1"        | Read containing barcode and UMI ["R1"/"R2"]                                                |   |
| --useLanes                | "all"       | Sequencing lanes to use For example, use only lane 1,2,3: Enter "1,2,3"                    |   |
| --genWhitelist            | "yes"       | Generate barcode whitelist ["no"/"yes"]                                                    |   |
| ---umi-tools              | -           | Path to executable UMI-tools                                                               |   |
| --useUMItools             | "no"        | Use the UMI-tools branch                                                                   |   |
| --samtools                | -           | Path to executable SAMtools                                                                |   |
| --featureCounts           | -           | Path to executable featureCounts                                                           |   |
| --UMITOOLSwhitelist       | -           | Path to barcode whitelist for UMI-tools branch                                             |   |
||||
| --nGenes                  | 100         | Minimum number of expressed genes of a cell                                                |   |
| --nUMIs                   | 125         | Minimum number of UMI counts of a cell                                                     |   |
| --MAD                     | 5           | nmads used scater's isOutlier function                                                     |   |
| --thresholdMT             | 1           | Maximal fraction of total UMI counts coming from MT genes                                  |   |
| --filterGenes             | 0.001       | Remove sparsely expressed genes Fraction of cells a gene has to be expressed               |   |
| --normalize               | "yes"       | Log normalize gene-barcode matrix ["yes"/"no"]          
 
