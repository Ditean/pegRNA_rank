# Modified atgRNA_rank
This is a modified version of the deep learning MLP model to rank atgRNA sequence by the [Abudayyeh-Gootenberg Lab](https://github.com/abugoot-lab/pegRNA_rank). Please refer to and cite their paper for information about atgRNA. This fork is more intended for personal use rather than development for general use.

Yarnall, M. T., Ioannidi, E. I., Schmitt-Ulms, C., Krajeski, R. N., Lim, J., Villiger, L., ... & Gootenberg, J. S. (2023). 
[Drag-and-drop genome insertion of large sequences without double-strand DNA cleavage using CRISPR-directed integrases](https://www.nature.com/articles/s41587-022-01527-4). Nature Biotechnology, 41(4), 500-512.

## Getting started
pegRNA_rank has been modified to accept the use of amplicons and inserts directly from the command-line. This provides a convenient way to run and increase the throughput of the original pegRNA_rank.py code.

### ***Providing a single amplicon and insertion sequence***
pegRNA_rank now accepts user inputs for amplicons (--amplicon) and inserts (--insert). The current version of pegRNA_rank can only accept a single amplicon and insert per execution.

```
python3 rank.py --amplicon <amplicond sequence> --insert <insert sequence>
```

### ***Providing a list of amplicons and insertion sequences***
pegRNA_rank has been modified with the additon of auto_pegRNA_rank.sh. This bash script allows users to provide a separate file containing all of the amplicons that you want to test and the insertion sequences that you want to screen. Each sequence needs to be on a new line, and the amplicons and insertions must be in different files.

auto_pegRNA_rank.sh currently requires both a amplicon and insertion file to be presented to function.

```
auto_pegRNA_rank.sh --amplicon <amplicon file> --insert <insertion file>
```

## Installing

### Prerequisites

The current pipeline has been tested on the following package versions
* **Python 3.6** or later is **required**
* **numpy 1.26.2**
* **pandas 2.1.3**
* **torch 2.2.1**

### pip installation

pegRNA_rank instructions: Make sure you have Python>=3.6 and PyTorch>=1.5 installed. Then, install dependencies with:
```
pip install -r requirements.txt
```

If you need to type python3 to use python and you are having problems with dependencies not being available, please consider using pip3
```
pip3 install -r requirements.txt
```

## Alternative method of running pegRNA_rank
When you only have a single insert or amplicon that you want to screen across a much larger library, consider using the the following command to read and loop each line within a text file.

### Example 1
```
while IFS= read -r amplicon; do python3 rank.py --amplicon "$amplicon" --insert ATGATGCCATGCTGATCGATGCTGAC; done < sequences.txt
```
In the following example, the command will loop across each line in sequences.txt ($amplicon) and use the insertion sequence _ATGATGCCATGCTGATCGATGCTGAC_.

***
### Example 2
```
while IFS= read -r insert; do python3 rank.py --amplicon ATGATGCCATGCTGATCGATGCTGAC --insert "$insert"; done < insertion.txt
```
In the following example, the command will use the amplicon sequence _ATGATGCCATGCTGATCGATGCTGAC_ and loop line by line the file insertion.txt for the insertion sequences to use.
