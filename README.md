[![DOI](https://zenodo.org/badge/97020646.svg)](https://zenodo.org/badge/latestdoi/97020646)

# SalmID
Rapid tool to check taxonomic ID of single isolate samples. Currently only IDs Salmonella species and subspecies, and some common contaminants (Listeria, Escherichia).

## Requirements:
Python 3

## Installation:
The easy way with homebrew ([Linux](http://linuxbrew.sh/) or [MacOS](https://brew.sh/)):
```
brew install brewsci/bio/salmid
```
Big thanks to [Torsten Seemann](https://tseemann.github.io/) for including this in homebrew!

Alernatively git clone to your machine:
```
git clone --recursive https://github.com/hcdenbakker/SalmID.git
```

Make SalmID executable:
```
cd SalmID
```

```
chmod +x SalmID.py
```


Add the SalmID folder to your path

To execute:
```
SalmID.py -e .fastq.gz
```
This will perform a SalmID run on all fastq.gz files in the current directory.
```
SalmID.py -i your_fastq_gz.fastq.gz
```
This will perform a SalmID run on an individual file (i.e., your_fastq_gz.fastq.gz)
```
SalmID.py -d directory_with_data -e _1.fastq.gz
```
This will perform a SalmID run on all files in directory 'directory_with_data' with extension '_1.fastq.gz'

## Todo's and thoughts for future releases:
- Provide coverage estimates for genomes in sample based on kmer frequencies
- Write code to use SalmID on long read (minion, pacbio) platforms
