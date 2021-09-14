# PeakCNV

## A Multi-Feature Ranking Algorithm-based tool for Genome-wide Copy Number Variation-association study 

PeakCNV consists of three steps: (1) Building CNVRs; (2) Cluster CNVRs; (3) Selection.  PeakCNV first computes the genome-wide probability of each base pair to be deleted or duplicated in case and control samples using the Fisher exact test resulting in identification of statistically significant CNVRs. Then, PeakCNV groups statistically significant CNVRs into different clusters based on two combined features. The first feature (importance) shows the number of samples in each region after removing the common samples between each two CNVRs. The second one is the distance between CNVRs. Finally, the best CNVR(s) from each cluster will be report as the representative of the cluster



## Table of Contents

- [Installation](#installation)
- [PeakCNV analysis workflow](#peakCNV-analysis-workflow)
- [Arguments](#arguments)
- [Usage example](#Usage-example)
- [Reference](#reference)
- [Author Info](#author-info)
- [License](#license)


## Installation

MutSpot runs on R (requires at least 3.2.0). Install the package from Github using the following R commands.

```
install.packages("devtools")
library(devtools)
install_github("mahdieh1/PeakCNV")
```
Alternatively, the package may downloaded from Github and installed in R:
```
# Clone/download MutSpot into the current working dirctory with the following command: git clone https://github.com/mahdieh1/PeakCNV.git
```


## PeakCNV analysis workflow

The full PeakCNV workflow includes the following 3 steps:
1. Building CNVRs
2. Clustering process
3. Selection process

By default, the PeakCNV() function runs the entire workflow. However, it is possible to run specific steps of the workflow by specifiying the run.to parameter (see full documentation).


## Arguments

| Parameter | Type | Description | Default |
| :---: | :---: | :---: | :---: |
| run.to | Numeric | Steps to run | 1,2,3 |
| working.dir	| String | Working directory | current working directory


## Usage example
By default, PeakCNV runs in the current working directory unless specified by the user. To run PeakCNV on genome assembly hg18 instead of hg19, please change test.genome in the working directory. By default, results will be saved in the working directory. In clustering step, for each chromosome, PeakCNV asks you the eps value based on the k nearest neighbors(knn) plot. The optimal value is an elbow, where a sharp change in the distance occurs.
```
library("PeakCNV")
```
Download the test data sets from https://github.com/mahdieh1/PeakCNV/tree/main/test-data into your working directory.

#### Input files ####
1. Case CNVs:

| Chr | Start | End | Sample-Id |
| :---: | :---: | :---: | :---: |
| 1 | 6742281 | 6742903 | SP7890 |

2. Control CNVs

| Chr | Start | End | Sample-Id |
| :---: | :---: | :---: | :---: |
| 1 | 6742281 | 6742903 | as999 |

#### Output files: ####

1. CNVRs

| Chr | Start | End | #Sample-Id | #case | #control | p-value |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | 6742281 | 6742903 | SP7890 | 15 | 30 | 0.02 |

2. clustered CNVRs

| Chr | Start | End | #case | Cluster-NO |
| :---: | :---: | :---: | :---: | :---: |  
| 1 | 6742281 | 6742903 | 15 | 1 |

3. Selected CNVRs

| Score | #chr | Start | End | #case | Cluster-NO | 
| :---: | :---: | :---: | :---: | :---: | :---: | 
| 56 | 21 | 6742281 | 6742903 | 15 |	225.86 | 0 |

## Reference
```
Please consider citing the follow paper when you use this code.
  Title={},
  Authors={}
}
```
## License

This package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

A copy of the GNU General Public License, version 3, is available at https://www.r-project.org/Licenses/GPL-3

## Author Info

In case of queries, please email: mahdieh.labani@students.mq.edu.au






