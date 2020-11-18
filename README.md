# CStreet: a computed <ins>C</ins>ell <ins>S</ins>tates <ins>tr</ins>ajectory inf<ins>e</ins>r<ins>e</ins>nce method for <ins>t</ins>ime-series single-cell RNA-seq data

**| [Overview](#overview) | [Installation](#installation) | [Quick Start](#quick-start) | [Parameter Details](#parameter-details) | [Run CStreet in python interface](#run-cstreet-in-python-interface) | [Citation](#citation) |**


![Figure1](https://github.com/yw-Hua/MarkdownPicture/raw/master/CStreet/Fig1.png)

## Overview

CStreet is a cell state trajectory inference method for time-series single-cell RNA-seq data. It is written in Python (Python 3.6 or higher) and is available as a command line tool and a Python library to meet the needs of different users.

CStreet uses time-series information to construct the k-nearest neighbors connections within and between time points. Then, CStreet calculates the connection probabilities of cell states and visualizes the trajectory, which may include multiple starting points and paths, using a force-directed layout method.

## Installation

CStreet has been packaged and uploaded to [PyPI](https://pypi.org/project/cstreet/). Before your installation, ensure that you have pip available. pip3 is the package installer for Python. If you do not have pip3 on your machine, click [here](https://pip.pypa.io/en/stable/) to install it. Then, CStreet and its relevant packages can be installed using a single command.

   ```shell
   $ pip3 install cstreet 
   ```


Type the command below to check whether CStreet has been installed successfully.

   ```shell
   $ CStreet -h
   ```

## Quick Start

### Step 1. CStreet installation following the above tutorial.

### Step 2. Input preparation

CStreet utilizes time-series expression levels in tab-delimited format or AnnData format as input. A small test dataset containing normalized expression levels at three time points can be downloaded here.

The cell state information can be generated using the built-in clustering function of CStreet or input by the user. The state information for the test dataset can be downloaded [here]().

### Step 3. Operation of CStreet

**Commandline**:

   ```shell
   $ CStreet -i ExpressionMatrix_t1.txt ExpressionMatrix_t1.txt ExpressionMatrix_t1.txt -s CellStates_t1.txt CellStates_t2.txt CellStates_t3.txt -n ProjectName
   ```

### Step 4. Output

The contents of the output directory in tree format will be displayed as described below, including the clustered cell state information if it is not provided by users, the connection probabilities of the cell states and a visualization of the inferred cell state trajectories.

```shell
PATH/ProjectName
├── ProjectName_CStreetTopology.pdf
├── SupplementaryFigures
│   ├── ProjectName_t1_ForceDirectedLayout.pdf
│   ├── ProjectName_t1_LouvainUMAPClustering.pdf
│   ├── ProjectName_t1_LouvainUMAPClusteringCoordinates.txt
│   ├── ProjectName_t2_ForceDirectedLayout.pdf
│   ├── ProjectName_t2_LouvainUMAPClustering.pdf
│   ├── ProjectName_t2_LouvainUMAPClusteringCoordinates.txt
│   └── ...
└── SupplementaryResults
    ├── ProjectName_BetweenTimePoints_CellStatesConnectionProbabilities.txt
    ├── ProjectName_t1_FilteredCellInfo.txt
    ├── ProjectName_t1_FilteredGeneInfo.txt
    ├── ProjectName_t1_CellStatesConnectionProbabilities.txt
    ├── ProjectName_t2_FilteredCellInfo.txt
    ├── ProjectName_t2_FilteredGeneInfo.txt
    ├── ProjectName_t2_CellStatesConnectionProbabilities.txt
    └── ...
```


   ![tiny_data_result.pdf](https://github.com/yw-Hua/MarkdownPicture/raw/master/CStreet/tiny_data_result.png)

## Parameter Details

The parameter details of CStreet are as follows:

-i --Input_ExpMatrix 
This indicates expression matrixes, which contain the time-series expression level as read counts or normalized values in tab-delimited format. 

-T --Input_CellonCol
This determines whether the expression level of one cell is displayed as a column in the expression matrixes. For example, '-T' indicates that the gene expression levels of one cell are displayed in a column and the expression levels of one gene across all cells are displayed in a row. 

-o --Output_Dir
The project directory, which is used to save all output files. DEFAULT: "./".

-n --Output_Name
The project name, which is used to generate output file names as a prefix. DEFAULT: "CStreet".

-s --Input_CellStates
An optional parameter that uses CStreet's built-in dimensionality reduction and clustering methods to perform clustering without knowing the cell states (DEFAULT) or accepts the user’s input. The input files should contain the cell state information sharing the same cell ID in the expression matrixes in tab-delimited format.

--CellClusterParam_PCAn
The number of principal components to use, which is enabled only if the cell state information is not provided. It can be set from 1 to the minimum dimension size of the expression matrixes. DEFAULT: 10.

--CellClusterParam_k
The number of nearest neighbors to be searched, which is enabled only if the cell state information is not provided. It should be in the range of 2 to 100 in general. DEFAULT: 15. 

--CellClusterParam_Resolution
The resolution of the Louvain algorithm, which is enabled only if the cell state information is not provided. A higher resolution means that more and smaller clusters are found. DEFAULT: 1.0. 

--Switch_DeadCellFilter
The switch for the dead cell filter, which filters cell outliers based on the count percent of mitochondrial genes. DEFAULT: "ON".

--Threshold_MitoPercent
The maximum count percent of mitochondrial genes needed for a cell to pass filtering, which is enabled only if '--Switch_DeadCellFilter' is "ON". DEFAULT: 0.2.

--Switch_LowCellNumGeneFilter
The switch for the low cell number gene filter, which retains genes that are expressed in at least a certain number of cells. DEFAULT: "ON".

--Threshold_LowCellNum
The minimum number of cells expressed that is required for a gene to pass filtering, which is enabled only if '-- Switch_LowCellNumGeneFilter' is "ON". DEFAULT: 3.

--Switch_LowGeneCellsFilter
The switch for the low gene number cell filter, which retains cells with at least a certain number of genes expressed. DEFAULT: "ON".

--Threshold_LowGeneNum
The minimum number of genes expressed that is required for a cell to pass filtering, which is enabled only if '-- Switch_LowGeneCellsFilter' is "ON". DEFAULT: 200.

--Switch_Normalize
The switch to enable total read count normalization for each cell. DEFAULT: "NO".

--Threshold_NormalizeBase
The base of normalization, which is enabled only if '--Switch_Normalize' is "ON". If the DEFAULT is chosen, it is CPM normalization. DEFAULT: 1e6.

--Switch_LogTransform
The switch to logarithmize the expression matrix. DEFAULT: "NO".

--WithinTimePointParam_PCAn
The number of principal components to use within a timepoint. It can be set from 1 to the minimum dimension size of the expression matrixes. DEFAULT: 10.

--WithinTimePointParam_k
The number of nearest neighbors to be searched within a timepoint. It should be in the range of 2 to 100 in general. DEFAULT: 15.

--BetweenTimePointParam_PCAn
The number of principal components to use between timepoints. It can be set from 1 to the minimum dimension size of the expression matrixes. DEFAULT: 10. 

--BetweenTimePointParam_k
The number of nearest neighbors to be searched between timepoints. It should be in the range of 2 to 100 in general. DEFAULT: 15.

--ProbParam_SamplingSize
The number of repeated sampling trials used to estimate the connection probability. DEFAULT: 5.

--ProbParam_RandomSeed
The random seed of repeated sampling, which makes the connection probability reproducible. DEFAULT: 0.

--FigureParam_FigureSize
The size of the resulting figure. For example: '--FigureParam_FigureSize 6 7' means width is 6 and height is 7. DEFAULT: 6 7 .)

--FigureParam_LabelBoxWidth
The width of the label box in the result figure. For example, '--FigureParam_LabelBoxWidth 10' means that 10 characters will be shown in the label box of the resulting figure. DEFAULT: 10.

--Threshold_MaxOutDegree
The maximum number of outdegrees for each cell state that is displayed, which will only be used for visualization. DEFAULT: 10.

--Threshold_MinCellNumofStates
The minimum cell number of each cell state that is displayed, which will only be used for visualization. DEFAULT: 50.



## Run CStreet in python interface

CStreet can also be used step by step in the Python interface and easily integrated into custom scripts. [Here]() is a tutorial written using Jupyter Notebook.

## Citation

   > 

