# CStreet: a computed <ins>C</ins>ell <ins>S</ins>tates <ins>tr</ins>ajectory inf<ins>e</ins>r<ins>e</ins>nce method for <ins>t</ins>ime-series single-cell RNA-seq data

**| [Overview](#overview) | [Installation](#installation) | [Quick Start](#quick-start) | [Parameter Details](#parameter-details) | [Run CStreet in python interface](#run-cstreet-in-python-interface) | [Citation](#citation) |**


![Figure1](https://github.com/yw-Hua/MarkdownPicture/raw/master/CStreet/Fig1.png)

## Overview

CStreet is a cell states trajectory inference method for time-series single-cell RNA-seq data. It is written in python (python 3.6 or higher) and is available as a commend line tool and a python library to meet the needs of different users.

CStreet takes advantage of time-series information to construct the *k*-nearest neighbors connections within and between time points. Then CStreet calculated the connection probabilities of cell states and visualized the trajectory which may include multiple starting points and paths using a force-directed layout method. 

## Installation

CStreet has been packaged and uploaded to [PyPI](https://pypi.org/project/cstreet/). Before your installation, you'll make sure you have pip available. The pip3 is the package installer for Python. If you don't have pip3 on your machine, try [click here](https://pip.pypa.io/en/stable/) to install it. Then CStreet and its relevant packages can be installed using one single commands as follows.

   ```shell
   $ pip3 install cstreet 
   ```


Type the following command to check whether CStreet has been installed successfully.

   ```shell
   $ CStreet -h
   ```

## Quick Start

**Input**: 

   - Expression data: Expression matrix containing the time-series expression level as reads counts or normalized values in tab delimited format, and anndata format are accepted as the input of CStreet. (For example: [ExpressionMatrix_t1.txt]() [ExpressionMatrix_t2.txt]() [ExpressionMatrix_t3.txt]())
   - Cell states info: The cell states information can be inputted by the user or generated using the internal clustering function of CStreet. (For example: [CellStates_t1.txt]() [CellStates_t2.txt]() [CellStates_t3.txt]())

**Commandline**:

   ```shell
   $ CStreet -i ExpressionMatrix_t1.txt ExpressionMatrix_t1.txt ExpressionMatrix_t1.txt -s CellStates_t1.txt CellStates_t2.txt CellStates_t3.txt -n ProjectName
   ```

**Output**: 

The contents of the output directory in tree format will be displayed as follows: 

```shell
PATH/ProjectName
├── cstreet_result.pdf
├── figures
│   ├── timepoint1_fa.pdf
│   ├── timepoint1_louvain_umap_cord.txt
│   ├── timepoint1_louvain_umap.pdf
│   ├── timepoint2_fa.pdf
│   ├── timepoint2_louvain_umap_cord.txt
│   ├── timepoint2_louvain_umap.pdf
│   ├── timepoint3_fa.pdf
│   ├── timepoint3_louvain_umap_cord.txt
│   └── timepoint3_louvain_umap.pdf
└── results
    ├── alltimepoint_link_cluster_graph.txt
    ├── timepoint1_cell_info.txt
    ├── timepoint1_gene_info.txt
    ├── timepoint1_inner_cluster_graph.txt
    ├── timepoint2_cell_info.txt
    ├── timepoint2_gene_info.txt
    ├── timepoint2_inner_cluster_graph.txt
    ├── timepoint3_cell_info.txt
    ├── timepoint3_gene_info.txt
    └── timepoint3_inner_cluster_graph.txt
```

   - The clustered cell states information if not provided by users. (timepoint*_cell_info.txt)
   - The connection probabilities of cell states. (alltimepoint_link_cluster_graph.txt, timepoint*_inner_cluster_graph.txt)
   - An visulization of inferred cell states trajectory. (cstreet_result.pdf)



   ![tiny_data_result.pdf](https://github.com/yw-Hua/MarkdownPicture/raw/master/CStreet/tiny_data_result.png)

## Parameter Details

```
usage: CStreet [-h] <-i ExpMatrix1 ExpMatrix2 ExpMatrix3 ...> [-s CellStates1 CellStates2 CellStates3 ...] [-n ProjectName] [-o OutputDir] [options]

CStreet is a cell states trajectory inference method for time-series single-cell RNA-seq data.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_EXPMATRIX [INPUT_EXPMATRIX ...], --Input_ExpMatrix INPUT_EXPMATRIX [INPUT_EXPMATRIX ...]
                        Expression matrixes, which will contain the time-
                        series expression level as reads counts or normalized
                        values in tab delimited format. For example: '-i
                        ExpressionMatrix_t1.txt ExpressionMatrix_t2.txt
                        ExpressionMatrix_t3.txt' means the input of 3
                        timepoints expression matrixes.
  --Input_CellonCol {True,False}, -T {True,False}
                        Whether the cells are arranged on rows or columns in
                        the expression matrixes. For example: '-T True' means
                        cells on columns and genes on rows. DEFAULT: False.
  --Output_Dir OUTPUT_DIR, -o OUTPUT_DIR
                        Project directory, which will be used to save all
                        output files. DEFAULT: "./"
  --Output_Name OUTPUT_NAME, -n OUTPUT_NAME
                        Project name, which will be used to generate output
                        file names. DEFAULT: "cstreet_project"
  --Input_CellStates [INPUT_CELLSTATES [INPUT_CELLSTATES ...]], -s [INPUT_CELLSTATES [INPUT_CELLSTATES ...]]
                        Cell states information, which must contain a columns
                        named "state" and the same cell ID with expression
                        matrixes in tab delimited format. Cell states
                        information can be inputted by the user or generated
                        by the internal clustering function of CStreet using
                        the Louvain algorithm. For example: 'CellStates_t1.txt
                        CellStates_t2.txt CellStates_t3.txt' means the cell
                        states information of 3 timepoints expression
                        matrixes.
  --CellClusterParam_PCAn CELLCLUSTERPARAM_PCAN
                        Number of principal components to use, which will be
                        enabled ONLY if cell states information is not
                        provided. It can be set to 1 - minimum dimension size
                        of expression matrixes. DEFAULT: 10
  --CellClusterParam_kNNn CELLCLUSTERPARAM_KNNN
                        Number of nearest neighbors to be searched, which will
                        be enabled ONLY if cell states information is not
                        provided. It should be in the range 2 to 100 in
                        general. DEFAULT: 15
  --CellClusterParam_Resolution CELLCLUSTERPARAM_RESOLUTION
                        Resolution of the Louvain algorithm, which will be
                        enabled ONLY if cell states information is not
                        provided. Higher resolution means finding more and
                        smaller clusters. DEFAULT: 1.0
  --Switch_DeadCellFilter {ON,OFF}
                        A switch of dead cell filter, which filter cell
                        outliers based on counts percent of Mitochondrial
                        gene. DEFAULT: "ON"
  --Threshold_MitoPercent THRESHOLD_MITOPERCENT
                        Maximum counts percent of Mitochondrial gene for a
                        cell to pass filtering, which will be enabled ONLY if
                        '--Switch_DeadCellFilter' is "ON". DEFAULT: 0.2
  --Switch_LowCellNumGeneFilter {ON,OFF}
                        A switch of low cell number gene filter, which keep
                        genes that are expressed in at least a number of
                        cells. DEFAULT: "ON"
  --Threshold_LowCellNum THRESHOLD_LOWCELLNUM
                        Minimum number of cells expressed required for a gene
                        to pass filtering, which will be enabled ONLY if '--
                        Switch_LowCellNumGeneFilter' is "ON". DEFAULT: 3
  --Switch_LowGeneCellsFilter {ON,OFF}
                        A switch of low gene number cell filter, which keep
                        cells with at least a number of genes expressed.
                        DEFAULT: "ON"
  --Threshold_LowGeneNum THRESHOLD_LOWGENENUM
                        Minimum number of genes expressed required for a cell
                        to pass filtering, which will be enabled ONLY if '--
                        Switch_LowGeneCellsFilter' is "ON". DEFAULT: 200
  --Switch_Normalize {ON,OFF}
                        A switch of total read count normalization for per
                        cell. DEFAULT: "NO"
  --Threshold_NormalizeBase THRESHOLD_NORMALIZEBASE
                        Normalize Base of normalization, which will be enabled
                        ONLY if '--Switch_Normalize' is "ON". If choosing
                        DEFAULT, it is CPM normalization. DEFAULT: 1e6
  --Switch_LogTransform {ON,OFF}
                        A switch of logarithmizing the expression matrix.
                        DEFAULT: "NO"
  --WithinTimePointParam_PCAn WITHINTIMEPOINTPARAM_PCAN
                        Number of principal components to use, which will be
                        used within a timepoint. It can be set to 1 - minimum
                        dimension size of expression matrixes. DEFAULT: 10
  --WithinTimePointParam_kNNn WITHINTIMEPOINTPARAM_KNNN
                        Number of nearest neighbors to be searched, which will
                        be used within a timepoint. It should be in the range
                        2 to 100 in general. DEFAULT: 15
  --BetweenTimePointParam_PCAn BETWEENTIMEPOINTPARAM_PCAN
                        Number of principal components to use, which will be
                        used between timepoints. It can be set to 1 - minimum
                        dimension size of expression matrixes. DEFAULT: 10
  --BetweenTimePointParam_kNNn BETWEENTIMEPOINTPARAM_KNNN
                        Number of nearest neighbors to be searched, which will
                        be used between timepoints. It should be in the range
                        2 to 100 in general. DEFAULT: 15
  --ProbParam_SamplingSize PROBPARAM_SAMPLINGSIZE
                        Number of repeated sampling trials to estimate the
                        connection probability. DEFAULT: 5
  --ProbParam_RandomSeed PROBPARAM_RANDOMSEED
                        Random seed of repeated sampling, which will make the
                        connection probability is reproducible. DEFAULT: 0
  --FigureParam_FigureSize FIGUREPARAM_FIGURESIZE
                        Figure size of the result figure. Format is (width,
                        height). DEFAULT:(6,7)
  --FigureParam_LabelBoxWidth FIGUREPARAM_LABELBOXWIDTH
                        Width of the label box in the result figure. For
                        example: '--FigureParam_LabelBoxWidth 10' means 10
                        characters will be showed in label box of result
                        figure. DEFAULT: 10
  --Threshold_MaxOutDegree THRESHOLD_MAXOUTDEGREE
                        Maximum number of outdegree for each cell state will
                        be displayed, which will ONLY be used for
                        visualization. DEFAULT: 10
  --Threshold_MinCellNumofStates THRESHOLD_MINCELLNUMOFSTATES
                        Minimum cell number of each cell state will be
                        displayed, which will ONLY be used for visualization.
                        DEFAULT: 50

```



## Run CStreet in python interface

CStreet can run directly or step by step in [Jupyter Notebook](https://jupyter.org/). 
The tutorial  is [here]().

## Citation

   > 

