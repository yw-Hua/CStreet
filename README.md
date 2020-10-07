# CStreet Overview

CStreet is a python script (python 3.6 or higher) for cell states trajectory construction by using *k*-nearest neighbors graph algorithm for time-series single-cell RNA-seq data. It is a **developmental version**.

# Installation

1. Prepare required packages
   CStreet depends on a number of `python3` packages available on pypi and all dependencies can be installed using  `pip3` commands :

   ```shell
   $ pip3 install scanpy
   $ pip3 install anndata
   $ pip3 install networkx
   $ pip3 install fa2
   $ pip3 install retrying
   ```

2. Download CStreet from github
   CStreet can be download using `git` command:

   ```shell
   $ cd /PATH/ # here you can replace "/PATH/" with any location you want
   $ git clone git://github.com/TongjiZhanglab/CStreet.git
   ```

3. Import the main class

   ```python
   import sys
   sys.path.append("/PATH/CStreet/") # here you should replace "/PATH/" with the location where CStreet has been installed at
   from cstreet import *
   ```

   



# Quick Start

**Input file**: Only expression matrix containing the time-series expression level as reads counts or normalized values for this developmental version.

**Output file**: An inferenced cell states trajectory.

1. Add new time-series single cell RNA-seq data.

   ```python
   import numpy as np
   import pandas as pd
   # Read single cell data as DataFrame
   data_t1=pd.read_table('data_t1.txt',header=0, sep="\t",index_col=0) 
   data_t1=pd.read_table('data_t2.txt',header=0, sep="\t",index_col=0)
   data_t2=pd.read_table('data_t3.txt',header=0, sep="\t",index_col=0)
   data_t3=pd.read_table('data_t4.txt',header=0, sep="\t",index_col=0)
   # Create a new CStreet object
   cdata=CStreetData()
   # add data into CStreet object
   cdata.add_new_timepoint_scdata(data_t1)
   cdata.add_new_timepoint_scdata(data_t2)
   cdata.add_new_timepoint_scdata(data_t3)
   ```

2. Customize parameters.

   ```python
   #Step1:cell_cluster
   cdata.params.cell_cluster_pca_n=10
   cdata.params.cell_cluster_knn_n=15
   cdata.params.cell_cluster_resolution=0.1
   
   #Step2:gene and cell filter
   cdata.params.filter_dead_cell=True
   cdata.params.percent_mito_cutoff=0.2
   cdata.params.filter_lowcell_gene=True
   cdata.params.min_cells=3
   cdata.params.filter_lowgene_cells=True
   cdata.params.min_genes=200
   
   #Step3:normalize
   cdata.params.normalize=True
   cdata.params.normalize_base=10000
   cdata.params.log_transform=True
   
   #Step4:get HVG
   cdata.params.highly_variable_genes=False
   
   #Step5:get_graph
   cdata.params.inner_graph_pca_n=10
   cdata.params.inner_graph_knn_n=15
   cdata.params.link_graph_pca_n=10
   cdata.params.link_graph_knn_n=15
   cdata.params.max_outgoing=10
   cdata.params.min_score=0.1
   cdata.params.min_cell_number=50
   ```

3. Run CStreet

   ```python
   cdata.run_cstreet()
   ```

   

# Result

**An example of inferenced cell trajectory**:

![results.png](https://github.com/yw-Hua/MarkdownPicture/blob/master/CStreet/results.png?raw=true)
