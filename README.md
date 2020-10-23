# CStreet Overview
CStreet is a python script (python 3.6 or higher) for cell states trajectory construction by using *k*-nearest neighbors graph algorithm for time-series single-cell RNA-seq data. It is a **developmental version**.

# Installation

1. Install CStreet by `pip3`

   CStreet can be installed directly by using  `pip3` commands :

   ```shell
   $ pip3 install cstreet
   ```




# Quick Start

**Input file**: Only expression matrix containing the time-series expression level as reads counts or normalized values for this developmental version.

**Output file**: An inferenced cell states trajectory.

1. Add new time-series single cell RNA-seq data.

   ```python
   from cstreet import *
   import pandas as pd
   # Read single cell data as DataFrame
   data_t1=pd.read_table('data_t1.txt',header=0, sep="\t",index_col=0) 
   data_t2=pd.read_table('data_t2.txt',header=0, sep="\t",index_col=0)
   data_t3=pd.read_table('data_t3.txt',header=0, sep="\t",index_col=0)
   # Create a new CStreet object
   cdata=CStreetData()
   # add data into CStreet object
   cdata.add_new_timepoint_scdata(data_t1)
   cdata.add_new_timepoint_scdata(data_t2)
   cdata.add_new_timepoint_scdata(data_t3)
   ```
   
2. Customize parameters.

   ```python
   #Step0:basic parameters
   cdata.params.output_dir="./"
   cdata.params.output_name="cstreet_project"
   
   
   #Step1:cell cluster
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
   
   #Step5:get graph
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

![results.png](https://github.com/yw-Hua/MarkdownPicture/blob/master/CStreet/results2.png?raw=true)

