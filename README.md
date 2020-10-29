# CStreet Overview
CStreet is a cell states trajectory construction method for time-series single-cell RNA-seq data. It is written in python (python 3.6 or higher) and is available as a commend line tool and a python library to meet the needs of different users.

[Figure1]

CStreet takes advantage of time-series information to construct the connections of *k*-nearest neighbors within and between time points. Then CStreet calculated the connection probabilities of cell states and visualized the trajectory which may include multiple starting points and paths using a force-directed layout method. 

# Installation
CStreet has been packaged and upload to PyPI. CStreet and its relevant packages can be installed using one single commands as follows.
   ```shell
   $ pip3 install cstreet # pip3 is the package installer for Python. If you don't have pip3 on your machine, try install it [https://pip.pypa.io/en/stable/].
   ```
Type the following command to check whether CStreet has been installed successfully.
   ```shell
   $ CStreet -h
   ```

# Quick Start
**Input**: 
- Expression data: Expression matrix containing the time-series expression level as reads counts or normalized values in tab delimited format, and anndata format are accepted as the input of CStreet.
- Cell states info: The cell states information can be inputted by the user or generated using the internal clustering function of CStreet.
**Output**: 
- An visulization of inferred cell states trajectory.
- The clustered cell states information if not provided by users.
- The connection probabilities of cell states.

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

   

# Example Results

[SFig1]

**An example of inferenced cell trajectory**:

![results.png](https://github.com/yw-Hua/MarkdownPicture/blob/master/CStreet/results2.png?raw=true)

