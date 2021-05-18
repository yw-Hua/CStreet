#!python3.6
#v1.1.1
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy
import scipy.optimize as op
from sklearn.decomposition import PCA 
from sklearn.neighbors import NearestNeighbors
from scipy.stats import zscore
import networkx as nx
from fa2 import ForceAtlas2
import matplotlib.pyplot as plt
import matplotlib as mpl
from functools import wraps
import time
from datetime import datetime
import traceback
from retrying import retry
import warnings
import random
from scipy.stats import t
warnings.filterwarnings("ignore")
from scipy.spatial.distance import correlation

#定义cstreet对象

class CStreetData(object):
    """docstring for CStreetData"""
    class params_object:
        """docstring for params_object"""
        __slots__=("__Output_Dir","__Output_Name","__CellClusterParam_PCAn","__CellClusterParam_k","__CellClusterParam_Resolution",
            "__Switch_DeadCellFilter","__Threshold_MitoPercent","__Switch_LowCellNumGeneFilter","__Threshold_LowCellNum","__Switch_LowGeneCellsFilter",
            "__Threshold_LowGeneNum","__Switch_Normalize","__Threshold_NormalizeBase","__Switch_LogTransform",
            "__WithinTimePointParam_PCAn","__WithinTimePointParam_k","__BetweenTimePointParam_PCAn","__BetweenTimePointParam_k",
            "__Threshold_MaxOutDegree","__Threshold_MinCellNumofStates","__ProbParam_RandomSeed","__ProbParam_SamplingSize",
            "__FigureParam_FigureSize","__FigureParam_LabelBoxWidth","__Threshold_MinProbability","__KNNParam_metric")
        
        def __init__(self):
            #Step0:basic params#
            self.__Output_Dir="./"
            self.__Output_Name="CStreet"
            #Step1:cell_cluster# 
            self.__CellClusterParam_PCAn=10
            self.__CellClusterParam_k=15
            self.__CellClusterParam_Resolution=0.1

            #Step2:gene and cell filter#
            self.__Switch_DeadCellFilter=True
            self.__Threshold_MitoPercent=0.2
            
            self.__Switch_LowCellNumGeneFilter=True
            self.__Threshold_LowCellNum=3
            
            self.__Switch_LowGeneCellsFilter=True
            self.__Threshold_LowGeneNum=200
            
            #Step3:normalize#
            self.__Switch_Normalize=True
            self.__Threshold_NormalizeBase=1000000
            self.__Switch_LogTransform=True

            #Step5:get_graph#
            self.__KNNParam_metric="euclidean"
            self.__WithinTimePointParam_PCAn=10
            self.__WithinTimePointParam_k=15

            self.__BetweenTimePointParam_PCAn=10
            self.__BetweenTimePointParam_k=15

            #Step6: calculate probability
            self.__ProbParam_RandomSeed=0
            self.__ProbParam_SamplingSize=5
            #Step7:plot graph#
            self.__FigureParam_FigureSize=(6,7)
            self.__FigureParam_LabelBoxWidth=10
            self.__Threshold_MinProbability="OTSU"
            self.__Threshold_MaxOutDegree=10
            self.__Threshold_MinCellNumofStates=0
        
        def __str__(self):
            s=""
            s+=f"\n#Step0:basic params# \n"
            s+=f"Output_Dir={self.__Output_Dir}\n"
            s+=f"Output_Name={self.__Output_Name}\n"            
            
            s+=f"\n#Step1:cell_cluster# \n"
            s+=f"CellClusterParam_PCAn={self.__CellClusterParam_PCAn}\n"
            s+=f"CellClusterParam_k={self.__CellClusterParam_k}\n"
            s+=f"CellClusterParam_Resolution={self.__CellClusterParam_Resolution}\n"

            s+=f"\n#Step2:gene and cell filter#\n"
            s+=f"Switch_DeadCellFilter={self.__Switch_DeadCellFilter}\n"
            s+=f"Threshold_MitoPercent={self.__Threshold_MitoPercent}\n"

            s+=f"Switch_LowCellNumGeneFilter={self.__Switch_LowCellNumGeneFilter}\n"
            s+=f"Threshold_LowCellNum={self.__Threshold_LowCellNum}\n"

            s+=f"Switch_LowGeneCellsFilter={self.__Switch_LowGeneCellsFilter}\n"
            s+=f"Threshold_LowGeneNum={self.__Threshold_LowGeneNum}\n"

            s+=f"\n#Step3:normalize#\n"
            s+=f"Switch_Normalize={self.__Switch_Normalize}\n"
            s+=f"Threshold_NormalizeBase={self.__Threshold_NormalizeBase}\n"
            s+=f"Switch_LogTransform={self.__Switch_LogTransform}\n"

            s+=f"\n#Step4:get_graph#\n"
            s+=f"KNNParam_metric={self.__KNNParam_metric}\n"
            s+=f"WithinTimePointParam_PCAn={self.__WithinTimePointParam_PCAn}\n"
            s+=f"WithinTimePointParam_k={self.__WithinTimePointParam_k}\n"

            s+=f"BetweenTimePointParam_PCAn={self.__BetweenTimePointParam_PCAn}\n"
            s+=f"BetweenTimePointParam_k={self.__BetweenTimePointParam_k}\n"
            
            s+=f"\n#Step5: calculate probability#\n"
            s+=f"ProbParam_RandomSeed={self.__ProbParam_RandomSeed}\n"
            s+=f"ProbParam_SamplingSize={self.__ProbParam_SamplingSize}\n"
            
            s+=f"\n#Step6:plot graph#\n"
            s+=f"FigureParam_FigureSize={self.__FigureParam_FigureSize}\n"
            s+=f"FigureParam_LabelBoxWidth={self.__FigureParam_LabelBoxWidth}\n"
            s+=f"Threshold_MinProbability={self.__Threshold_MinProbability}\n"
            s+=f"Threshold_MaxOutDegree={self.__Threshold_MaxOutDegree}\n"
            s+=f"Threshold_MinCellNumofStates={self.__Threshold_MinCellNumofStates}\n"
            return s
        def __repr__(self):
            return self.__str__()
        @property
        def Output_Dir(self):
            return self.__Output_Dir
        @Output_Dir.setter
        def Output_Dir(self,value):
            if not isinstance(value,str):
                raise ValueError('Output_Dir must be a string')
            if value[-1] != "/":
                value+="/"
            self.__Output_Dir = value

        @property
        def Output_Name(self):
            return self.__Output_Name
        @Output_Name.setter
        def Output_Name(self,value):
            if not isinstance(value,str):
                raise ValueError('Output_Name must be a string')
            self.__Output_Name = value

        @property
        def CellClusterParam_PCAn(self):
            return self.__CellClusterParam_PCAn
        @CellClusterParam_PCAn.setter
        def CellClusterParam_PCAn(self,value):
            if not isinstance(value,int):
                raise ValueError('CellClusterParam_PCAn must be an integer')
            if value <= 0:
                raise ValueError('CellClusterParam_PCAn must be bigger than 0')
            self.__CellClusterParam_PCAn = value

        @property
        def CellClusterParam_k(self):
            return self.__CellClusterParam_k
        @CellClusterParam_k.setter
        def CellClusterParam_k(self,value):
            if not isinstance(value,int):
                raise ValueError('CellClusterParam_k must be an integer')
            if value <= 0:
                raise ValueError('CellClusterParam_k must be bigger than 0')
            self.__CellClusterParam_k = value

        @property
        def CellClusterParam_Resolution(self):
            return self.__CellClusterParam_Resolution
        @CellClusterParam_Resolution.setter
        def CellClusterParam_Resolution(self,value):
            if not isinstance(value,(int,float)):
                raise ValueError('CellClusterParam_Resolution must be numeric')
            if value <= 0:
                raise ValueError('CellClusterParam_Resolution must be bigger than 0')        
            self.__CellClusterParam_Resolution = value

        @property
        def Switch_DeadCellFilter(self):
            return self.__Switch_DeadCellFilter
        @Switch_DeadCellFilter.setter
        def Switch_DeadCellFilter(self,value):
            if not isinstance(value,bool):
                raise ValueError('Switch_DeadCellFilter must be True or False')
            self.__Switch_DeadCellFilter = value

        @property
        def Threshold_MitoPercent(self):
            return self.__Threshold_MitoPercent
        @Threshold_MitoPercent.setter
        def Threshold_MitoPercent(self,value):
            if not isinstance(value,(int,float)):
                raise ValueError('Threshold_MitoPercent must be numeric')
            if value < 0 or value > 1:
                raise ValueError('Threshold_MitoPercent must be between 0.0 and 1.0') 
            self.__Threshold_MitoPercent = value

        @property
        def Switch_LowCellNumGeneFilter(self):
            return self.__Switch_LowCellNumGeneFilter
        @Switch_LowCellNumGeneFilter.setter
        def Switch_LowCellNumGeneFilter(self,value):
            if not isinstance(value,bool):
                raise ValueError('Switch_LowCellNumGeneFilter must be True or False')
            self.__Switch_LowCellNumGeneFilter = value

        @property
        def Threshold_LowCellNum(self):
            return self.__Threshold_LowCellNum
        @Threshold_LowCellNum.setter
        def Threshold_LowCellNum(self,value):
            if not isinstance(value,int):
                raise ValueError('Threshold_LowCellNum must be an integer')
            if value <= 0:
                raise ValueError('Threshold_LowCellNum must be bigger than 0')
            self.__Threshold_LowCellNum = value

        @property
        def Switch_LowGeneCellsFilter(self):
            return self.__Switch_LowGeneCellsFilter
        @Switch_LowGeneCellsFilter.setter
        def Switch_LowGeneCellsFilter(self,value):
            if not isinstance(value,bool):
                raise ValueError('Switch_LowGeneCellsFilter must be True or False')
            self.__Switch_LowGeneCellsFilter = value

        @property
        def Threshold_LowGeneNum(self):
            return self.__Threshold_LowGeneNum
        @Threshold_LowGeneNum.setter
        def Threshold_LowGeneNum(self,value):
            if not isinstance(value,int):
                raise ValueError('Threshold_LowGeneNum must be an integer')
            if value <= 0:
                raise ValueError('Threshold_LowGeneNum must be bigger than 0')
            self.__Threshold_LowGeneNum = value

        @property
        def Switch_Normalize(self):
            return self.__Switch_Normalize
        @Switch_Normalize.setter
        def Switch_Normalize(self,value):
            if not isinstance(value,bool):
                raise ValueError('Switch_Normalize must be True or False')
            self.__Switch_Normalize = value

        @property
        def Threshold_NormalizeBase(self):
            return self.__Threshold_NormalizeBase
        @Threshold_NormalizeBase.setter
        def Threshold_NormalizeBase(self,value):
            if not isinstance(value,int):
                raise ValueError('Threshold_NormalizeBase must be an integer')
            if value <= 0:
                raise ValueError('Threshold_NormalizeBase must be bigger than 0')
            self.__Threshold_NormalizeBase = value

        @property
        def Switch_LogTransform(self):
            return self.__Switch_LogTransform
        @Switch_LogTransform.setter
        def Switch_LogTransform(self,value):
            if not isinstance(value,bool):
                raise ValueError('Switch_LogTransform must be True or False')
            self.__Switch_LogTransform = value

        @property
        def WithinTimePointParam_PCAn(self):
            return self.__WithinTimePointParam_PCAn
        @WithinTimePointParam_PCAn.setter
        def WithinTimePointParam_PCAn(self,value):
            if not isinstance(value,int):
                raise ValueError('WithinTimePointParam_PCAn must be an integer')
            if value <= 0:
                raise ValueError('WithinTimePointParam_PCAn must be bigger than 0')
            self.__WithinTimePointParam_PCAn = value

        @property
        def WithinTimePointParam_k(self):
            return self.__WithinTimePointParam_k
        @WithinTimePointParam_k.setter
        def WithinTimePointParam_k(self,value):
            if not isinstance(value,int):
                raise ValueError('WithinTimePointParam_k must be an integer')
            if value <= 0:
                raise ValueError('WithinTimePointParam_k must be bigger than 0')
            self.__WithinTimePointParam_k = value

        @property
        def BetweenTimePointParam_PCAn(self):
            return self.__BetweenTimePointParam_PCAn
        @BetweenTimePointParam_PCAn.setter
        def BetweenTimePointParam_PCAn(self,value):
            if not isinstance(value,int):
                raise ValueError('BetweenTimePointParam_PCAn must be an integer')
            if value <= 0:
                raise ValueError('BetweenTimePointParam_PCAn must be bigger than 0')
            self.__BetweenTimePointParam_PCAn = value

        @property
        def BetweenTimePointParam_k(self):
            return self.__BetweenTimePointParam_k
        @BetweenTimePointParam_k.setter
        def BetweenTimePointParam_k(self,value):
            if not isinstance(value,int):
                raise ValueError('BetweenTimePointParam_k must be an integer')
            if value <= 0:
                raise ValueError('BetweenTimePointParam_k must be bigger than 0')
            self.__BetweenTimePointParam_k = value

        @property
        def ProbParam_SamplingSize(self):
            return self.__ProbParam_SamplingSize
        @ProbParam_SamplingSize.setter
        def ProbParam_SamplingSize(self,value):
            if not isinstance(value,int):
                raise ValueError('ProbParam_SamplingSize must be be an integer')
            if value <= 0:
                raise ValueError('ProbParam_SamplingSize must be bigger than 0')
            self.__ProbParam_SamplingSize = value

        @property
        def ProbParam_RandomSeed(self):
            return self.__ProbParam_RandomSeed
        @ProbParam_RandomSeed.setter
        def ProbParam_RandomSeed(self,value):
            if not isinstance(value,int):
                raise ValueError('ProbParam_RandomSeed must be be an integer')
            if value < 0:
                raise ValueError('ProbParam_RandomSeed must be bigger than 0')
            self.__ProbParam_RandomSeed = value

        @property
        def FigureParam_FigureSize(self):
            return self.__FigureParam_FigureSize
        @FigureParam_FigureSize.setter
        def FigureParam_FigureSize(self,value):
            if not isinstance(value,tuple):
                raise ValueError('FigureParam_FigureSize must be be a tuple')
            self.__FigureParam_FigureSize = value

        @property
        def FigureParam_LabelBoxWidth(self):
            return self.__FigureParam_LabelBoxWidth
        @FigureParam_LabelBoxWidth.setter
        def FigureParam_LabelBoxWidth(self,value):
            if not isinstance(value,int):
                raise ValueError('FigureParam_LabelBoxWidth must be be a tuple')
            if value <= 0:
                raise ValueError('ProbParam_RandomSeed must be bigger than 0')
            self.__FigureParam_LabelBoxWidth = value

        @property
        def Threshold_MaxOutDegree(self):
            return self.__Threshold_MaxOutDegree
        @Threshold_MaxOutDegree.setter
        def Threshold_MaxOutDegree(self,value):
            if not isinstance(value,int):
                raise ValueError('Threshold_MaxOutDegree must be be an integer')
            if value <= 0:
                raise ValueError('Threshold_MaxOutDegree must be bigger than 0')
            self.__Threshold_MaxOutDegree = value

        @property
        def Threshold_MinCellNumofStates(self):
            return self.__Threshold_MinCellNumofStates
        @Threshold_MinCellNumofStates.setter
        def Threshold_MinCellNumofStates(self,value):
            if not isinstance(value,int):
                raise ValueError('Threshold_MinCellNumofStates must be an integer')
            if value < 0:
                raise ValueError('Threshold_MinCellNumofStates must be bigger than 0')
            self.__Threshold_MinCellNumofStates = value

        @property
        def Threshold_MinProbability(self):
            return self.__Threshold_MinProbability
        @Threshold_MinProbability.setter
        def Threshold_MinProbability(self,value):
            if not isinstance(value,(int,float)) and value != "OTSU":
                raise ValueError('Threshold_MinProbability must be numeric')
            if value != "OTSU":
                if value < 0 or value > 1 :
                    raise ValueError('Threshold_MinProbability must be between 0.0 and 1.0') 
            self.__Threshold_MinProbability = value
        
        @property
        def KNNParam_metric(self):
            return self.__KNNParam_metric
        @KNNParam_metric.setter
        def KNNParam_metric(self,value):
            if not isinstance(value,str):
                raise ValueError('KNNParam_metric must be "euclidean" or "correlation"')
            self.__KNNParam_metric = value

    params=None
    timepoint_scdata_dict={}
    between_knn_graph=None
    between_cluster_graph=None
    filtered_cluster_node=None
    between_G=None
    __timepoint_scdata_num=1

    def __init__(self): 
        super(CStreetData, self).__init__()
        self.params=self.params_object()
    def __str__(self):
        s=""
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            s+=f"timepoint:{timepoint}\n{adata}\n\n"
        return s
    def __repr__(self):
        s=""
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            s+=f"timepoint:{timepoint}\n{adata}\n\n"
        return s
    # 定义装饰器
    def function_timer(function):
        @wraps(function)
        def function_timer(*args, **kwargs):
            print(f'\n[Function: {function.__name__} start...]\n')
            t0 = time.time()
            result = function(*args, **kwargs)
            t1 = time.time()
            print(f'\n[Function: {function.__name__} finished, spent time: {(t1 - t0):.2f}s]\n')
            print('='*80)
            return result
        return function_timer

    def except_output(function):
        @wraps(function)
        def execept_raise(*args, **kwargs):
            try:
                return function(*args, **kwargs)
            except Exception as e:
                sign = '=' * 80 + '\n'
                print(f'{sign}>>>Error Time:\t{datetime.now()}\n>>>Error Func:\t{function.__name__}\n>>>Error Info:\t{e}')
                print(f'{sign}{traceback.format_exc()}{sign}')
                raise e
        return execept_raise

    def params_filter(params_list):
        def params_filter(function):
            @wraps(function)
            def wrapper(*args, **kwargs):
                kwargs_filter={}
                for k, v in kwargs.items():
                    if k not in params_list:
                        print(f"Parameter '{k}' is invalid ,and it will be ignored.\n")
                    else:
                        kwargs_filter[k]=v
                return function(*args, **kwargs_filter)
            return wrapper
        return params_filter
    
    # 定义类方法
    @except_output
    @params_filter(["timepoint_scdata","timepoint_scdata_cluster"])
    def add_new_timepoint_scdata(self,timepoint_scdata,timepoint_scdata_cluster=None):
        data=pd.DataFrame(timepoint_scdata)
        data=data.fillna(0)
        ind=(data.std()>0.01).to_list()
        data=data.loc[:,ind]
        self.timepoint_scdata_dict[self.__timepoint_scdata_num]=ad.AnnData(data)
        if timepoint_scdata_cluster != None:
            self.timepoint_scdata_dict[self.__timepoint_scdata_num].obs["scdata_cluster"]=[f"timepoint{self.__timepoint_scdata_num}_{c}" for c in list(timepoint_scdata_cluster)]

            self.timepoint_scdata_dict[self.__timepoint_scdata_num].uns["cluster_flag"]=True
        else:
            self.timepoint_scdata_dict[self.__timepoint_scdata_num].obs["scdata_cluster"]=[0]*self.timepoint_scdata_dict[self.__timepoint_scdata_num].n_obs
            self.timepoint_scdata_dict[self.__timepoint_scdata_num].uns["cluster_flag"]=False
        self.__timepoint_scdata_num+=1
    
    def __create_folder(self):
        Output_Dir = self.params.Output_Dir
        Output_Name = self.params.Output_Name
        if not os.path.exists(Output_Dir):
            raise ValueError(f'{Output_Dir} : No such directory')
        elif not os.path.exists(Output_Dir+Output_Name):
            os.makedirs(Output_Dir+Output_Name)
        else:
            print(f"The result folder {Output_Dir+Output_Name} exists! CStreet overwrite it. To avoid the overwriting, try the -o parameter.")

        if not os.path.exists(Output_Dir+Output_Name+"/"+"SupplementaryFigures"):
            os.makedirs(Output_Dir+Output_Name+"/"+"SupplementaryFigures")
        if not os.path.exists(Output_Dir+Output_Name+"/"+"SupplementaryResults"):
            os.makedirs(Output_Dir+Output_Name+"/"+"SupplementaryResults")


    @except_output
    @function_timer
    @params_filter(['CellClusterParam_PCAn','CellClusterParam_k','CellClusterParam_Resolution','Switch_Normalize','Switch_LogTransform'])
    def cell_clusters(self,**kwargs):

        self.__create_folder()

        CellClusterParam_PCAn = self.params.CellClusterParam_PCAn = kwargs.setdefault('CellClusterParam_PCAn', self.params.CellClusterParam_PCAn)
        CellClusterParam_k = self.params.CellClusterParam_k = kwargs.setdefault('CellClusterParam_k', self.params.CellClusterParam_k)
        CellClusterParam_Resolution = self.params.CellClusterParam_Resolution=kwargs.setdefault("CellClusterParam_Resolution",self.params.CellClusterParam_Resolution)
        normalize_flag = self.params.Switch_Normalize=kwargs.setdefault("Switch_Normalize",self.params.Switch_Normalize)
        log_flag = self.params.Switch_LogTransform=kwargs.setdefault("Switch_LogTransform",self.params.Switch_LogTransform)
        
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            print(f"timepoint:{timepoint}")
            if adata.uns['cluster_flag'] :
                print(f"clusters have been given")
            else:
                adata_copy=adata.copy()
                adata_copy.obs_names_make_unique()
                adata_copy.var_names_make_unique()
                # MT pct
                mito_genes = adata_copy.var_names.str.startswith('MT-')
                adata_copy.obs['percent_mito'] = np.sum(adata_copy[:, mito_genes].X, axis=1) / np.sum(adata_copy.X, axis=1)
                adata_copy.obs['n_counts'] = adata_copy.X.sum(axis=1)
                # normalize log-transform
                if normalize_flag:
                    sc.pp.normalize_per_cell(adata_copy, counts_per_cell_after=1e4)
                if log_flag:
                    sc.pp.log1p(adata_copy)
                # high variable genes
                sc.pp.highly_variable_genes(adata_copy, min_mean=0.0125, max_mean=3, min_disp=0.5)
                adata_high = adata_copy[:, adata_copy.var['highly_variable']]
                # linear regression
                sc.pp.regress_out(adata_high, ['n_counts', 'percent_mito'])
                sc.pp.scale(adata_high, max_value=10)
                # pca
                sc.tl.pca(adata_high, n_comps=CellClusterParam_PCAn, svd_solver='arpack')
                # knn
                sc.pp.neighbors(adata_high, n_neighbors=CellClusterParam_k, n_pcs=CellClusterParam_PCAn)
                sc.tl.louvain(adata_high, resolution=CellClusterParam_Resolution)
                adata_high.obs["louvain"]=[str(int(i)+1) for i in adata_high.obs["louvain"]]
                sc.tl.umap(adata_high)
                umap_cord=pd.DataFrame(adata_high.obsm["X_umap"])
                umap_cord.index=adata_high.obs.index
                umap_cord["louvain"]=adata_high.obs["louvain"]
                umap_cord.columns=["x","y","louvain"]

                fig=sc.pl.umap(adata_high, color='louvain',return_fig=True)
                plt.show(block=False)
                plt.pause(1.0)
                plt.close()
                ouput_path=self.params.Output_Dir+self.params.Output_Name+"/SupplementaryFigures/"
                output_name=self.params.Output_Name
                fig.savefig(ouput_path+f"{output_name}_tp{timepoint}_LouvainUMAPClustering.pdf", bbox_inches='tight')
                umap_cord.to_csv(ouput_path+f"{output_name}_tp{timepoint}_LouvainUMAPClusteringCoordinates.txt",sep="\t")
                
                adata.obs["scdata_cluster"]=[f"timepoint{timepoint}_cluster{int(c)}" for c in adata_high.obs["louvain"].tolist()]
                adata.uns["cluster_flag"]=True
            # 按字符串排序
            node_cluster=adata.obs["scdata_cluster"]
            #node_cluster.index=adata.obs["cell_id"]
            cluster_set=node_cluster.unique().tolist()
            cluster_set.sort()
            adata.uns["cluster_set"]=cluster_set
            adata.uns["cluster_counts"]=node_cluster.value_counts()

    @except_output
    @function_timer
    @params_filter(["Threshold_MitoPercent"])
    def filter_dead_cell(self,**kwargs):
        Threshold_MitoPercent=self.params.Threshold_MitoPercent=kwargs.setdefault("Threshold_MitoPercent",self.params.Threshold_MitoPercent)

        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
            adata.obs['percent_mito'] = np.sum(adata[:, adata.var['mt']].X, axis=1) / np.sum(adata.X, axis=1)
            adata.obs['n_counts'] = adata.X.sum(axis=1)
            raw_cell_num=adata.n_obs
            adata=adata[adata.obs['percent_mito']<Threshold_MitoPercent,:]
            filter_cell_num=adata.n_obs
            self.timepoint_scdata_dict[timepoint]=adata
            print(f'timepoint:{timepoint}')
            print(f'filtered out {raw_cell_num-filter_cell_num} cells that are detected in more than {Threshold_MitoPercent} mito percent')
            print()

    @except_output
    @function_timer
    @params_filter(["Threshold_LowCellNum"])
    def filter_lowcell_gene(self,**kwargs):
        Threshold_LowCellNum=self.params.Threshold_LowCellNum=kwargs.setdefault("Threshold_LowCellNum",self.params.Threshold_LowCellNum)
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            raw_gene_num=adata.n_vars
            sc.pp.filter_genes(adata, min_cells=Threshold_LowCellNum)
            filter_gene_num=adata.n_vars
            print(f'timepoint:{timepoint}')
            print(f'filtered out {raw_gene_num-filter_gene_num} genes that are detected in less than {Threshold_LowCellNum} cells')
            print()
    
    @except_output
    @function_timer
    @params_filter(["Threshold_LowGeneNum"])
    def filter_lowgene_cells(self,**kwargs):
        Threshold_LowGeneNum=self.params.Threshold_LowGeneNum=kwargs.setdefault("Threshold_LowGeneNum",self.params.Threshold_LowGeneNum)

        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            raw_cell_num=adata.n_obs
            sc.pp.filter_cells(adata, min_genes=Threshold_LowGeneNum)
            filter_cell_num=adata.n_obs
            print(f'timepoint:{timepoint}')
            print(f'filtered out {raw_cell_num-filter_cell_num} cells that are detected in less than {Threshold_LowGeneNum} genes')
            print()

    @except_output
    @function_timer
    @params_filter(["Threshold_NormalizeBase"])
    def normalize_data(self,**kwargs):
        Threshold_NormalizeBase=self.params.Threshold_NormalizeBase=kwargs.setdefault("Threshold_NormalizeBase",self.params.Threshold_NormalizeBase)
        print(f'Normalize data to {Threshold_NormalizeBase} count ...')
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            sc.pp.normalize_total(adata,target_sum=Threshold_NormalizeBase)

    @except_output
    @function_timer
    def log_transform(self):
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            sc.pp.log1p(adata)
    
    @except_output
    @function_timer
    @params_filter(["WithinTimePointParam_PCAn","WithinTimePointParam_k"])
    def get_knn_within(self,**kwargs):
        pca_n=self.params.WithinTimePointParam_PCAn=kwargs.setdefault("WithinTimePointParam_PCAn",self.params.WithinTimePointParam_PCAn)
        k=self.params.WithinTimePointParam_k=kwargs.setdefault("WithinTimePointParam_k",self.params.WithinTimePointParam_k)
        KNNParam_metric=self.params.KNNParam_metric=kwargs.setdefault("KNNParam_metric",self.params.KNNParam_metric)
        if KNNParam_metric != "euclidean" :
            KNNParam_metric = correlation
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            print(f"timepoint:{timepoint}")
            df=pd.DataFrame(adata.X)
            df.index=adata.obs_names
            df.columns=adata.var_names
            # zscore
            df_zscore = df.apply(zscore)
            
            # PCA
            if not pca_n:
                n_components = min(min(df.shape), 100)
            else:
                n_components = pca_n
            pca = PCA(n_components=n_components, svd_solver='full')
            df_zscore_pca = pca.fit_transform(df_zscore)
            sample_name = np.array(['timepoint%s_cell%s'%(timepoint, j) for j in range(df.shape[0])])

            # KNN graph
            nbrs = NearestNeighbors(n_neighbors=k+1, metric=KNNParam_metric, n_jobs=-2)
            nbrs.fit(df_zscore_pca)
            dists, ind = nbrs.kneighbors(df_zscore_pca) # n*k_init dist & index
            adj = nbrs.kneighbors_graph(df_zscore_pca, mode='distance') # n*n dist
            adata.obs["cell_id"]=sample_name
            adata.obsm["within_dists"]=dists[:, 1:]
            adata.obsm["within_ind"]=sample_name[ind[:, 1:]]
            adata.obsm['within_distance']=adj
            print()

    @except_output
    @function_timer
    @params_filter(["BetweenTimePointParam_PCAn","BetweenTimePointParam_k"])
    def get_knn_between(self,**kwargs):
        pca_n=self.params.BetweenTimePointParam_PCAn=kwargs.setdefault("BetweenTimePointParam_PCAn",self.params.BetweenTimePointParam_PCAn)
        k=self.params.BetweenTimePointParam_k=kwargs.setdefault("BetweenTimePointParam_k",self.params.BetweenTimePointParam_k)
        KNNParam_metric=self.params.KNNParam_metric=kwargs.setdefault("KNNParam_metric",self.params.KNNParam_metric)
        
        if KNNParam_metric != "euclidean" :
            KNNParam_metric = correlation
        for timepoint in list(self.timepoint_scdata_dict.keys())[:-1]:
            print(f"timepoint between {timepoint} and {timepoint+1} ")


            adata_start=self.timepoint_scdata_dict[timepoint]
            df_start=pd.DataFrame(adata_start.X)
            df_start.index=adata_start.obs["cell_id"]
            df_start.columns=adata_start.var_names
            
            adata_end=self.timepoint_scdata_dict[timepoint+1]
            df_end=pd.DataFrame(adata_end.X)
            df_end.index=adata_end.obs["cell_id"]
            df_end.columns=adata_end.var_names

            common_genes = df_end.columns.intersection(df_start.columns)
            df_start=df_start.loc[:,common_genes]
            df_end=df_end.loc[:,common_genes]
            # zscore
            df_start = df_start.apply(zscore)
            df_end = df_end.apply(zscore)
            
            # PCA for data t+1, and projection data t to t+1(need to define PC number)
            if not pca_n:
                n_components = min(min(df_end.shape), 100)
            else:
                n_components = pca_n
            pca = PCA(n_components=n_components, svd_solver='full')
            df_end_pca = pca.fit_transform(df_end)
            df_start_pca = pca.transform(df_start)

            # KNN distance (start)
            nbrs = NearestNeighbors(n_neighbors=k+1, metric=KNNParam_metric, n_jobs=-2)
            nbrs.fit(df_end_pca)
            dists, ind = nbrs.kneighbors(df_start_pca) # n*k_init dist & index
            adj = nbrs.kneighbors_graph(df_start_pca, mode='distance') # n*n dist
            adata_start.obsm["between_later_dists"]=dists[:, 1:]
            adata_start.obsm["between_later_ind"]=adata_end.obs["cell_id"][ind[:, 1:]]
            adata_start.obsm['between_later_distance']=adj
            # KNN distance (end)
            nbrs = NearestNeighbors(n_neighbors=k+1, metric=KNNParam_metric, n_jobs=-2)
            nbrs.fit(df_start_pca)
            dists, ind = nbrs.kneighbors(df_end_pca) # n*k_init dist & index
            adj = nbrs.kneighbors_graph(df_end_pca, mode='distance') # n*n dist
            adata_end.obsm["between_previous_dists"]=dists[:, 1:]
            adata_end.obsm["between_previous_ind"]=adata_start.obs["cell_id"][ind[:, 1:]]
            adata_end.obsm['between_previous_distance']=adj

    @except_output
    @function_timer
    def get_knn_graph(self):
        #within graph
        for (timepoint,adata) in self.timepoint_scdata_dict.items():

            print(f"timepoint:{timepoint}")
            node1 = np.repeat(adata.obs["cell_id"], adata.obsm["within_dists"].shape[1])
            node2 = adata.obsm["within_ind"].flatten()
            Distance = adata.obsm["within_dists"].flatten()
            Distance_normalize = (adata.obsm["within_dists"]/adata.obsm["within_dists"][:,0][:,None]).flatten()
            Distance_zscore = zscore(Distance)
            
            df_graph=pd.DataFrame({
                'Node1':node1,
                'Node2':node2,
                'Distance':Distance,
                'Distance_normalize':Distance_normalize,
                'Distance_zscore':Distance_zscore
                })
            adata.uns['within_knn_graph']=df_graph


        #between graph
        for timepoint in list(self.timepoint_scdata_dict.keys())[:-1]:
            print(f"timepoint between {timepoint} and {timepoint+1} ")
            adata_start=self.timepoint_scdata_dict[timepoint]
            node1 = np.repeat(adata_start.obs["cell_id"], adata_start.obsm["between_later_dists"].shape[1])
            node2 = adata_start.obsm["between_later_ind"].flatten()
            Distance = adata_start.obsm["between_later_dists"].flatten()
            Distance_normalize = (adata_start.obsm["between_later_dists"]/adata_start.obsm["between_later_dists"][:,0][:,None]).flatten()
            Distance_zscore = zscore(Distance)

            df_graph=pd.DataFrame({
                'Node1':node1,
                'Node2':node2,
                'Distance':Distance,
                'Distance_normalize':Distance_normalize,
                'Distance_zscore':Distance_zscore
                })
            adata_start.uns['between_later_knn_graph']=df_graph.reset_index(drop=True)

            adata_end=self.timepoint_scdata_dict[timepoint+1]
            node1 = np.repeat(adata_end.obs["cell_id"], adata_end.obsm["between_previous_dists"].shape[1])
            node2 = adata_end.obsm["between_previous_ind"].flatten()
            Distance = adata_end.obsm["between_previous_dists"].flatten()
            Distance_normalize = (adata_end.obsm["between_previous_dists"]/adata_end.obsm["between_previous_dists"][:,0][:,None]).flatten()
            Distance_zscore = zscore(Distance)
            df_graph=pd.DataFrame({
                'Node1':node2,
                'Node2':node1,
                'Distance':Distance,
                'Distance_normalize':Distance_normalize,
                'Distance_zscore':Distance_zscore                
                })
            adata_end.uns['between_previous_knn_graph']=df_graph

    @except_output
    @function_timer
    def filter_knn_graph(self):
        knn_graph=pd.DataFrame()
        #within graph
        def Maximum_Score(df):
            Q1=df.quantile(0.25)
            Q3=df.quantile(0.75)
            return Q3+(Q3-Q1)*1.5

        for (timepoint,adata) in self.timepoint_scdata_dict.items():

            print(f"timepoint:{timepoint}")
            df_graph=adata.uns['within_knn_graph']
            within_max_normalize_keep_id=df_graph.loc[:,'Distance_normalize']<=Maximum_Score(df_graph.loc[:,'Distance_normalize'])
            within_max_zscore_keep_id=df_graph.loc[:,'Distance_zscore']<=Maximum_Score(df_graph.loc[:,'Distance_zscore'])
            keep_id=within_max_normalize_keep_id & within_max_zscore_keep_id
            df_graph=df_graph.loc[keep_id,:].reset_index(drop=True)

            # df_graph_edge=df_graph.loc[:,['Node1', 'Node2']].apply(lambda x:x.sort_values().str.cat(sep="->"), axis=1)
            # df_graph.loc[:,['Node1', 'Node2']]=df_graph_edge.str.split("->", expand=True)
            
            
            # duplicate_id=df_graph.loc[:,['Node1', 'Node2']].apply(lambda x:x.sort_values().str.cat(), axis=1).duplicated()
            
            # if False:
            #     df_graph=df_graph.loc[duplicate_id,:].reset_index(drop=True)
            # else:
            #     df_graph=df_graph.loc[~duplicate_id,:].reset_index(drop=True)

            adata.uns['within_knn_graph']=df_graph

        #between graph
        # between_max_normalize=3
        # between_max_zscore=0

        for timepoint in list(self.timepoint_scdata_dict.keys())[:-1]:
            print(f"timepoint between {timepoint} and {timepoint+1} ")
            adata_start=self.timepoint_scdata_dict[timepoint]
            df_graph=adata_start.uns['between_later_knn_graph']
            between_max_normalize_keep_id=df_graph.loc[:,'Distance_normalize']<=Maximum_Score(df_graph.loc[:,'Distance_normalize'])
            between_max_zscore_keep_id=df_graph.loc[:,'Distance_zscore']<=Maximum_Score(df_graph.loc[:,'Distance_zscore'])
            keep_id=between_max_normalize_keep_id & between_max_zscore_keep_id
            df_graph=df_graph.loc[keep_id,:].reset_index(drop=True)


            adata_end=self.timepoint_scdata_dict[timepoint+1]
            df_graph2=adata_end.uns['between_previous_knn_graph']
            between_max_normalize_keep_id=df_graph2.loc[:,'Distance_normalize']<=Maximum_Score(df_graph2.loc[:,'Distance_normalize'])
            between_max_zscore_keep_id=df_graph2.loc[:,'Distance_zscore']<=Maximum_Score(df_graph2.loc[:,'Distance_zscore'])
            keep_id=between_max_normalize_keep_id & between_max_zscore_keep_id
            df_graph2=df_graph2.loc[keep_id,:].reset_index(drop=True)

            df_graph_new=pd.concat([df_graph2,df_graph])
            # if False :
            #     df_graph_new=pd.concat([df_graph2,df_graph])
            #     df_graph_new=df_graph_new.loc[df_graph_new.duplicated(['Node1','Node2']),:].reset_index(drop=True)
            # else:
            #     df_graph_new=pd.concat([df_graph2,df_graph])
            #     df_graph_new=df_graph_new.loc[~df_graph_new.duplicated(['Node1','Node2']),:].reset_index(drop=True)

            knn_graph=pd.concat([knn_graph,df_graph_new])
        self.knn_graph=knn_graph

    def __get_cgKNN(self,knn_graph,node_cluster,within=True,ProbParam_RandomSeed=6886,k=15,ProbParam_SamplingSize=10):
        knn_graph.loc[:,'Node1_cluster'] = knn_graph.loc[:,'Node1'].map(node_cluster["scdata_cluster"].to_dict())
        knn_graph.loc[:,'Node2_cluster'] = knn_graph.loc[:,'Node2'].map(node_cluster["scdata_cluster"].to_dict())
        
        #去除相同类型的细胞之间的连接
        if within :
            nonself_idx = knn_graph.apply(lambda x:x['Node1_cluster']!=x['Node2_cluster'], axis=1)
            knn_graph=knn_graph.loc[nonself_idx,:]

        if knn_graph.shape[0]==0:
            return knn_graph

        knn_graph["total_edge_num"]=1
        knn_graph_egde=knn_graph.groupby(['Node1_cluster','Node2_cluster'],as_index=False)['total_edge_num'].count()
        knn_graph_egde['Node1_cluster_cell_num']=knn_graph_egde['Node1_cluster'].map(node_cluster["scdata_cluster"].value_counts().to_dict())
        knn_graph_egde['Node2_cluster_cell_num']=knn_graph_egde['Node2_cluster'].map(node_cluster["scdata_cluster"].value_counts().to_dict())
        #knn_graph_egde['rate']=knn_graph_egde["total_edge_num"]/(knn_graph_egde['Node2_cluster_cell_num']*knn_graph_egde['Node1_cluster_cell_num'])

        def get_edge_score(df,node_cluster,ProbParam_RandomSeed,k,ProbParam_SamplingSize):
            Node1_cluster=df["Node1_cluster"]
            Node2_cluster=df["Node2_cluster"]
            
            Node1_df=node_cluster.loc[node_cluster["scdata_cluster"]==Node1_cluster,:]
            Node2_df=node_cluster.loc[node_cluster["scdata_cluster"]==Node2_cluster,:]
            

            def tmp_filter(df,Node1_cluster_sample_cells,Node2_cluster_sample_cells):
                if df["Node1"].isin(Node1_cluster_sample_cells) & df["Node2"].isin(Node2_cluster_sample_cells):
                    return 1
                else:
                    return 0

            results_list=[]
            random.seed(ProbParam_RandomSeed)
            for i in range(ProbParam_SamplingSize):

                if random.choice([True, False]):
                    Node1_cluster_sample_cells=list(Node1_df.sample(frac=0.5, random_state=ProbParam_RandomSeed+i, axis=0)["cell_id"])
                    Node2_cluster_sample_cells=list(Node2_df.sample(frac=1, random_state=ProbParam_RandomSeed+i, axis=0)["cell_id"])
                else:
                    Node1_cluster_sample_cells=list(Node1_df.sample(frac=1, random_state=ProbParam_RandomSeed+i, axis=0)["cell_id"])
                    Node2_cluster_sample_cells=list(Node2_df.sample(frac=0.5, random_state=ProbParam_RandomSeed+i, axis=0)["cell_id"])
                # !!!spend too many time!!!
                tmp_edge_num=sum(knn_graph["Node1"].isin(Node1_cluster_sample_cells) & knn_graph["Node2"].isin(Node2_cluster_sample_cells))
                # !!!!!!!!!!!!!!!!!!!!!!!!!
                a=len(Node1_cluster_sample_cells)
                b=len(Node2_cluster_sample_cells)
                tmp_edge_num_max=(a+b)*k
                results_list.append(tmp_edge_num/tmp_edge_num_max)
            
            mean=np.mean(results_list)
            sd=np.std(results_list)
            n=ProbParam_SamplingSize
            dof = n - 1
            alpha = 1.0 - 0.95
            conf_interval = t.ppf(1-alpha/2., dof) * sd/np.sqrt(n)
            conf_interval = f"({mean-conf_interval:.5f},{mean+conf_interval:.5f})"
            return [mean,conf_interval]+results_list
        knn_graph_egde['probability_rep']=knn_graph_egde.apply(get_edge_score,args=(node_cluster,ProbParam_RandomSeed,k,ProbParam_SamplingSize),axis=1)
        knn_graph_egde[["probability"]+["95_conf_interval"]+[f"probability_rep{i+1}" for i in range(ProbParam_SamplingSize)]]=knn_graph_egde['probability_rep'].apply(pd.Series)
        knn_graph_egde=knn_graph_egde.drop(columns='probability_rep')

        return knn_graph_egde



    @except_output
    @function_timer
    @params_filter(["ProbParam_SamplingSize","ProbParam_RandomSeed"])
    def get_state_trajectory(self,**kwargs):
        ProbParam_SamplingSize=self.params.ProbParam_SamplingSize=kwargs.setdefault("ProbParam_SamplingSize",self.params.ProbParam_SamplingSize)
        ProbParam_RandomSeed=self.params.ProbParam_RandomSeed=kwargs.setdefault("ProbParam_RandomSeed",self.params.ProbParam_RandomSeed)        
        within_k=self.params.WithinTimePointParam_k
        between_k=self.params.BetweenTimePointParam_k
        #within graph
        all_node_cluster=pd.DataFrame()
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            print(f"timepoint:{timepoint}")
            
            knn_graph=adata.uns["within_knn_graph"]
            node_cluster=adata.obs[["cell_id","scdata_cluster"]]
            node_cluster.index=adata.obs["cell_id"]
            adata.uns["cluster_graph"]=self.__get_cgKNN(knn_graph,node_cluster,within=True,ProbParam_RandomSeed=ProbParam_RandomSeed,k=within_k,ProbParam_SamplingSize=ProbParam_SamplingSize)
            
            all_node_cluster=pd.concat([all_node_cluster,node_cluster])

        #between graph
        knn_graph_egde=self.__get_cgKNN(self.knn_graph,all_node_cluster,within=False,ProbParam_RandomSeed=ProbParam_RandomSeed,k=between_k,ProbParam_SamplingSize=ProbParam_SamplingSize)

        self.cluster_graph=knn_graph_egde

    def __get_nxG(self,knn_graph_egde, allnodes,min_score):
        G=nx.DiGraph()
        for node in allnodes:
            G.add_node(node,timepoint=str(node).split("_")[0])
        for index,row in knn_graph_egde.iterrows():
            if (row["Node1_cluster"] in allnodes) & (row["Node2_cluster"] in allnodes) & (row["probability"]>=min_score) :
                G.add_edge(row["Node1_cluster"],row["Node2_cluster"],probability=row["probability"],start_timepoint=G.nodes[row["Node1_cluster"]]["timepoint"])
        return G

    def __if_ZeroDivisionError(exception):
        return isinstance(exception, ZeroDivisionError)

    @retry(retry_on_exception=__if_ZeroDivisionError)
    def __force_directed_layout(self,G,timepoint,verbose=True, iterations=2000):
        """" Function to compute force directed layout from the G
        :param G: networkx graph object converted from sparse matrix representing affinities between cells
        :param cell_names: pandas Series object with cell names
        :param verbose: Verbosity for force directed layout computation
        :param iterations: Number of iterations used by ForceAtlas 
        :return: Pandas data frame representing the force directed layout
        """

        forceatlas2 = ForceAtlas2(
            # Behavior alternatives
            outboundAttractionDistribution=False,  
            linLogMode=False,  
            adjustSizes=False,  
            edgeWeightInfluence=1.0,
            # Performance
            jitterTolerance=1.0,  
            barnesHutOptimize=True,
            barnesHutTheta=1.2,
            multiThreaded=False,  
            # Tuning
            scalingRatio=2.0,
            strongGravityMode=False,
            gravity=1.0,
            # Log
            verbose=verbose)

        ## use affinity construct KNN graph
        if len(list(G.nodes))==1:
            positions={f"{list(G.nodes)[0]}":(5,5)}
        else:
            positions = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=iterations)
            x_max=pd.DataFrame(np.array([list(positions[i]) for i in positions.keys()]))[0].max()
            x_min=pd.DataFrame(np.array([list(positions[i]) for i in positions.keys()]))[0].min()
            
            y_max=pd.DataFrame(np.array([list(positions[i]) for i in positions.keys()]))[1].max()
            y_min=pd.DataFrame(np.array([list(positions[i]) for i in positions.keys()]))[1].min()
            for node,(x,y) in positions.items():
                positions[node]=((x-x_min)/(x_max-x_min)*10,(y-y_min)/(y_max-y_min)*10)

        ## plot KNN graph
        f = plt.figure(figsize=(20,15))
        ax=f.add_subplot(111)
        ax.set(xlim=[-1, 11],ylim=[-1, 11],title='ForceAtlas')
        nx.draw_networkx_nodes(G, positions, ax=ax, node_size=800, node_color="blue", alpha=0.5)
        nx.draw_networkx_edges(G, positions, ax=ax, width=10,edge_color="green", alpha=0.5)
        nx.draw_networkx_labels(G, positions, ax=ax, font_size=20)
        plt.tight_layout()
        plt.show(block=False)
        plt.pause(0.05)
        plt.close()

        ouput_path=self.params.Output_Dir+self.params.Output_Name+"/SupplementaryFigures/"
        output_name=self.params.Output_Name
        f.savefig(ouput_path+f"{output_name}_tp{timepoint}_ForceDirectedLayout.pdf", bbox_inches='tight')

        return positions

    @except_output
    @function_timer
    @params_filter(["Threshold_MaxOutDegree","Threshold_MinCellNumofStates","Threshold_MinProbability"])
    def get_knn_nxG(self,**kwargs):
        topN=self.params.Threshold_MaxOutDegree=kwargs.setdefault("Threshold_MaxOutDegree",self.params.Threshold_MaxOutDegree)
        #min_score=self.params.min_score=kwargs.setdefault("min_score",self.params.min_score)
        Threshold_MinProbability=self.params.Threshold_MinProbability=kwargs.setdefault("Threshold_MinProbability",self.params.Threshold_MinProbability)
        Threshold_MinCellNumofStates=self.params.Threshold_MinCellNumofStates=kwargs.setdefault("Threshold_MinCellNumofStates",self.params.Threshold_MinCellNumofStates)

        def otsu(p_list):
            p_list_unique=np.unique(p_list)
            threshold_p = 0
            max_g = 0
            for p in p_list_unique:
                n0 = p_list[p_list>=p]
                n1 = p_list[p_list<p]
                w0 = len(n0)/len(p_list)
                w1 = len(n1)/len(p_list)
                u0 = np.mean(n0) if len(n0) > 0 else 0.
                u1 = np.mean(n1) if len(n0) > 0 else 0.
                
                g = w0 * w1 * (u0 - u1) ** 2
                if g > max_g:
                    max_g = g
                    threshold_p = p
            return threshold_p

        filtered_node_cluster=np.array([])
        node_cluster_num=[]
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            print(f"timepoint:{timepoint}")
            knn_graph_egde=adata.uns['cluster_graph'].loc[:,["Node1_cluster","Node2_cluster","probability"]]
            allnodes=adata.uns["cluster_set"]
            node_cluster_num.append(len(allnodes))
            
            p_list=knn_graph_egde["probability"]
            zero_n=(len(allnodes)**2-len(allnodes))-len(p_list)
            p_list=np.append(p_list,[0.0]*zero_n)
            if Threshold_MinProbability =="OTSU":
                min_score=otsu(p_list)
            else:
                min_score=Threshold_MinProbability
            adata.uns["min_score"]=min_score

            allnodes_counts=adata.uns["cluster_counts"]
            filtered_nodes=[node for node in allnodes if allnodes_counts[node]>Threshold_MinCellNumofStates]
            if len(filtered_nodes)==0 :
                raise ValueError(f"No cell state is retained in timepoint {timepoint} because of a high Threshold_MinCellNumofStates value")
            filtered_node_cluster=np.hstack([filtered_node_cluster,filtered_nodes])
            
            within_G = self.__get_nxG(knn_graph_egde, filtered_nodes,min_score).to_undirected()
            adata.uns["within_G"]=within_G
            
            fa_cord_all = self.__force_directed_layout(G=within_G,timepoint=timepoint)
            adata.uns["fa_cord"]=fa_cord_all

        knn_graph_egde=self.cluster_graph.sort_values(by='probability', ascending=False).groupby(["Node1_cluster"],as_index=False).head(topN)
        self.filtered_cluster_node=filtered_node_cluster
        p_list=knn_graph_egde["probability"]

        edge_n=0
        for i in range(len(node_cluster_num)-1):
            edge_n+=node_cluster_num[i]*node_cluster_num[i+1]
        zero_n=edge_n-len(p_list)
        p_list=np.append(p_list,[0.0]*zero_n)
        if Threshold_MinProbability =="OTSU":
            min_score=otsu(p_list)
        else:
            min_score=Threshold_MinProbability
        self.min_score=min_score
        between_G=self.__get_nxG(knn_graph_egde, filtered_node_cluster,min_score)
        self.between_G=between_G
    
    @except_output
    @function_timer
    @params_filter(["FigureParam_LabelBoxWidth","FigureParam_FigureSize"])
    def draw_nxG(self,**kwargs):
        FigureParam_FigureSize=self.params.FigureParam_FigureSize=kwargs.setdefault("FigureParam_FigureSize",self.params.FigureParam_FigureSize)
        FigureParam_LabelBoxWidth=self.params.FigureParam_LabelBoxWidth=kwargs.setdefault("FigureParam_LabelBoxWidth",self.params.FigureParam_LabelBoxWidth)
        def draw_networkx_edges(
            G,
            pos,
            pos2,
            edgelist=None,
            width=1.0,
            edge_color="k",
            style="solid",
            alpha=None,
            arrowstyle="-|>",
            arrowsize=10,
            edge_cmap=None,
            edge_vmin=None,
            edge_vmax=None,
            ax=None,
            arrows=True,
            label=None,
            node_size=300,
            nodelist=None,
            node_shape="o",
            connectionstyle=None,
            min_source_margin=0,
            min_target_margin=0,
        ):
            try:
                import matplotlib.pyplot as plt
                from matplotlib.colors import colorConverter, Colormap, Normalize
                from matplotlib.collections import LineCollection
                from matplotlib.patches import FancyArrowPatch
                import numpy as np
            except ImportError as e:
                raise ImportError("Matplotlib required for draw()") from e
            except RuntimeError:
                print("Matplotlib unable to open display")
                raise

            if ax is None:
                ax = plt.gca()

            if edgelist is None:
                edgelist = list(G.edges())

            if len(edgelist) == 0:  # no edges!
                if not G.is_directed() or not arrows:
                    return LineCollection(None)
                else:
                    return []

            if nodelist is None:
                nodelist = list(G.nodes())

            # FancyArrowPatch handles color=None different from LineCollection
            if edge_color is None:
                edge_color = "k"

            # set edge positions
            edge_pos = np.asarray([(pos[e[0]], pos2[e[1]]) for e in edgelist])

            # Check if edge_color is an array of floats and map to edge_cmap.
            # This is the only case handled differently from matplotlib
            if (
                np.iterable(edge_color)
                and (len(edge_color) == len(edge_pos))
                and np.alltrue([isinstance(c, (int,float)) for c in edge_color])
            ):
                if edge_cmap is not None:
                    assert isinstance(edge_cmap, Colormap)
                else:
                    edge_cmap = plt.get_cmap()
                if edge_vmin is None:
                    edge_vmin = min(edge_color)
                if edge_vmax is None:
                    edge_vmax = max(edge_color)
                color_normal = Normalize(vmin=edge_vmin, vmax=edge_vmax)
                edge_color = [edge_cmap(color_normal(e)) for e in edge_color]

            if not G.is_directed() or not arrows:
                edge_collection = LineCollection(
                    edge_pos,
                    colors=edge_color,
                    linewidths=width,
                    antialiaseds=(1,),
                    linestyle=style,
                    transOffset=ax.transData,
                    alpha=alpha,
                )

                edge_collection.set_cmap(edge_cmap)
                edge_collection.set_clim(edge_vmin, edge_vmax)

                edge_collection.set_zorder(1)  # edges go behind nodes
                edge_collection.set_label(label)
                ax.add_collection(edge_collection)

                return edge_collection

            arrow_collection = None

            if G.is_directed() and arrows:
                # Note: Waiting for someone to implement arrow to intersection with
                # marker.  Meanwhile, this works well for polygons with more than 4
                # sides and circle.

                def to_marker_edge(marker_size, marker):
                    if marker in "s^>v<d":  # `large` markers need extra space
                        return np.sqrt(2 * marker_size) / 2
                    else:
                        return np.sqrt(marker_size) / 2

                # Draw arrows with `matplotlib.patches.FancyarrowPatch`
                arrow_collection = []
                mutation_scale = arrowsize  # scale factor of arrow head

                # FancyArrowPatch doesn't handle color strings
                arrow_colors = colorConverter.to_rgba_array(edge_color, alpha)
                for i, (src, dst) in enumerate(edge_pos):
                    x1, y1 = src
                    x2, y2 = dst
                    shrink_source = 0  # space from source to tail
                    shrink_target = 0  # space from  head to target
                    if np.iterable(node_size):  # many node sizes
                        source, target = edgelist[i][:2]
                        source_node_size = node_size[nodelist.index(source)]
                        target_node_size = node_size[nodelist.index(target)]
                        shrink_source = to_marker_edge(source_node_size, node_shape)
                        shrink_target = to_marker_edge(target_node_size, node_shape)
                    else:
                        shrink_source = shrink_target = to_marker_edge(node_size, node_shape)

                    if shrink_source < min_source_margin:
                        shrink_source = min_source_margin

                    if shrink_target < min_target_margin:
                        shrink_target = min_target_margin

                    if len(arrow_colors) == len(edge_pos):
                        arrow_color = arrow_colors[i]
                    elif len(arrow_colors) == 1:
                        arrow_color = arrow_colors[0]
                    else:  # Cycle through colors
                        arrow_color = arrow_colors[i % len(arrow_colors)]

                    if np.iterable(width):
                        if len(width) == len(edge_pos):
                            line_width = width[i]
                        else:
                            line_width = width[i % len(width)]
                    else:
                        line_width = width

                    arrow = FancyArrowPatch(
                        (x1, y1),
                        (x2, y2),
                        arrowstyle=arrowstyle,
                        shrinkA=shrink_source,
                        shrinkB=shrink_target,
                        mutation_scale=mutation_scale,
                        color=arrow_color,
                        linewidth=line_width,
                        connectionstyle=connectionstyle,
                        linestyle=style,
                        zorder=1,
                    )  # arrows go behind nodes

                    # There seems to be a bug in matplotlib to make collections of
                    # FancyArrowPatch instances. Until fixed, the patches are added
                    # individually to the axes instance.
                    arrow_collection.append(arrow)
                    ax.add_patch(arrow)

            # update view
            minx = np.amin(np.ravel(edge_pos[:, :, 0]))
            maxx = np.amax(np.ravel(edge_pos[:, :, 0]))
            miny = np.amin(np.ravel(edge_pos[:, :, 1]))
            maxy = np.amax(np.ravel(edge_pos[:, :, 1]))

            w = maxx - minx
            h = maxy - miny
            padx, pady = 0.05 * w, 0.05 * h
            corners = (minx - padx, miny - pady), (maxx + padx, maxy + pady)
            ax.update_datalim(corners)
            ax.autoscale_view()

            ax.tick_params(
                axis="both",
                which="both",
                bottom=False,
                left=False,
                labelbottom=False,
                labelleft=False,
            )

            return arrow_collection

        def draw_networkx_labels(
            G,
            pos,
            labels=None,
            font_size=12,
            font_color="k",
            font_family="sans-serif",
            font_weight="normal",
            alpha=None,
            bbox=None,
            horizontalalignment="center",
            verticalalignment="center",
            ax=None,
            zorder=3
        ):
            try:
                import matplotlib.pyplot as plt
            except ImportError as e:
                raise ImportError("Matplotlib required for draw()") from e
            except RuntimeError:
                print("Matplotlib unable to open display")
                raise

            if ax is None:
                ax = plt.gca()

            if labels is None:
                labels = {n: n for n in G.nodes()}

            text_items = {}  # there is no text collection so we'll fake one
            for n, label in labels.items():
                (x, y) = pos[n]
                if not isinstance(label, str):
                    label = str(label)  # this makes "1" and 1 labeled the same
                t = ax.text(
                    x,
                    y,
                    label,
                    size=font_size,
                    color=font_color,
                    family=font_family,
                    weight=font_weight,
                    alpha=alpha,
                    horizontalalignment=horizontalalignment,
                    verticalalignment=verticalalignment,
                    transform=ax.transData,
                    bbox=bbox,
                    clip_on=True,
                    zorder=zorder
                )
                text_items[n] = t

            ax.tick_params(
                axis="both",
                which="both",
                bottom=False,
                left=False,
                labelbottom=False,
                labelleft=False,
            )

            return text_items
        ZORDER=0.01
        TRANS=25 
        SQ=np.array([[-2.5,12.5],[12.5,12.5],[12.5,-2.5],[-2.5,-2.5],[-2.5,12.5]]).T
        T=np.array([[1,0],[1,1.5]])
        SQ=T@SQ
        CCCOLOR=["#CE0013","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD","#16557A"]
        number_of_scdata=len(self.timepoint_scdata_dict)
        fig = plt.figure(figsize=FigureParam_FigureSize)
        ax = fig.add_subplot(1,1,1)
        
        all_pos={}
        all_pos2={}
        all_pos3={}
        
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            ax.fill(SQ[0]+TRANS*timepoint,SQ[1],facecolor="#CCCDCF",lw=3,alpha=0.8,zorder=ZORDER*timepoint)
            #ax.plot(SQ[0]+TRANS*timepoint,SQ[1],'black',lw=3,zorder=ZORDER*timepoint)
            plt.text(5+TRANS*timepoint,-5, f'$ t_{timepoint} $',size=15, color="black",weight="bold",verticalalignment="center",horizontalalignment="center")
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)

            G=adata.uns["within_G"]
            pos=adata.uns["fa_cord"].copy()
            pos2={}
            pos3={}
            pos4={}
            pos5={}
            width_list=[]
            for (u, v, wt) in G.edges.data("probability"):
                width_list.append(wt)
            if len(width_list)!=0:
                width_list=[3*i/max(width_list) for i in width_list]

            pos_array=np.array(list(pos.values())).T
            pos_array=(T@pos_array)
            pos_df=pd.DataFrame(pos_array.T)
            pos_df.index=list(pos.keys())

            for node in pos.keys():
                x=pos_df.loc[node,0]
                y=pos_df.loc[node,1]
                pos[node]=(x+TRANS*timepoint,y)
                
            labels={}
            labels2={}
            for i,nodes in enumerate(list(G.nodes())):
                labels[nodes]=str(i+1)
                labels2[nodes]=(" "*2+"_".join(str(nodes).split("_")[1:])+" "*100)[:FigureParam_LabelBoxWidth+2]
                x_pos=TRANS*timepoint
                y_pos=np.linspace(-10,-60,(G.number_of_nodes()+2))[1:-1][i]
                pos2[nodes]=(x_pos,y_pos)
                pos3[nodes]=(x_pos-1.2,y_pos) #label-bbox
                pos4[nodes]=(x_pos+FigureParam_LabelBoxWidth*1.3,y_pos) #start pos
                pos5[nodes]=(x_pos-1.5,y_pos) #end pos
            
            all_pos={**all_pos,**pos}
            all_pos2={**all_pos2,**pos4}
            all_pos3={**all_pos3,**pos5}
            #pos
            nodes=nx.draw_networkx_nodes(G, pos, ax=ax, node_size=100, node_color=CCCOLOR[timepoint], alpha=1)
            nodes.set_zorder((ZORDER)*timepoint+0.0005)
            edges=nx.draw_networkx_edges(G, pos, ax=ax, width=width_list,edge_color="black", alpha=1)
            if G.number_of_edges() != 0:
                edges.set_zorder((ZORDER)*timepoint+0.00015)
            nx.draw_networkx_labels(G, pos,labels, ax=ax, font_size=8, font_color="white",font_weight="bold")
            
            bbox={
                'facecolor': '#EEEEEE', #填充色
               'edgecolor':'#EEEEEE',#外框色
               'alpha': 1, #框透明度
               'pad': 1,#本文与框周围距离 
              }
            #pos2
            nodes=nx.draw_networkx_nodes(G, pos2, ax=ax, node_size=100, node_color=CCCOLOR[timepoint], alpha=1)
            nodes.set_zorder(3)
            
            nx.draw_networkx_labels(G, pos2,labels, ax=ax, font_size=8, font_color="white",font_weight="bold")
            draw_networkx_labels(G, pos3,labels2,horizontalalignment="left", ax=ax, font_size=12,bbox=bbox,font_family="monospace",zorder=1)
        
        between_G=self.between_G

        width_list=[]
        for (u, v, wt) in between_G.edges.data("probability"):
            width_list.append(wt)
        if len(width_list)!=0:
            width_list=[3*i/max(width_list) for i in width_list]
        timepoint_list=[]
        for (u, v, tp) in between_G.edges.data('start_timepoint'):
            timepoint_list.append(tp)    
        
        between_edges=nx.draw_networkx_edges(between_G, all_pos, ax=ax,width=width_list,edge_color="#0060A8",arrowstyle="->", arrowsize=20,alpha=1,min_source_margin=0)
        for i in range(len(timepoint_list)):
            tp=timepoint_list[i].split("timepoint")[1]
            between_edges[i].set_zorder((ZORDER)*int(tp)+0.00025)
        
        between_edges2=draw_networkx_edges(between_G, all_pos2,all_pos3, ax=ax,width=width_list,edge_color="#0060A8",arrowstyle="->", arrowsize=20,alpha=1,min_source_margin=10,node_size=100)
        for i in range(len(timepoint_list)):
            tp=timepoint_list[i].split("timepoint")[1]
            between_edges2[i].set_zorder((ZORDER)*int(tp)+0.025)
        plt.tight_layout()
        plt.show(block=False)
        plt.pause(0.05)
        plt.close()
        ouput_path=self.params.Output_Dir+self.params.Output_Name+"/"
        output_name=self.params.Output_Name
        fig.savefig(ouput_path+f"{output_name}_CStreetTopology.pdf", bbox_inches='tight')
    
    @except_output
    @function_timer
    def output_results(self):
        cell_num=self.params.Threshold_MinCellNumofStates
        ouput_path=self.params.Output_Dir+self.params.Output_Name+"/SupplementaryResults/"
        output_name=self.params.Output_Name
        all_cluster_graph=self.cluster_graph
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            print(f"timepoint:{timepoint}")
            min_score=adata.uns["min_score"]
            with open(ouput_path+f"{output_name}_tp{timepoint}_CellStatesConnectionProbabilities.txt","w") as fw:
                s=f'''#Only connections and cell states with following conditions have been displayed in the result figure:
# Connection probability >= {min_score:.5f}
# Total cell number of a state >= {cell_num}
'''
                fw.write(s)
            adata.obs.to_csv(ouput_path+f"{output_name}_tp{timepoint}_FilteredCellInfo.txt",sep="\t")
            adata.var.to_csv(ouput_path+f"{output_name}_tp{timepoint}_FilteredGeneInfo.txt",sep="\t")
            cluster_graph_file=adata.uns["cluster_graph"].loc[:,["Node1_cluster","Node2_cluster","total_edge_num","Node1_cluster_cell_num","Node2_cluster_cell_num","probability","95_conf_interval"]]
            cluster_graph_file.to_csv(ouput_path+f"{output_name}_tp{timepoint}_CellStatesConnectionProbabilities.txt",mode='a',sep="\t",float_format='%.5f')
            all_cluster_graph=pd.concat([all_cluster_graph,adata.uns["cluster_graph"]])
        min_score=self.min_score
        with open(ouput_path+f"{output_name}_BetweenTimePoints_CellStatesConnectionProbabilities.txt","w") as fw:
            s=f'''#Only connections and cell states with following conditions have been displayed in the result figure:
# Connection probability >= {min_score:.5f}
# Total cell number of a state >= {cell_num}
'''
            fw.write(s)
        cluster_graph_file=self.cluster_graph.loc[:,["Node1_cluster","Node2_cluster","total_edge_num","Node1_cluster_cell_num","Node2_cluster_cell_num","probability","95_conf_interval"]]
        cluster_graph_file.to_csv(ouput_path+f"{output_name}_BetweenTimePoints_CellStatesConnectionProbabilities.txt",mode='a',sep="\t",float_format='%.5f')
        all_cluster_graph=all_cluster_graph.loc[:,["Node1_cluster","Node2_cluster","probability","95_conf_interval"]]
        all_cluster_graph.columns=["SourceNode","TargetNode","ConnectionProbabilities","95ConfidenceIntervals"]
        all_cluster_graph.to_csv(f"{self.params.Output_Dir}{self.params.Output_Name}/{output_name}_CellStatesConnCytoscape.txt",mode='w',sep="\t",float_format='%.5f',index=False)
    @function_timer
    def run_cstreet(self):

        if True :
            self.cell_clusters()
        
        if self.params.Switch_DeadCellFilter == True:
            self.filter_dead_cell()
        
        if self.params.Switch_LowCellNumGeneFilter == True:
            self.filter_lowcell_gene()
        
        if self.params.Switch_LowGeneCellsFilter == True:
            self.filter_lowgene_cells()
        
        if self.params.Switch_Normalize == True :
            self.normalize_data()

        if self.params.Switch_LogTransform == True :
            self.log_transform()
        
        if True :
            self.get_knn_within()

            self.get_knn_between()

            self.get_knn_graph()

            self.filter_knn_graph()
            
            self.get_state_trajectory()
            
            self.get_knn_nxG()
            
            self.draw_nxG()

            self.output_results()
