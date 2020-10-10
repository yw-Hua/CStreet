#!python3.6
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

warnings.filterwarnings("ignore")

#定义cstreet对象

class CStreetData(object):
    """docstring for CStreetData"""
    class params_object:
        """docstring for params_object"""
        __slots__=("__output_dir","__output_name","__cell_cluster_pca_n","__cell_cluster_knn_n","__cell_cluster_resolution",
            "__filter_dead_cell","__percent_mito_cutoff","__filter_lowcell_gene","__min_cells","__filter_lowgene_cells",
            "__min_genes","__normalize","__normalize_base","__log_transform","__highly_variable_genes",
            "__inner_graph_pca_n","__inner_graph_knn_n","__link_graph_pca_n","__link_graph_knn_n",
            "__max_outgoing","__min_score","__min_cell_number")
        
        def __init__(self):
            #Step0:basic params#
            self.__output_dir="./"
            self.__output_name="cstreet_project"
            #Step1:cell_cluster# 
            self.__cell_cluster_pca_n=10
            self.__cell_cluster_knn_n=15
            self.__cell_cluster_resolution=1

            #Step2:gene and cell filter#
            self.__filter_dead_cell=True
            self.__percent_mito_cutoff=0.2
            
            self.__filter_lowcell_gene=True
            self.__min_cells=3
            
            self.__filter_lowgene_cells=True
            self.__min_genes=200
            
            #Step3:normalize#
            self.__normalize=True
            self.__normalize_base=10000
            self.__log_transform=True

            #Step4:get HVG#
            self.__highly_variable_genes=False

            #Step5:get_graph#
            self.__inner_graph_pca_n=10
            self.__inner_graph_knn_n=15

            self.__link_graph_pca_n=10
            self.__link_graph_knn_n=15

            #Step6:plot graph#
            self.__max_outgoing=10
            self.__min_score=0.1
            self.__min_cell_number=50
        
        def __str__(self):
            s=""
            s+=f"\n#Step0:basic params# \n"
            s+=f"output_dir={self.__output_dir}\n"
            s+=f"output_name={self.__output_name}\n"            
            
            s+=f"\n#Step1:cell_cluster# \n"
            s+=f"cell_cluster_pca_n={self.__cell_cluster_pca_n}\n"
            s+=f"cell_cluster_knn_n={self.__cell_cluster_knn_n}\n"
            s+=f"cell_cluster_resolution={self.__cell_cluster_resolution}\n"

            s+=f"\n#Step2:gene and cell filter#\n"
            s+=f"filter_dead_cell={self.__filter_dead_cell}\n"
            s+=f"percent_mito_cutoff={self.__percent_mito_cutoff}\n"

            s+=f"filter_lowcell_gene={self.__filter_lowcell_gene}\n"
            s+=f"min_cells={self.__min_cells}\n"

            s+=f"filter_lowgene_cells={self.__filter_lowgene_cells}\n"
            s+=f"min_genes={self.__min_genes}\n"

            s+=f"\n#Step3:normalize#\n"
            s+=f"normalize={self.__normalize}\n"
            s+=f"normalize_base={self.__normalize_base}\n"
            s+=f"log_transform={self.__log_transform}\n"

            s+=f"\n#Step4:get HVG#\n"
            s+=f"highly_variable_genes={self.__highly_variable_genes}\n"

            s+=f"\n#Step5:get_graph#\n"
            s+=f"inner_graph_pca_n={self.__inner_graph_pca_n}\n"
            s+=f"inner_graph_knn_n={self.__inner_graph_knn_n}\n"

            s+=f"link_graph_pca_n={self.__link_graph_pca_n}\n"
            s+=f"link_graph_knn_n={self.__link_graph_knn_n}\n"

            s+=f"\n#Step6:plot graph#\n"
            s+=f"max_outgoing={self.__max_outgoing}\n"
            s+=f"min_score={self.__min_score}\n"
            s+=f"min_cell_number={self.__min_cell_number}\n"
            return s
        def __repr__(self):
            return self.__str__()
        @property
        def output_dir(self):
            return self.__output_dir
        @output_dir.setter
        def output_dir(self,value):
            if not isinstance(value,str):
                raise ValueError('output_dir must be a string')
            if value[-1] != "/":
                value+="/"
            self.__output_dir = value

        @property
        def output_name(self):
            return self.__output_name
        @output_name.setter
        def output_name(self,value):
            if not isinstance(value,str):
                raise ValueError('output_name must be a string')
            self.__output_name = value

        @property
        def cell_cluster_pca_n(self):
            return self.__cell_cluster_pca_n
        @cell_cluster_pca_n.setter
        def cell_cluster_pca_n(self,value):
            if not isinstance(value,int):
                raise ValueError('cell_cluster_pca_n must be an integer')
            if value <= 0:
                raise ValueError('cell_cluster_pca_n must be bigger than 0')
            self.__cell_cluster_pca_n = value

        @property
        def cell_cluster_knn_n(self):
            return self.__cell_cluster_knn_n
        @cell_cluster_knn_n.setter
        def cell_cluster_knn_n(self,value):
            if not isinstance(value,int):
                raise ValueError('cell_cluster_knn_n must be an integer')
            if value <= 0:
                raise ValueError('cell_cluster_knn_n must be bigger than 0')
            self.__cell_cluster_knn_n = value

        @property
        def cell_cluster_resolution(self):
            return self.__cell_cluster_resolution
        @cell_cluster_resolution.setter
        def cell_cluster_resolution(self,value):
            if not isinstance(value,(int,float)):
                raise ValueError('cell_cluster_resolution must be numeric')
            if value <= 0:
                raise ValueError('cell_cluster_resolution must be bigger than 0')        
            self.__cell_cluster_resolution = value

        @property
        def filter_dead_cell(self):
            return self.__filter_dead_cell
        @filter_dead_cell.setter
        def filter_dead_cell(self,value):
            if not isinstance(value,bool):
                raise ValueError('filter_dead_cell must be True or False')
            self.__filter_dead_cell = value

        @property
        def percent_mito_cutoff(self):
            return self.__percent_mito_cutoff
        @percent_mito_cutoff.setter
        def percent_mito_cutoff(self,value):
            if not isinstance(value,(int,float)):
                raise ValueError('percent_mito_cutoff must be numeric')
            if value <= 0 or value >= 1:
                raise ValueError('percent_mito_cutoff must be between 0.0 and 1.0') 
            self.__percent_mito_cutoff = value

        @property
        def filter_lowcell_gene(self):
            return self.__filter_lowcell_gene
        @filter_lowcell_gene.setter
        def filter_lowcell_gene(self,value):
            if not isinstance(value,bool):
                raise ValueError('filter_lowcell_gene must be True or False')
            self.__filter_lowcell_gene = value

        @property
        def min_cells(self):
            return self.__min_cells
        @min_cells.setter
        def min_cells(self,value):
            if not isinstance(value,int):
                raise ValueError('min_cells must be an integer')
            if value <= 0:
                raise ValueError('min_cells must be bigger than 0')
            self.__min_cells = value

        @property
        def filter_lowgene_cells(self):
            return self.__filter_lowgene_cells
        @filter_lowgene_cells.setter
        def filter_lowgene_cells(self,value):
            if not isinstance(value,bool):
                raise ValueError('filter_lowgene_cells must be True or False')
            self.__filter_lowgene_cells = value

        @property
        def min_genes(self):
            return self.__min_genes
        @min_genes.setter
        def min_genes(self,value):
            if not isinstance(value,int):
                raise ValueError('min_genes must be an integer')
            if value <= 0:
                raise ValueError('min_genes must be bigger than 0')
            self.__min_genes = value

        @property
        def normalize(self):
            return self.__normalize
        @normalize.setter
        def normalize(self,value):
            if not isinstance(value,bool):
                raise ValueError('normalize must be True or False')
            self.__normalize = value

        @property
        def normalize_base(self):
            return self.__normalize_base
        @normalize_base.setter
        def normalize_base(self,value):
            if not isinstance(value,int):
                raise ValueError('normalize_base must be an integer')
            if value <= 0:
                raise ValueError('normalize_base must be bigger than 0')
            self.__normalize_base = value

        @property
        def log_transform(self):
            return self.__log_transform
        @log_transform.setter
        def log_transform(self,value):
            if not isinstance(value,bool):
                raise ValueError('log_transform must be True or False')
            self.__log_transform = value

        @property
        def highly_variable_genes(self):
            return self.__highly_variable_genes
        @highly_variable_genes.setter
        def highly_variable_genes(self,value):
            if not isinstance(value,bool):
                raise ValueError('highly_variable_genes must be True or False')
            self.__highly_variable_genes = value

        @property
        def inner_graph_pca_n(self):
            return self.__inner_graph_pca_n
        @inner_graph_pca_n.setter
        def inner_graph_pca_n(self,value):
            if not isinstance(value,int):
                raise ValueError('inner_graph_pca_n must be an integer')
            if value <= 0:
                raise ValueError('inner_graph_pca_n must be bigger than 0')
            self.__inner_graph_pca_n = value

        @property
        def inner_graph_knn_n(self):
            return self.__inner_graph_knn_n
        @inner_graph_knn_n.setter
        def inner_graph_knn_n(self,value):
            if not isinstance(value,int):
                raise ValueError('inner_graph_knn_n must be an integer')
            if value <= 0:
                raise ValueError('inner_graph_knn_n must be bigger than 0')
            self.__inner_graph_knn_n = value

        @property
        def link_graph_pca_n(self):
            return self.__link_graph_pca_n
        @link_graph_pca_n.setter
        def link_graph_pca_n(self,value):
            if not isinstance(value,int):
                raise ValueError('link_graph_pca_n must be an integer')
            if value <= 0:
                raise ValueError('link_graph_pca_n must be bigger than 0')
            self.__link_graph_pca_n = value

        @property
        def link_graph_knn_n(self):
            return self.__link_graph_knn_n
        @link_graph_knn_n.setter
        def link_graph_knn_n(self,value):
            if not isinstance(value,int):
                raise ValueError('link_graph_knn_n must be an integer')
            if value <= 0:
                raise ValueError('link_graph_knn_n must be bigger than 0')
            self.__link_graph_knn_n = value

        @property
        def max_outgoing(self):
            return self.__max_outgoing
        @max_outgoing.setter
        def max_outgoing(self,value):
            if not isinstance(value,int):
                raise ValueError('max_outgoing must be be an integer')
            if value <= 0:
                raise ValueError('max_outgoing must be bigger than 0')
            self.__max_outgoing = value

        @property
        def min_score(self):
            return self.__min_score
        @min_score.setter
        def min_score(self,value):
            if not isinstance(value,(int,float)):
                raise ValueError('min_score must be numeric')
            if value <= 0 or value >= 1:
                raise ValueError('min_score must be between 0.0 and 1.0')  
            self.__min_score = value

        @property
        def min_cell_number(self):
            return self.__min_cell_number
        @min_cell_number.setter
        def min_cell_number(self,value):
            if not isinstance(value,int):
                raise ValueError('min_cell_number must be an integer')
            if value <= 0:
                raise ValueError('min_cell_number must be bigger than 0')
            self.__min_cell_number = value

    params=None
    timepoint_scdata_dict={}
    link_knn_graph=None
    link_cluster_graph=None
    filtered_cluster_node=None
    link_G=None
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
        self.timepoint_scdata_dict[self.__timepoint_scdata_num]=ad.AnnData(data)
        if timepoint_scdata_cluster != None:
            self.timepoint_scdata_dict[self.__timepoint_scdata_num].obs["scdata_cluster"]=[f"timepoint{self.__timepoint_scdata_num}_{c}" for c in list(timepoint_scdata_cluster)]

            self.timepoint_scdata_dict[self.__timepoint_scdata_num].uns["cluster_flag"]=True
        else:
            self.timepoint_scdata_dict[self.__timepoint_scdata_num].obs["scdata_cluster"]=[0]*self.timepoint_scdata_dict[self.__timepoint_scdata_num].n_obs
            self.timepoint_scdata_dict[self.__timepoint_scdata_num].uns["cluster_flag"]=False
        self.__timepoint_scdata_num+=1
    
    def __create_folder(self,folder_name):
        output_dir = self.params.output_dir
        output_name = self.params.output_name
        if not os.path.exists(output_dir):
            raise ValueError(f'{output_dir} : No such directory')
        elif not os.path.exists(output_dir+output_name):
            os.makedirs(output_dir+output_name)
        else:
            print(f"{output_dir+output_name} exists !")

        if not os.path.exists(output_dir+output_name+"/"+folder_name):
            os.makedirs(output_dir+output_name+"/"+folder_name)
        else:
            print(f"{output_dir+output_name+'/'+folder_name} exists !")



    @except_output
    @function_timer
    @params_filter(['pca_n','knn_n','resolution'])
    def cell_clusters(self,**kwargs):

        self.__create_folder("figures")

        pca_n = self.params.cell_cluster_pca_n = kwargs.setdefault('pca_n', self.params.cell_cluster_pca_n)
        knn_n = self.params.cell_cluster_knn_n = kwargs.setdefault('knn_n', self.params.cell_cluster_knn_n)
        resolution = self.params.cell_cluster_resolution=kwargs.setdefault("resolution",self.params.cell_cluster_resolution)
        
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
                sc.pp.normalize_per_cell(adata_copy, counts_per_cell_after=1e4)
                sc.pp.log1p(adata_copy)
                # high variable genes
                sc.pp.highly_variable_genes(adata_copy, min_mean=0.0125, max_mean=3, min_disp=0.5)
                adata_high = adata_copy[:, adata_copy.var['highly_variable']]
                # linear regression
                sc.pp.regress_out(adata_high, ['n_counts', 'percent_mito'])
                sc.pp.scale(adata_high, max_value=10)
                # pca
                sc.tl.pca(adata_high, n_comps=pca_n, svd_solver='arpack')
                # knn
                sc.pp.neighbors(adata_high, n_neighbors=knn_n, n_pcs=pca_n)
                sc.tl.louvain(adata_high, resolution=resolution)
                sc.tl.umap(adata_high)
                fig=sc.pl.umap(adata_high, color='louvain',return_fig=True)
                plt.show()
                ouput_path=self.params.output_dir+self.params.output_name+"/figures/"
                fig.savefig(ouput_path+f"timepoint{timepoint}_louvain_umap.pdf", bbox_inches='tight')

                adata.obs["scdata_cluster"]=[f"timepoint{timepoint}_cluster{int(c)+1}" for c in adata_high.obs["louvain"].tolist()]
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
    @params_filter(["percent_mito_cutoff"])
    def filter_dead_cell(self,**kwargs):
        percent_mito_cutoff=self.params.percent_mito_cutoff=kwargs.setdefault("percent_mito_cutoff",self.params.percent_mito_cutoff)

        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
            adata.obs['percent_mito'] = np.sum(adata[:, adata.var['mt']].X, axis=1) / np.sum(adata.X, axis=1)
            adata.obs['n_counts'] = adata.X.sum(axis=1)
            raw_cell_num=adata.n_obs
            adata=adata[adata.obs['percent_mito']<percent_mito_cutoff,:]
            filter_cell_num=adata.n_obs
            self.timepoint_scdata_dict[timepoint]=adata
            print(f'timepoint:{timepoint}')
            print(f'filtered out {raw_cell_num-filter_cell_num} cells that are detected in more than {percent_mito_cutoff} mito percent')
            print()

    @except_output
    @function_timer
    @params_filter(["min_cells"])
    def filter_lowcell_gene(self,**kwargs):
        min_cells=self.params.min_cells=kwargs.setdefault("min_cells",self.params.min_cells)
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            raw_gene_num=adata.n_vars
            sc.pp.filter_genes(adata, min_cells=min_cells)
            filter_gene_num=adata.n_vars
            print(f'timepoint:{timepoint}')
            print(f'filtered out {raw_gene_num-filter_gene_num} genes that are detected in less than {min_cells} cells')
            print()
    
    @except_output
    @function_timer
    @params_filter(["min_genes"])
    def filter_lowgene_cells(self,**kwargs):
        min_genes=self.params.min_genes=kwargs.setdefault("min_genes",self.params.min_genes)

        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            raw_cell_num=adata.n_obs
            sc.pp.filter_cells(adata, min_genes=min_genes)
            filter_cell_num=adata.n_obs
            print(f'timepoint:{timepoint}')
            print(f'filtered out {raw_cell_num-filter_cell_num} cells that are detected in less than {min_genes} genes')
            print()

    @except_output
    @function_timer
    @params_filter(["normalize_base"])
    def normalize_data(self,**kwargs):
        normalize_base=self.params.normalize_base=kwargs.setdefault("normalize_base",self.params.normalize_base)
        print(f'Normalize data to {normalize_base} count ...')
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            sc.pp.normalize_total(adata,target_sum=normalize_base)

    @except_output
    @function_timer
    def log_transform(self):
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            sc.pp.log1p(adata)

    @except_output
    @function_timer
    def highly_variable_genes(self):
        min_exp = 0
        min_exp_cnt = 10
        vargene = 2000
        minGeneCor = 0.2
        excludeGeneCor = 0.4
        cell_cycle = ['MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6','CDCA7','DTL','PRIM1','UHRF1','MLF1IP','HELLS','RFC2','RPA2','NASP','RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2','ATAD2','RAD51','RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2','USP1','CLSPN','POLA1','CHAF1B','BRIP1','E2F8','HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2','NUF2','CKS1B','MKI67','TMPO','CENPF','TACC3','FAM64A','SMC4','CCNB2','CKAP2L','CKAP2','AURKB','BUB1','KIF11','ANP32E','TUBB4B','GTSE1','KIF20B','HJURP','CDCA3','HN1','CDC20','TTK','CDC25C','KIF2C','RANGAP1','NCAPD2','DLGAP5','CDCA2','CDCA8','ECT2','KIF23','HMMR','AURKA','PSRC1','ANLN','LBR','CKAP5','CENPE','CTCF','NEK2','G2E3','GAS2L3','CBX5','CENPA']
        housekeep = ['RPS18','GAPDH','PGK1','PPIA','RPL13A','RPLP0','B2M','YWHAZ','SDHA','TFRC','RPA1','RPA2','RPAIN','RPE','RPL15','RPL15','RPL22','RPL32','RPL32','RPL35A','RPL4','RPL7L1','RPN1','RPN2','RPP30','RPP38','RPRD1A','RPS19BP1','RPS6KA5','RPS6KB1','RPUSD1','RPUSD2','RPUSD4']
        def __get_vscore(norm_data):
            '''
            '''
            min_mean = 0 # Exclude genes with average expression of this value or lower
            fit_percentile = 33 # Fit to this percentile of CV-vs-mean values
            error_wt = 1; # Exponent of fitting function. Value of 2 is least-squares (L2). Value of 1 is robust (L1).
            nBins = 50;

            # PREPARE DATA FOR FITTING
            mu_gene = norm_data.mean()
            idx = mu_gene>min_mean
            mu_gene = mu_gene[idx]
            FF_gene = norm_data.var()[idx]/mu_gene

            # Perform fit on log-transformed data:
            data_x = np.log(mu_gene)
            data_y = np.log(FF_gene/mu_gene)

            def runningquantile(x, y, p, nBins):
                x = x.sort_values()
                y = y[x.sort_values().index]
                dx = (x[-1]-x[0])/nBins
                xOut = np.arange(x[0]+dx/2, x[-1]-dx/2, dx)
                yOut = []
                for k,v in enumerate(xOut):
                    ind = (x>=(v-dx/2)) & (x<(v+dx/2))
                    if ind.sum()>0:
                        yOut.append(np.percentile(y[ind], p))
                    else:
                        if k>0:
                            yOut.append(yOut[k-1])
                        else:
                            yOut.append(np.nan)
                yOut = np.array(yOut)
                return(xOut, yOut)

            x,y = runningquantile(data_x, data_y, fit_percentile, nBins);
            
            FFhist, FFhist_vals = np.histogram(np.log10(FF_gene),200)
            c = np.exp(FFhist_vals[FFhist.argmax()]) # z = 1/( (1+a)(1+b) )
            c = max([1,c]) # Do not allow c to fall below 1.

            # FIT c=[CV_eff]^2
            errFun = lambda b: sum(np.abs(np.log(c*np.exp(-x) + b)-y)**error_wt)
            b = op.fmin(func=errFun, x0=0.1)[0]
            a = c/(1+b) - 1

            v_scores = FF_gene/((1+a)*(1+b) + b*mu_gene)
            v_scores = v_scores
            print('a:%.3f, b:%.3f, c:%.3f'%(a, b, c))
            return(v_scores)

        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            raw_gene_num=adata.n_vars
            adata.var["highly_variable_genes"]=False
            print(f"timepoint:{timepoint}")
            df=pd.DataFrame(adata.X)
            df.index=adata.obs_names
            df.columns=adata.var_names

            #基于threshold进行filter
            df = df.loc[:, df.apply(lambda x:sum(x>min_exp)>min_exp_cnt)]
            MTidx = df.columns.str.contains('MT-')
            df = df.loc[:, ~MTidx]
            # 基于CV进行filter
            vscore = __get_vscore(df)
            var_gene = vscore.sort_values(ascending=False)[:vargene].index.tolist()
            # 基于min correlation进行filter
            dfCor = df.loc[:, var_gene].corr()
            minCorIdx = (dfCor.abs()>minGeneCor).sum()>1
            # 删除housekeeping和cell cycle基因
            ccIdx = ~(dfCor.columns.str.upper()).isin(cell_cycle)
            hkIdx = ~(dfCor.columns.str.upper()).isin(housekeep)
            ccCorIdx = ~((dfCor.loc[:, ~ccIdx].abs()>excludeGeneCor).sum(axis=1)>0)
            hkCorIdx = ~((dfCor.loc[:, ~hkIdx].abs()>excludeGeneCor).sum(axis=1)>0)
            # keep idx
            keepIdx = minCorIdx & ccIdx & hkIdx & ccCorIdx & hkCorIdx
            adata.var["highly_variable_genes"][var_gene]=keepIdx
            self.timepoint_scdata_dict[timepoint]=adata[:,adata.var["highly_variable_genes"]]
            hv_gene_num=self.timepoint_scdata_dict[timepoint].n_vars
            print(f"filtered out {hv_gene_num} highly variable genes")
            print()
    
    @except_output
    @function_timer
    @params_filter(["pca_n","k"])
    def get_knn_inner(self,**kwargs):
        pca_n=self.params.inner_graph_pca_n=kwargs.setdefault("pca_n",self.params.inner_graph_pca_n)
        k=self.params.inner_graph_knn_n=kwargs.setdefault("k",self.params.inner_graph_knn_n)

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
            nbrs = NearestNeighbors(n_neighbors=k+1, metric='correlation', n_jobs=-2)
            nbrs.fit(df_zscore_pca)
            dists, ind = nbrs.kneighbors(df_zscore_pca) # n*k_init dist & index
            adj = nbrs.kneighbors_graph(df_zscore_pca, mode='distance') # n*n dist
            adata.obs["cell_id"]=sample_name
            adata.obsm["inner_dists"]=dists[:, 1:]
            adata.obsm["inner_ind"]=sample_name[ind[:, 1:]]
            adata.obsm['inner_distance']=adj
            print()

    @except_output
    @function_timer
    @params_filter(["pca_n","k"])
    def get_knn_link(self,**kwargs):
        pca_n=self.params.link_graph_pca_n=kwargs.setdefault("pca_n",self.params.link_graph_pca_n)
        k=self.params.link_graph_knn_n=kwargs.setdefault("k",self.params.link_graph_knn_n)
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
            nbrs = NearestNeighbors(n_neighbors=k+1, metric='correlation', n_jobs=-2)
            nbrs.fit(df_end_pca)
            dists, ind = nbrs.kneighbors(df_start_pca) # n*k_init dist & index
            adj = nbrs.kneighbors_graph(df_start_pca, mode='distance') # n*n dist
            adata_start.obsm["link_later_dists"]=dists[:, 1:]
            adata_start.obsm["link_later_ind"]=adata_end.obs["cell_id"][ind[:, 1:]]
            adata_start.obsm['link_later_distance']=adj
            # KNN distance (end)
            nbrs = NearestNeighbors(n_neighbors=k+1, metric='correlation', n_jobs=-2)
            nbrs.fit(df_start_pca)
            dists, ind = nbrs.kneighbors(df_end_pca) # n*k_init dist & index
            adj = nbrs.kneighbors_graph(df_end_pca, mode='distance') # n*n dist
            adata_end.obsm["link_previous_dists"]=dists[:, 1:]
            adata_end.obsm["link_previous_ind"]=adata_start.obs["cell_id"][ind[:, 1:]]
            adata_end.obsm['link_previous_distance']=adj

    @except_output
    @function_timer
    def get_knn_graph(self):
        #inner graph
        for (timepoint,adata) in self.timepoint_scdata_dict.items():

            print(f"timepoint:{timepoint}")
            node1 = np.repeat(adata.obs["cell_id"], adata.obsm["inner_dists"].shape[1])
            node2 = adata.obsm["inner_ind"].flatten()
            Distance = adata.obsm["inner_dists"].flatten()
            Distance_normalize = (adata.obsm["inner_dists"]/adata.obsm["inner_dists"][:,0][:,None]).flatten()
            Distance_zscore = zscore(Distance)
            
            df_graph=pd.DataFrame({
                'Node1':node1,
                'Node2':node2,
                'Distance':Distance,
                'Distance_normalize':Distance_normalize,
                'Distance_zscore':Distance_zscore
                })
            adata.uns['inner_knn_graph']=df_graph


        #link graph
        for timepoint in list(self.timepoint_scdata_dict.keys())[:-1]:
            print(f"timepoint between {timepoint} and {timepoint+1} ")
            adata_start=self.timepoint_scdata_dict[timepoint]
            node1 = np.repeat(adata_start.obs["cell_id"], adata_start.obsm["link_later_dists"].shape[1])
            node2 = adata_start.obsm["link_later_ind"].flatten()
            Distance = adata_start.obsm["link_later_dists"].flatten()
            Distance_normalize = (adata_start.obsm["link_later_dists"]/adata_start.obsm["link_later_dists"][:,0][:,None]).flatten()
            Distance_zscore = zscore(Distance)

            df_graph=pd.DataFrame({
                'Node1':node1,
                'Node2':node2,
                'Distance':Distance,
                'Distance_normalize':Distance_normalize,
                'Distance_zscore':Distance_zscore
                })
            adata_start.uns['link_later_knn_graph']=df_graph.reset_index(drop=True)

            adata_end=self.timepoint_scdata_dict[timepoint+1]
            node1 = np.repeat(adata_end.obs["cell_id"], adata_end.obsm["link_previous_dists"].shape[1])
            node2 = adata_end.obsm["link_previous_ind"].flatten()
            Distance = adata_end.obsm["link_previous_dists"].flatten()
            Distance_normalize = (adata_end.obsm["link_previous_dists"]/adata_end.obsm["link_previous_dists"][:,0][:,None]).flatten()
            Distance_zscore = zscore(Distance)
            df_graph=pd.DataFrame({
                'Node1':node2,
                'Node2':node1,
                'Distance':Distance,
                'Distance_normalize':Distance_normalize,
                'Distance_zscore':Distance_zscore                
                })
            adata_end.uns['link_previous_knn_graph']=df_graph

    @except_output
    @function_timer
    def filter_knn_graph(self):
        knn_graph=pd.DataFrame()
        #inner graph
        def Maximum_Score(df):
            Q1=df.quantile(0.25)
            Q3=df.quantile(0.75)
            return Q3+(Q3-Q1)*1.5

        for (timepoint,adata) in self.timepoint_scdata_dict.items():

            print(f"timepoint:{timepoint}")
            df_graph=adata.uns['inner_knn_graph']
            inner_max_normalize_keep_id=df_graph.loc[:,'Distance_normalize']<=Maximum_Score(df_graph.loc[:,'Distance_normalize'])
            inner_max_zscore_keep_id=df_graph.loc[:,'Distance_zscore']<=Maximum_Score(df_graph.loc[:,'Distance_zscore'])
            keep_id=inner_max_normalize_keep_id & inner_max_zscore_keep_id
            df_graph=df_graph.loc[keep_id,:].reset_index(drop=True)

            duplicate_id=df_graph.loc[:,['Node1', 'Node2']].apply(lambda x:x.sort_values().str.cat(), axis=1).duplicated()
            
            if True:
                df_graph=df_graph.loc[duplicate_id,:].reset_index(drop=True)
            else:
                df_graph=df_graph.loc[~duplicate_id,:].reset_index(drop=True)

            adata.uns['inner_knn_graph']=df_graph

        #link graph
        # link_max_normalize=3
        # link_max_zscore=0

        for timepoint in list(self.timepoint_scdata_dict.keys())[:-1]:
            print(f"timepoint between {timepoint} and {timepoint+1} ")
            adata_start=self.timepoint_scdata_dict[timepoint]
            df_graph=adata_start.uns['link_later_knn_graph']
            link_max_normalize_keep_id=df_graph.loc[:,'Distance_normalize']<=Maximum_Score(df_graph.loc[:,'Distance_normalize'])
            link_max_zscore_keep_id=df_graph.loc[:,'Distance_zscore']<=Maximum_Score(df_graph.loc[:,'Distance_zscore'])
            keep_id=link_max_normalize_keep_id & link_max_zscore_keep_id
            df_graph=df_graph.loc[keep_id,:].reset_index(drop=True)


            adata_end=self.timepoint_scdata_dict[timepoint+1]
            df_graph2=adata_end.uns['link_previous_knn_graph']
            link_max_normalize_keep_id=df_graph2.loc[:,'Distance_normalize']<=Maximum_Score(df_graph2.loc[:,'Distance_normalize'])
            link_max_zscore_keep_id=df_graph2.loc[:,'Distance_zscore']<=Maximum_Score(df_graph2.loc[:,'Distance_zscore'])
            keep_id=link_max_normalize_keep_id & link_max_zscore_keep_id
            df_graph2=df_graph2.loc[keep_id,:].reset_index(drop=True)

            if True :
                df_graph_new=pd.concat([df_graph2,df_graph])
                df_graph_new=df_graph_new.loc[df_graph_new.duplicated(['Node1','Node2']),:].reset_index(drop=True)
            else:
                df_graph_new=pd.concat([df_graph2,df_graph])
                df_graph_new=df_graph_new.loc[~df_graph_new.duplicated(['Node1','Node2']),:].reset_index(drop=True)

            knn_graph=pd.concat([knn_graph,df_graph_new])
        self.knn_graph=knn_graph

    def __get_inner_cgKNN(self,adata):
        knn_graph=adata.uns["inner_knn_graph"]
        node_cluster=adata.obs["scdata_cluster"]
        node_cluster.index=adata.obs["cell_id"]


        knn_graph.loc[:,'Node1_cluster'] = knn_graph.loc[:,'Node1'].map(node_cluster.to_dict())
        knn_graph.loc[:,'Node2_cluster'] = knn_graph.loc[:,'Node2'].map(node_cluster.to_dict())
        #去除相同类型的细胞之间的连接
        nonself_idx = knn_graph.apply(lambda x:x['Node1_cluster']!=x['Node2_cluster'], axis=1)
        knn_graph=knn_graph.loc[nonself_idx,:]
        if knn_graph.shape[0]!=0:
            #设定边的唯一id
            knn_graph["edge"]=knn_graph.loc[:,["Node1_cluster","Node2_cluster"]].apply(lambda x:x.sort_values().str.cat(sep="->"), axis=1)
            knn_graph=knn_graph["edge"].str.split("->", expand=True)
            knn_graph.columns=["Node1_cluster","Node2_cluster"]
            
            knn_graph["edge_num"]=1
            knn_graph_egde=knn_graph.groupby(['Node1_cluster','Node2_cluster'],as_index=False)['edge_num'].count()
            knn_graph_egde['Node1_cluster_cell_num']=knn_graph_egde['Node1_cluster'].map(adata.obs['scdata_cluster'].value_counts().to_dict())
            knn_graph_egde["edge_num_per_cell"]=knn_graph_egde["edge_num"]/knn_graph_egde["Node1_cluster_cell_num"]

            tmp_dict=knn_graph_egde.groupby(['Node1_cluster'],as_index=False)['edge_num_per_cell'].sum()
            tmp_dict.index=tmp_dict["Node1_cluster"]
            knn_graph_egde['outgoing_edge_num']=knn_graph_egde['Node1_cluster'].map(tmp_dict['edge_num_per_cell'].to_dict())
            knn_graph_egde["score"]=knn_graph_egde["edge_num_per_cell"]/knn_graph_egde["outgoing_edge_num"]
        else:
            knn_graph_egde=knn_graph

        adata.uns["cluster_graph"]=knn_graph_egde



    def __get_link_cgKNN(self,knn_graph,node_cluster,node_cluster_cell_num):
        knn_graph.loc[:,'Node1_cluster'] = knn_graph.loc[:,'Node1'].map(node_cluster.to_dict())
        knn_graph.loc[:,'Node2_cluster'] = knn_graph.loc[:,'Node2'].map(node_cluster.to_dict())
        
        knn_graph["edge"]=knn_graph.loc[:,["Node1_cluster","Node2_cluster"]].apply(lambda x:x.sort_values().str.cat(sep="->"), axis=1)
        knn_graph=knn_graph["edge"].str.split("->", expand=True)
        knn_graph.columns=["Node1_cluster","Node2_cluster"]

        knn_graph["edge_num"]=1
        knn_graph_egde=knn_graph.groupby(['Node1_cluster','Node2_cluster'],as_index=False)['edge_num'].count()
        knn_graph_egde['Node1_cluster_cell_num']=knn_graph_egde['Node1_cluster'].map(node_cluster_cell_num)
        knn_graph_egde["edge_num_per_cell"]=knn_graph_egde["edge_num"]/knn_graph_egde["Node1_cluster_cell_num"]

        tmp_dict=knn_graph_egde.groupby(['Node1_cluster'],as_index=False)['edge_num_per_cell'].sum()
        tmp_dict.index=tmp_dict["Node1_cluster"]
        knn_graph_egde['outgoing_edge_num']=knn_graph_egde['Node1_cluster'].map(tmp_dict['edge_num_per_cell'].to_dict())
        knn_graph_egde["score"]=knn_graph_egde["edge_num_per_cell"]/knn_graph_egde["outgoing_edge_num"]

        return knn_graph_egde

    @except_output
    @function_timer
    def get_cluster_trajectory(self):
        #inner graph
        node_cluster=pd.DataFrame()
        node_cluster_cell_num={}
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            print(f"timepoint:{timepoint}")
            node_cluster=pd.concat([node_cluster,adata.obs[["cell_id","scdata_cluster"]]])
            node_cluster.index=node_cluster["cell_id"]
            node_cluster_cell_num={**node_cluster_cell_num,**adata.obs['scdata_cluster'].value_counts().to_dict()}
            self.__get_inner_cgKNN(adata)

        #link graph
        knn_graph_egde=self.__get_link_cgKNN(self.knn_graph,node_cluster["scdata_cluster"],node_cluster_cell_num)

        self.cluster_graph=knn_graph_egde.sort_values(by='score',ascending=False).reset_index(drop=True)

    def __get_nxG(self,knn_graph_egde, allnodes,min_score):
        G=nx.DiGraph()
        for node in allnodes:
            G.add_node(node,timepoint=str(node).split("_")[0])
        for index,row in knn_graph_egde.iterrows():
            if (row["Node1_cluster"] in allnodes) & (row["Node2_cluster"] in allnodes) & (row["score"]>min_score) :
                G.add_edge(row["Node1_cluster"],row["Node2_cluster"],score=row["score"],start_timepoint=G.nodes[row["Node1_cluster"]]["timepoint"])
        return G

    def if_ZeroDivisionError(exception):
        return isinstance(exception, ZeroDivisionError)

    @retry(retry_on_exception=if_ZeroDivisionError)
    def __force_directed_layout(self,G,timepoint,verbose=True, iterations=50):
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
            scalingRatio=1,
            strongGravityMode=False,
            gravity=15.0,
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
        plt.show()

        ouput_path=self.params.output_dir+self.params.output_name+"/figures/"
        f.savefig(ouput_path+f"timepoint{timepoint}_fa.pdf", bbox_inches='tight')

        return positions

    @except_output
    @function_timer
    @params_filter(["topN","min_score","min_cell_number"])
    def get_knn_nxG(self,**kwargs):
        topN=self.params.max_outgoing=kwargs.setdefault("topN",self.params.max_outgoing)
        min_score=self.params.min_score=kwargs.setdefault("min_score",self.params.min_score)
        min_cell_number=self.params.min_cell_number=kwargs.setdefault("min_cell_number",self.params.min_cell_number)

        filtered_node_cluster=np.array([])
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            print(f"timepoint:{timepoint}")
            knn_graph_egde=adata.uns['cluster_graph'].loc[:,["Node1_cluster","Node2_cluster","score"]]
            allnodes=adata.uns["cluster_set"]
            allnodes_counts=adata.uns["cluster_counts"]
            filtered_nodes=[node for node in allnodes if allnodes_counts[node]>min_cell_number]

            filtered_node_cluster=np.hstack([filtered_node_cluster,filtered_nodes])
            
            inner_G = self.__get_nxG(knn_graph_egde, filtered_nodes,min_score).to_undirected()
            adata.uns["inner_G"]=inner_G
            
            fa_cord_all = self.__force_directed_layout(G=inner_G,timepoint=timepoint)
            adata.uns["fa_cord"]=fa_cord_all

        knn_graph_egde=self.cluster_graph.groupby(["Node1_cluster"],as_index=False).head(topN)
        self.filtered_cluster_node=filtered_node_cluster
        link_G=self.__get_nxG(knn_graph_egde, filtered_node_cluster,min_score)
        self.link_G=link_G
    
    @except_output
    @function_timer
    #@params_filter()
    def draw_nxG(self,**kwargs):

        ZORDER=0.1
        TRANS=25 
        SQ=np.array([[-2.5,12.5],[12.5,12.5],[12.5,-2.5],[-2.5,-2.5],[-2.5,12.5]]).T
        T=np.array([[1,0],[1,1.5]])
        SQ=T@SQ
        CCCOLOR=["#CE0013","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD","#16557A"]
        number_of_scdata=len(self.timepoint_scdata_dict)
        fig = plt.figure(figsize=(number_of_scdata*10,30))
        ax = fig.add_subplot(1,1,1)
        all_pos={}
        all_pos2={}
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            
            ax.fill(SQ[0]+TRANS*timepoint,SQ[1],facecolor="#979a9a",edgecolor='gray',lw=3,alpha=0.65,zorder=ZORDER*timepoint)
            ax.plot(SQ[0]+TRANS*timepoint,SQ[1],'black',lw=3,zorder=ZORDER*timepoint)
            plt.text(5+TRANS*timepoint,-5, f'$ t_{timepoint} $',size=40, color="black",weight="bold",verticalalignment="center",horizontalalignment="center")
            
            G=adata.uns["inner_G"]
            pos=adata.uns["fa_cord"].copy()
            pos2={}
            pos3={}
            pos4={}
            width_list=[]
            for (u, v, wt) in G.edges.data('score'):
                width_list.append(wt*5)
                
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
                labels2[nodes]=str(nodes).split("_")[1]
                x_pos=TRANS*timepoint-2
                y_pos=np.linspace(-10,-60,(G.number_of_nodes()+2))[1:-1][i]
                pos2[nodes]=(x_pos,y_pos)
                pos3[nodes]=(x_pos+1,y_pos)
                pos4[nodes]=(x_pos-0.3,y_pos)
            all_pos2={**all_pos2,**pos4}
            all_pos={**all_pos,**pos}
            #pos
            nodes=nx.draw_networkx_nodes(G, pos, ax=ax, node_size=800, node_color=CCCOLOR[timepoint], alpha=1)
            nodes.set_zorder((ZORDER+0.05)*timepoint)
            edges=nx.draw_networkx_edges(G, pos, ax=ax, width=width_list,edge_color="black", alpha=1)
            if G.number_of_edges() != 0:
                edges.set_zorder((ZORDER+0.015)*timepoint)
            nx.draw_networkx_labels(G, pos,labels, ax=ax, font_size=20, font_color="white",font_weight="bold")
            
            #pos2
            nodes=nx.draw_networkx_nodes(G, pos2, ax=ax, node_size=1500, node_color=CCCOLOR[timepoint], alpha=1)
            nodes.set_zorder((ZORDER+0.05)*timepoint)
            nx.draw_networkx_labels(G, pos2,labels, ax=ax, font_size=20, font_color="white",font_weight="bold")
            nx.draw_networkx_labels(G, pos3,labels2,horizontalalignment="left", ax=ax, font_size=30)
        
        link_G=self.link_G

        width_list=[]
        for (u, v, wt) in link_G.edges.data('score'):
            width_list.append(wt*5)
        timepoint_list=[]
        for (u, v, tp) in link_G.edges.data('start_timepoint'):
            timepoint_list.append(tp)    
        link_edges=nx.draw_networkx_edges(link_G, all_pos, ax=ax,width=width_list,edge_color="#0060A8",arrowstyle="->", arrowsize=50,alpha=1)
        for i in range(len(timepoint_list)):
            tp=timepoint_list[i].split("timepoint")[1]
            link_edges[i].set_zorder((ZORDER+0.025)*int(tp))
        link_edges2=nx.draw_networkx_edges(link_G, all_pos2, ax=ax,width=width_list,edge_color="#0060A8",arrowstyle="->", arrowsize=50,alpha=1)
        for i in range(len(timepoint_list)):
            tp=timepoint_list[i].split("timepoint")[1]
            link_edges2[i].set_zorder((ZORDER+0.025)*int(tp))
        plt.tight_layout()
        plt.show()
        ouput_path=self.params.output_dir+self.params.output_name+"/"
        fig.savefig(ouput_path+"cstreet_result.pdf", bbox_inches='tight')

    def run_cstreet(self):

        if True :
            self.cell_clusters()
        
        if self.params.filter_dead_cell == True:
            self.filter_dead_cell()
        
        if self.params.filter_lowcell_gene == True:
            self.filter_lowcell_gene()
        
        if self.params.filter_lowgene_cells == True:
            self.filter_lowgene_cells()
        
        if self.params.normalize == True :
            self.normalize_data()

        if self.params.log_transform == True :
            self.log_transform()

        if self.params.highly_variable_genes == True :
            self.highly_variable_genes()
        
        if True :
            self.get_knn_inner()

            self.get_knn_link()

            self.get_knn_graph()

            self.filter_knn_graph()
            
            self.get_cluster_trajectory()
            
            self.get_knn_nxG()
            
            self.draw_nxG()










