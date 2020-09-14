#!python3.6
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

#定义cstreet对象

class CStreetData(object):
    """docstring for CStreetData"""
    params={
    
    #Step1:cell_cluster# 
    
    "cell_cluster":True,
    "cell_cluster_pca_n":10,
    "cell_cluster_knn_n":15,
    "cell_cluster_resolution":1,

    #Step2:gene and cell filter#

    "filter_dead_cell":True,
    "percent_mito_cutoff":0.2,
    
    "filter_lowcell_gene":True,
    "min_cells":3,
    
    "filter_lowgene_cells":True,
    "min_genes":200,
    
    #Step3:normalize#
    "normalize":True,
    "normalize_base":10000,
    "log_transform":True,

    #Step4:get HVG#
    "highly_variable_genes":True,

    #Step5:get_graph#
    "inner_graph_pca_n":10,
    "inner_graph_knn_n":15,

    "link_graph_pca_n":10,
    "link_graph_knn_n":15,

    #Step6:plot graph#
    "max_outgoing":10,
    "min_score":0.1
    }

    timepoint_scdata_dict={}
    link_knn_graph=None
    link_cluster_graph=None
    all_cluster_node=None
    link_G=None
    __timepoint_scdata_num=0

    def __init__(self): 
        super(CStreetData, self).__init__()
        
    def add_new_timepoint_scdata(self,timepoint_scdata,timepoint_scdata_cluster=None):
        self.timepoint_scdata_dict[self.__timepoint_scdata_num]=ad.AnnData(pd.DataFrame(timepoint_scdata))
        if timepoint_scdata_cluster != None:
            self.timepoint_scdata_dict[self.__timepoint_scdata_num].obs["scdata_cluster"]=[f"timepoint{self.__timepoint_scdata_num}_{c}" for c in list(timepoint_scdata_cluster)]

            self.timepoint_scdata_dict[self.__timepoint_scdata_num].uns["cluster_flag"]=True
        else:
            self.timepoint_scdata_dict[self.__timepoint_scdata_num].obs["scdata_cluster"]=[0]*self.timepoint_scdata_dict[self.__timepoint_scdata_num].n_obs
            self.timepoint_scdata_dict[self.__timepoint_scdata_num].uns["cluster_flag"]=False
        self.__timepoint_scdata_num+=1

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

    def cell_clusters(self,**kwargs):

        pca_n = kwargs.setdefault('pca_n', self.params["cell_cluster_pca_n"])
        knn_n=kwargs.setdefault('knn_n', self.params["cell_cluster_knn_n"])
        resolution=kwargs.setdefault("resolution",self.params["cell_cluster_resolution"])

        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            print(f"timepoint:{timepoint}")
            if adata.uns['cluster_flag'] :
                print(f"clusters have been gave")
                continue
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
            sc.pl.umap(adata_high, color='louvain')
            adata.obs["scdata_cluster"]=[f"timepoint{timepoint}_cluster{c}" for c in adata_high.obs["louvain"].tolist()]
            adata.uns["cluster_flag"]=True


    def filter_dead_cell(self,**kwargs):
        percent_mito_cutoff=kwargs.setdefault("percent_mito_cutoff",self.params["percent_mito_cutoff"])

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

    def filter_lowcell_gene(self,**kwargs):
        min_cells=kwargs.setdefault("min_cells",self.params["min_cells"])
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            raw_gene_num=adata.n_vars
            sc.pp.filter_genes(adata, min_cells=min_cells)
            filter_gene_num=adata.n_vars
            print(f'timepoint:{timepoint}')
            print(f'filtered out {raw_gene_num-filter_gene_num} genes that are detected in less than {min_cells} cells')
            print()
    def filter_lowgene_cells(self,**kwargs):
        min_genes=kwargs.setdefault("min_genes",self.params["min_genes"])

        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            raw_cell_num=adata.n_obs
            sc.pp.filter_cells(adata, min_genes=min_genes)
            filter_cell_num=adata.n_obs
            print(f'timepoint:{timepoint}')
            print(f'filtered out {raw_cell_num-filter_cell_num} cells that are detected in less than {min_genes} genes')
            print()

    def normalize_data(self,**kwargs):
        normalize_base=kwargs.setdefault("normalize_base",self.params["normalize_base"])
        print(f'Normalize data to {normalize_base} count ...')
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            sc.pp.normalize_total(adata,target_sum=normalize_base)

    def log_transform(self):
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            sc.pp.log1p(adata)

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
    def get_knn_inner(self,**kwargs):
        pca_n=kwargs.setdefault("pca_n",self.params["inner_graph_pca_n"])
        k=kwargs.setdefault("k",self.params["inner_graph_knn_n"])

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

    def get_knn_link(self,**kwargs):
        pca_n=kwargs.setdefault("pca_n",self.params["link_graph_pca_n"])
        k=kwargs.setdefault("k",self.params["link_graph_knn_n"])
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
        adata.uns["cluster_set"]=node_cluster.unique()
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

    def __get_nxG(self,knn_graph_egde, allnodes):
        G=nx.DiGraph()
        for node in allnodes:
            G.add_node(node,timepoint=str(node).split("_")[0])
        for index,row in knn_graph_egde.iterrows():
            G.add_edge(row["Node1_cluster"],row["Node2_cluster"],score=row["score"],start_timepoint=G.nodes[row["Node1_cluster"]]["timepoint"])
        return G

    def __force_directed_layout(self,G,verbose=True, iterations=50):
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
        ax.set(xlim=[-1, 11],ylim=[-1, 11],title='An Example Axes')
        nx.draw_networkx_nodes(G, positions, ax=ax, node_size=800, node_color="blue", alpha=0.5)
        nx.draw_networkx_edges(G, positions, ax=ax, width=10,edge_color="green", alpha=0.5)
        nx.draw_networkx_labels(G, positions, ax=ax, font_size=20)
        plt.show()
        #f.savefig(pdf_name, bbox_inches='tight')

        return positions

    def get_knn_nxG(self):
        all_node_cluster=np.array([])
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            print(f"timepoint:{timepoint}")
            knn_graph_egde=adata.uns['cluster_graph'].loc[:,["Node1_cluster","Node2_cluster","score"]]
            allnodes=adata.uns["cluster_set"]
            all_node_cluster=np.hstack([all_node_cluster,allnodes])
            
            inner_G = self.__get_nxG(knn_graph_egde, allnodes).to_undirected()
            adata.uns["inner_G"]=inner_G
            
            fa_cord_all = self.__force_directed_layout(G=inner_G)
            adata.uns["fa_cord"]=fa_cord_all

        knn_graph_egde=self.cluster_graph
        self.all_cluster_node=all_node_cluster
        link_G=self.__get_nxG(knn_graph_egde, all_node_cluster)
        self.link_G=link_G
    
    def draw_nxG(self,**kwargs):
        topN=kwargs.setdefault("",self.params["max_outgoing"])
        min_score=kwargs.setdefault("min_score",self.params["min_score"])

        TRANS=16 
        sq=np.array([[-2.5,12.5],[12.5,12.5],[12.5,-2.5],[-2.5,-2.5],[-2.5,12.5]]).T
        T=np.array([[1,0],[1,1.5]])
        sq=T@sq
        hhcolor=["orangered","green","red","orangered","red","green"]

        fig = plt.figure(figsize=(40,60))
        ax = fig.add_subplot(1,1,1)
        all_pos={}
        all_pos2={}
        for (timepoint,adata) in self.timepoint_scdata_dict.items():
            
            ax.fill(sq[0]+TRANS*timepoint,sq[1],facecolor="gray",edgecolor='black',lw=3,alpha=0.55,zorder=timepoint)
            G=adata.uns["inner_G"]
            pos=adata.uns["fa_cord"].copy()
            pos2={}
            pos3={}
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
            for i,nodes in enumerate(adata.uns["cluster_set"]):
                labels[nodes]=str(i)
                labels2[nodes]=str(nodes).split("_")[1]
                pos2[nodes]=(5+TRANS*timepoint,np.linspace(-10,-70,G.number_of_nodes())[i])
                pos3[nodes]=(5+TRANS*timepoint,np.linspace(-10,-70,G.number_of_nodes())[i]-2)
            all_pos2={**all_pos2,**pos2}
            all_pos={**all_pos,**pos}
            #pos
            nodes=nx.draw_networkx_nodes(G, pos, ax=ax, node_size=800, node_color=hhcolor[timepoint], alpha=1)
            nodes.set_zorder(timepoint+0.5)
            edges=nx.draw_networkx_edges(G, pos, ax=ax, width=width_list,edge_color="black", alpha=1)
            if G.number_of_edges() != 0:
                edges.set_zorder(timepoint+0.15)
            nx.draw_networkx_labels(G, pos,labels, ax=ax, font_size=20)
            #pos2
            nodes=nx.draw_networkx_nodes(G, pos2, ax=ax, node_size=1200, node_color=hhcolor[timepoint], alpha=1)
            nodes.set_zorder(timepoint+0.5)
            nx.draw_networkx_labels(G, pos2,labels, ax=ax, font_size=20)
            nx.draw_networkx_labels(G, pos3,labels2,verticalalignment="bottom", ax=ax, font_size=30)

        knn_graph_egde=self.cluster_graph.groupby(["Node1_cluster"],as_index=False).head(topN)
        knn_graph_egde=knn_graph_egde.loc[knn_graph_egde["score"]>=min_score,:]

        all_node_cluster=self.all_cluster_node
        link_G=self.__get_nxG(knn_graph_egde, all_node_cluster)

        width_list=[]
        for (u, v, wt) in link_G.edges.data('score'):
            width_list.append(wt*5)
        timepoint_list=[]
        for (u, v, tp) in link_G.edges.data('start_timepoint'):
            timepoint_list.append(tp)    
        link_edges=nx.draw_networkx_edges(link_G, all_pos, ax=ax,width=width_list,edge_color="blue",arrowstyle="->", arrowsize=25,alpha=1)
        for i in range(len(timepoint_list)):
            tp=timepoint_list[i].split("timepoint")[1]
            link_edges[i].set_zorder(int(tp)+0.25)
        link_edges2=nx.draw_networkx_edges(link_G, all_pos2, ax=ax,width=width_list,edge_color="blue",arrowstyle="->", arrowsize=25,alpha=1)
        for i in range(len(timepoint_list)):
            tp=timepoint_list[i].split("timepoint")[1]
            link_edges2[i].set_zorder(int(tp)+0.25)
        plt.show()

    def run_cstreet(self):
        params_dict=self.params

        if params_dict["cell_cluster"] == True :
            print("cell clusters:")
            self.cell_clusters()
            print("---------done!---------")
        
        if params_dict['filter_dead_cell'] == True:
            print("filter_dead_cell")
            self.filter_dead_cell()
            print("---------done!---------")
        if params_dict["filter_lowcell_gene"]:
            print("filter_lowcell_gene")
            self.filter_lowcell_gene()
            print("---------done!---------")
        if params_dict["filter_lowgene_cells"] == True:
            print("filter_lowgene_cells")
            self.filter_lowgene_cells()
            print("---------done!---------")
        
        if params_dict["normalize"] == True :
            self.normalize_data()
        if params_dict["log_transform"] == True :
            self.log_transform()

        if params_dict["highly_variable_genes"] == True :
            self.highly_variable_genes()
        if True :
            print("get_knn_inner")
            self.get_knn_inner()
            print("---------done!---------")

            print("get_knn_link")
            self.get_knn_link()
            print("---------done!---------")

            print("get_knn_graph")
            self.get_knn_graph()
            print("---------done!---------")            

            print("filter_knn_graph")
            self.filter_knn_graph()
            print("---------done!---------")
            
            print("fa2")
            self.get_cluster_trajectory()
            print("---------done!---------")
            
            print("draw_nxG")
            self.get_knn_nxG()
            self.draw_nxG()
            print("---------done!---------")










