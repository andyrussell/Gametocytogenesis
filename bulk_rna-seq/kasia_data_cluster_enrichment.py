import os
import pandas as pd
from scipy.stats import hypergeom ## for hypergeomtric distributions
import numpy as np
import statsmodels.stats.multitest as mtest ## for multiple hypothesis testing
current_folder = os.getcwd() ## current directory path



def apply_enrich(tmp,sex_gene_list,outfile):
    ## groupby along Phenotyepes
    #First argumnet: temperory dataframe
    #Second argument: list of male-specific genes
    #Third argumnet: file where we store enrichment analysis


    M=len(tmp.index) # total number of genes
    pvals=[] # list for storing pvalues
    n_list=[] # store number of genes in cluster
    x_list=[] ## store number of sex-specific genes in cluster
    genes=[] ## we store genes (sex-specifc gene ids)
    clusters=[] ## store cluster number

    ## group clusters
    for k,v in tmp.groupby([ 'Cluster']).indices.items():
        ## print(k, len(v))
        clusters.append(k) ## store cluster number
        clust_genes=set(tmp.index[v].to_list()) ## genes in the cluster

        sex_gene_in_clust= set(sex_gene_list)& clust_genes ## how many sex-specific genes are found in cluster
        ## how many are in phenoype
        M=M ## total genes
        N=len(sex_gene_list) ## sex-specific genes
        n=len(clust_genes) ## total number of cluster genes

        rv = hypergeom(M, n, N) ## hypergeometric  distributions
        x=len(sex_gene_in_clust) ## number of sex-specific genes in cluster

        ## for p-value we need to get summation of all probalities greater than > x
        pval= 1-rv.cdf(x)
        pvals.append(pval)
        n_list.append(n)
        x_list.append(x)
        genes.append(', '.join(sex_gene_in_clust))
    fdr=mtest.multipletests(pvals, alpha=0.05, method='fdr_bh')
    enrich={'Cluster':clusters, 'pvals':pvals,'cluster_size':n_list,'# of genes found in cluster':x_list,'genes':genes,'FDR':fdr[1]}

    ## create enrichment table
    enrich_df=pd.DataFrame.from_dict(enrich)
    enrich_soreted=enrich_df.sort_values(by=['pvals'],ascending=True)

    enrich_soreted.to_csv(outfile,sep='\t',index=None)






def enrich_cluster_final():
    ''' Calculate p-value and FDR using enrichment analysis. We used hypergeometric test'''

    ### input file (provided by Oliver): this file has 3 columns
    # First column: 'new gene ID'
    # Second column: 'Cluster'
    # Third column: 'cons phenotype'
    ##
    input_file=current_folder+'/Kasia new cluster enrichment for Vikash.txt'
    df=pd.read_csv(input_file,sep='\t')## input file as a dataframe
    df=df.fillna('NA') ## nan is fill with NA
    # df['consensus phenotype']=df['consensus phenotype'].str.replace('None','NA')
    tmp=df.copy() ## this is the temperory data frame same as a original dataframe
    tmp.set_index('new gene ID',inplace=True) ## make gene ids as index

    ## seprate the dataframes according sex-specific phenotypes
    male_df=tmp[(tmp['cons phenotype']== 'Male')]
    female_df=tmp[(tmp['cons phenotype']== 'Female')]
    both_df=tmp[(tmp['cons phenotype']== 'Both')]


    ## enrichment analysis for male
    outfile=current_folder+'/Male_Cluster_enrichment_final_new.txt'
    #First argumnet: temperory dataframe
    #Second argument: list of male-specific genes
    #Third argumnet: file where we store enrichment analysis
    apply_enrich(tmp,male_df.index.to_list(),outfile)

    ## enrichment analysis for female: similar as male

    outfile=current_folder+'/Female_Cluster_enrichment_final_new.txt'
    apply_enrich(tmp ,female_df.index.to_list(),outfile)

    ## enrichment analysis for both sexes: : similar as male

    outfile=current_folder+'/Both_Cluster_enrichment_with_phenotype_final_new.txt'
    apply_enrich(tmp,both_df.index.to_list(),outfile)
if __name__ == '__main__':

   enrich_cluster_final()
