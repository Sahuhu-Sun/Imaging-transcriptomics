
import pandas as pd
from scipy.stats import pearsonr

# Background gene
df1 = pd.read_excel(path+'expression_L.xlsx')

# PLS1 significant genes
df2 = pd.read_excel(path+'PLS_res.xlsx')

# AD-related genes from GeneCards
df3 = pd.read_csv(path+'GeneCards-SearchResults.csv')
df3 = df3[df3['Relevance score']>7]

bgge = df1.columns
opge0 = list( set(bgge)&set(list(df3['Gene Symbol'])) )
opge = list( set(opge0)&set(list(df2['Genes'])) )

print("-Number of genes overlapping with background genes：",len(opge0))
print('-Number of genes overlapping with PLS1 significant genes',len(opge))
print(df2[df2['Genes'].isin(opge)])

# t-map correlation analysis
df4 = pd.read_excel(path+'t_map_pval.xlsx')
t_stats = list(df4['t_stats'])
t_stats = t_stats[0:41]

exp_opge = df1[opge]
corr_res = pd.DataFrame()
ps=[];cors=[]
for ge in opge:
    corr, p_value = pearsonr(exp_opge[ge], t_stats)
    ps.append(p_value)
    cors.append(corr)

corr_res['Genes'] = opge
corr_res['Correlation'] = cors
corr_res['p-value'] = ps
corr_res = corr_res.sort_values(by='Correlation')
print(corr_res)

