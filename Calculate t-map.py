import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

df = pd.read_excel(path + '/subjects.xlsx', index_col=0,sheet_name=1) 
df['Dx'] = df['Dx'].astype('category')

# Take out the column corresponding to the SUVR value of the brain region.
ya = df.iloc[:, :83] 

t_stats=[];p_val=[];age_p=[];sex_p=[];edu_p=[]
# Regression covariate
for i in range(83):
    x = df[['age','sex','edu','Dx']].copy()
    x.loc[:,'suvr'] = list(ya.loc[:, i])
    
    model = smf.ols('suvr ~ age + sex + edu + C(Dx, Treatment("CN"))', data=x).fit()
    
    # Extract t-statistic and p-value of group variables
    t_stats.append(model.tvalues.iloc[1])
    p_val.append(model.pvalues.iloc[1])
    age_p.append(model.pvalues.iloc[2])
    sex_p.append(model.pvalues.iloc[3])
    edu_p.append(model.pvalues.iloc[4])
    
    print(model.rsquared)
    print(model.summary())

# Save statistics and significance
df_new=pd.DataFrame()
df_new['t_stats']=t_stats
df_new['p_val']=p_val
df_new['age_p']=age_p
df_new['sex_p']=sex_p
df_new['edu_p']=edu_p

# Attach brain area information
dkatlas = pd.read_csv(path + '/atlas-desikankilliany.csv')
df_new['label']=dkatlas['label']
df_new['hemisphere']=dkatlas['hemisphere']

df_new.to_excel(path + '/t_map_pval.xlsx')
