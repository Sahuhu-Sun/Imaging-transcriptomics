import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

from sklearn.cross_decomposition import PLSRegression
from statsmodels.stats.multitest import multipletests

from scipy import stats
from scipy.stats import pearsonr, norm

from sklearn.utils import resample

df_t = pd.read_excel(path + '/t_map_pval.xlsx')
expression_L = pd.read_excel(path + '/expression_L.xlsx')
Gene = expression_L.columns

t_stats = list(df_t['t_stats'])
t_stats_L = t_stats[0:41]

# input
X = expression_L  
Y = np.array(t_stats_L)  

genes = np.array(Gene)
geneindex = np.arange(len(genes))

# Standardize x and y
scaler_X = StandardScaler()
X_scaled = scaler_X.fit_transform(X)

scaler_Y = StandardScaler()
Y_scaled = scaler_Y.fit_transform(Y.reshape(-1, 1)).ravel()

# PLS 
# choose two principal components here after several attempts
pls = PLSRegression(n_components=2)  
pls.fit(X_scaled, Y_scaled)

# Extracting the weights of PLS1 and PLS2
pls_weights1 = pls.x_weights_[:, 0]  # PLS1
pls_weights2 = pls.x_weights_[:, 1]  # PLS2

# Sorting gene sequence in advance
PLS1_order = np.argsort(-pls_weights1)  
PLS2_order = np.argsort(-pls_weights2)

PLS1w = pls_weights1[PLS1_order]
PLS1ids = genes[PLS1_order]
geneindex1 = geneindex[PLS1_order]

PLS2w = pls_weights2[PLS2_order]
PLS2ids = genes[PLS2_order]
geneindex2 = geneindex[PLS2_order]

# Bootstrapping
n_iterations = 1000
PLS1weights = np.zeros((n_iterations, X_scaled.shape[1]))
PLS2weights = np.zeros((n_iterations, X_scaled.shape[1]))

for i in range(n_iterations):
    X_resampled, Y_resampled = resample(X_scaled, Y_scaled, random_state=i)
    pls.fit(X_resampled, Y_resampled)
    
    temp1 = pls.x_weights_[:, 0]
    temp2 = pls.x_weights_[:, 1]
    
    # Sort in original order
    newW1 = temp1[PLS1_order]
    newW2 = temp2[PLS2_order]
    
    # Make sure the direction is the same
    if np.corrcoef(PLS1w, newW1)[0, 1] < 0:
        newW1 *= -1
    if np.corrcoef(PLS2w, newW2)[0, 1] < 0:
        newW2 *= -1

    PLS1weights[i, :] = newW1
    PLS2weights[i, :] = newW2

# Calculation standard error
PLS1sw = np.std(PLS1weights, axis=0)
PLS2sw = np.std(PLS2weights, axis=0)

# Calculate Bootstrap weight
PLS1z = PLS1w / PLS1sw
PLS2z = PLS2w / PLS2sw

p_values1 = 2 * (1 - norm.cdf(np.abs(PLS1z)))
p_values2 = 2 * (1 - norm.cdf(np.abs(PLS2z)))

df_pls1 = pd.DataFrame({'Genes': PLS1ids, 'Weights': PLS1w, 'p_value': p_values1,'Z_score': PLS1z})
df_pls2 = pd.DataFrame({'Genes': PLS2ids, 'Weights': PLS2w, 'p_value': p_values2,'Z_score': PLS2z})

# FDR 
df_pls1["FDR_BH"] = multipletests(df_pls1["p_value"], method="fdr_bh")[1]
df_pls2["FDR_BH"] = multipletests(df_pls2["p_value"], method="fdr_bh")[1]

significant_genes1 = df_pls1[df_pls1["FDR_BH"] <= 0.05]  

PLS1_pos = significant_genes1[significant_genes1["Weights"] > 0]
PLS1_neg = significant_genes1[significant_genes1["Weights"] < 0]

print(f"After FDR correction, the number of PLS1 genes significantly correlated with p<0.05: {len(significant_genes1)}")
