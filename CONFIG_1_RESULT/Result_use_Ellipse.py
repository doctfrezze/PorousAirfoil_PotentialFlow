import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from COMPUTATION.Richardson_extrapolation import extrapolate_matrices
from PLOT import plot_convergence, plot_extrapolated_vs_aoa

#%% Case 1
# Lecture du fichier CSV
df = pd.read_csv("resultats_simulations_Ellipse4.csv")

# Tri des valeurs uniques de AoA et NumPan
aoa_values = sorted(df['AoA'].unique())
numpan_values = sorted(df['NumPan'].unique())

# Initialisation des matrices
CL_matrix = np.zeros((len(aoa_values), len(numpan_values)))
CD_matrix = np.zeros((len(aoa_values), len(numpan_values)))

# Remplissage des matrices
for i, aoa in enumerate(aoa_values):
    for j, numpan in enumerate(numpan_values):
        subset = df[(df['AoA'] == aoa) & (df['NumPan'] == numpan)]
        if not subset.empty:
            CL_matrix[i, j] = subset['CL'].values[0]
            CD_matrix[i, j] = subset['CD'].values[0]

CL_matrix_extrapolated1, CD_matrix_extrapolated1,h = extrapolate_matrices(CL_matrix,CD_matrix,numpan_values)

# h_list = 1 / NumPan_values
h_list = [max(numpan_values) / n for n in numpan_values]

plot_convergence(h_list, CL_matrix, CD_matrix, CL_matrix_extrapolated1, CD_matrix_extrapolated1, aoa_values)


#%% Case 2
# Lecture du fichier CSV
df = pd.read_csv("resultats_simulations_Ellipse5.csv")

# Tri des valeurs uniques de AoA et NumPan
aoa_values = sorted(df['AoA'].unique())
numpan_values = sorted(df['NumPan'].unique())

# Initialisation des matrices
CL_matrix = np.zeros((len(aoa_values), len(numpan_values)))
CD_matrix = np.zeros((len(aoa_values), len(numpan_values)))

# Remplissage des matrices
for i, aoa in enumerate(aoa_values):
    for j, numpan in enumerate(numpan_values):
        subset = df[(df['AoA'] == aoa) & (df['NumPan'] == numpan)]
        if not subset.empty:
            CL_matrix[i, j] = subset['CL'].values[0]
            CD_matrix[i, j] = subset['CD'].values[0]

CL_matrix_extrapolated2, CD_matrix_extrapolated2,h = extrapolate_matrices(CL_matrix,CD_matrix,numpan_values)


# h_list = 1 / NumPan_values
h_list = [max(numpan_values) / n for n in numpan_values]

plot_convergence(h_list, CL_matrix, CD_matrix, CL_matrix_extrapolated2, CD_matrix_extrapolated2, aoa_values)

print('Case 3 : ')


#%% Case 3
# Lecture du fichier CSV
df = pd.read_csv("resultats_simulations_Ellipse6.csv")

# Tri des valeurs uniques de AoA et NumPan
aoa_values = sorted(df['AoA'].unique())
numpan_values = sorted(df['NumPan'].unique())

# Initialisation des matrices
CL_matrix = np.zeros((len(aoa_values), len(numpan_values)))
CD_matrix = np.zeros((len(aoa_values), len(numpan_values)))

# Remplissage des matrices
for i, aoa in enumerate(aoa_values):
    for j, numpan in enumerate(numpan_values):
        subset = df[(df['AoA'] == aoa) & (df['NumPan'] == numpan)]
        if not subset.empty:
            CL_matrix[i, j] = subset['CL'].values[0]
            CD_matrix[i, j] = subset['CD'].values[0]

CL_matrix_extrapolated3, CD_matrix_extrapolated3,h = extrapolate_matrices(CL_matrix,CD_matrix,numpan_values)


# h_list = 1 / NumPan_values
h_list = [max(numpan_values) / n for n in numpan_values]

plot_convergence(h_list, CL_matrix, CD_matrix, CL_matrix_extrapolated3, CD_matrix_extrapolated3, aoa_values)








plot_extrapolated_vs_aoa(aoa_values, CL_matrix_extrapolated1, CD_matrix_extrapolated1, CL_matrix_extrapolated2, CD_matrix_extrapolated2, CL_matrix_extrapolated3, CD_matrix_extrapolated3)