import numpy as np
import math as math
import matplotlib.pyplot as plt
import pandas as pd

def CSV_SAVE(CL_lists,CD_lists,file_name,AoA_sweep,h):                           

    # Initialisation du tableau final
    data = []

    # Parcours de tous les cas
    for case_num, (CL_mat, CD_mat) in enumerate(zip(CL_lists, CD_lists), start=1):
        for i, aoa in enumerate(AoA_sweep):
            for j, h_j in enumerate(h):
                cl = CL_mat[i, j]
                cd = CD_mat[i, j]
                data.append([case_num, aoa, h_j, cl, cd])

    # Création du DataFrame
    df = pd.DataFrame(data, columns=['Case', 'AoA', 'h', 'CL', 'CD'])

    # Sauvegarde en CSV
    df.to_csv(file_name+".csv", index=False)