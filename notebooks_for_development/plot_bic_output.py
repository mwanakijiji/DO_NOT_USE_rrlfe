#!/usr/bin/env python
# coding: utf-8

# This plots BIC values for different permutations of coefficients

# Created 2020 Mar 10 by E.S.

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import os
from dataclasses import dataclass
from typing import *



input_file_name = "bic_run_20220129_3_weighted_squares/results_bic_run_20220129_3_weighted_squares.csv"
output_plot_name = "junk.pdf"

df = pd.read_csv(input_file_name,names=["coeffs_array",
                                                    "bic",
                                                    "n_params",
                                                    "chi_squared",
                                                    "n_samples",
                                                    "a_coeff",
                                                    "b_coeff",
                                                    "c_coeff",
                                                    "d_coeff",
                                                    "f_coeff",
                                                    "g_coeff",
                                                    "h_coeff",
                                                    "k_coeff",
                                                    "m_coeff",
                                                    "n_coeff",
                                                    "res"],delimiter=";")


#BIC_baseline_logg2pt5only = 92.5231102439195
#BIC_baseline_logg3pt0only = 85.5973111789901
#BIC_baseline_logg2pt5and3pt0 = 168.61055978365414

BIC_baseline = 251.2375827551901 # from rough coefficients {a,b,c,d} at end of Feb 2021

# sort by BIC
df_sorted = df.sort_values("bic").reset_index()

# add column of del_BIC = BIC_i - BIC_baseline
## ## CHANGE LATER!
df_sorted["del_bic"] = np.subtract(df_sorted["bic"],BIC_baseline)


#plt.bar(np.arange(len(df_sorted["del_bic"])),df_sorted["del_bic"])
#plt.title("del_BIC")
#plt.xlabel("permutation number, sorted by BIC")
#plt.ylabel("BIC")
#plt.show()

ax = df_sorted["del_bic"].plot(kind='barh')
ax.set_title("$\Delta$BIC")
ax.set_ylabel("permutation, sorted by $\Delta$BIC")
ax.set_xlabel("BIC")
ax.set_yticklabels(list(zip(df_sorted["coeffs_array"],df_sorted["res"])), fontsize=3)
#ax.get_yaxis().set_visible(False)
ax.set_xlim([-200,600])

#for p in ax.patches:
#    ax.annotate(str(df_sorted["coeffs_array"][p]), (p.get_x() * 1.0, p.get_height() * 1.005))

plt.savefig(output_plot_name)
print("Saved BIC plot as " + output_plot_name)
#plt.show()
