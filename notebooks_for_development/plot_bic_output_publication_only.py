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

BIC_baseline = 458.1620541985783 # from Layden a,b,c,d model

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

#import ipdb; ipdb.set_trace()
# make list of strings for more elegant plotting
list_pre = list(df_sorted["coeffs_array"])

# flag the alphanumeric characters to extract them
char_string_list = list()
for combo_num in range(0,len(list_pre)):
    alpha_bool = [s.isalpha() for s in  list_pre[combo_num]]
    result = [x for x, y in zip(list_pre[combo_num], alpha_bool) if y]
    char_string = ''.join(result)
    char_string_list.append(char_string)

# replace the 'NaN' with 'Layden 94' for clarity
char_string_list_post = [p.replace("NaN","L94") for p in char_string_list]

ax = df_sorted["del_bic"].plot(kind='barh')
ax.set_ylabel("Coefficient permutation", fontsize=14)
ax.set_xlabel("$\Delta$BIC", fontsize=14)
ax.set_yticklabels(char_string_list_post, fontsize=7)
ax.set_xlim([-100,300])

# stagger the tick marks
for tick in ax.yaxis.get_major_ticks()[1::2]:
    tick.set_pad(30)

#for p in ax.patches:
#    ax.annotate(str(df_sorted["coeffs_array"][p]), (p.get_x() * 1.0, p.get_height() * 1.005))

plt.tight_layout()
plt.savefig(output_plot_name)
print("Saved BIC plot as " + output_plot_name)

# a second plot, to just look at models that fit better than Layden 1994
plt.clf()
ax = df_sorted["del_bic"][:22].plot(kind='barh')
ax.set_ylabel("Coefficient permutation", fontsize=14)
ax.set_xlabel("$\Delta$BIC", fontsize=14)
ax.set_yticklabels(char_string_list_post[:22], fontsize=12)
ax.set_xlim([-320,10])
# stagger the tick marks
for tick in ax.yaxis.get_major_ticks()[1::2]:
    tick.set_pad(30)
#for p in ax.patches:
#    ax.annotate(str(df_sorted["coeffs_array"][p]), (p.get_x() * 1.0, p.get_height() * 1.005))
plt.tight_layout()
plt.savefig("junk2.pdf", bbox_inches='tight')
print("Saved BIC plot as " + "junk2.pdf")
#plt.show()
