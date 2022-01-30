#!/usr/bin/env python
# coding: utf-8

# This reads in Balmer line widths and plots the scaling with respect to H-delta

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file_name = "/Users/bandari/Documents/git.repos/rrlfe/ew_products/all_data_input_mcmc_20220130_run_1.csv"

df = pd.read_csv(file_name)

plt.figure(figsize=(10,7))
plt.scatter(df["EW_Hdelta"], np.add(df["EW_Heps"],10), s=12, label="W"+r"$\epsilon$ + 10")
plt.scatter(df["EW_Hdelta"], np.add(df["EW_Hbeta"],5), s=12, label="W"+r"$\beta$ + 5")
plt.scatter(df["EW_Hdelta"], df["EW_Hgamma"], s=12, label="W"+r"$\gamma$")
plt.xlabel(r"$W_{\delta}$"+" "+"($\AA$)", fontsize=25)
plt.ylabel("$W$"+" "+"($\AA$)", fontsize=25)
plt.legend(fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tight_layout()
plt.savefig("junk.pdf")
