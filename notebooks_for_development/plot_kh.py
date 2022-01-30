# This makes a KH plot with isometallicity contours

# Created from parent 2022 Jan. 30 by E.S.


# #### In the following, we plot fits in KH space and write out data including the BIC to select
# #### the best model among variations that consider up to third-degree terms involving H (Balmer
# #### EW) and F (Fe/H), viz.
#
# #### $K = a + bH + cF + dHF + f(H^{2}) + g(F^{2}) + h(H^{2})F + kH(F^{2}) $
# #### $+ m(H^{3}) + n(F^{3}) $


import pandas as pd
import numpy as np
import astropy
import itertools
import multiprocessing
import random
import string
import seaborn as sns
from astropy import stats
from scipy import optimize
import matplotlib
matplotlib.use('Agg') # necessary in the cloud
import matplotlib.pyplot as plt

# set the coefficients of the model
coeff_array= np.array([2.44651703e+01,-2.61158619e+00 ,1.23783649e+01, -1.09623588e+00,
                    8.67772095e-02,1.55202172e+00,2.90194894e-02,-8.25723508e-02,0,0])
#coeff_array = np.array([1.,1.,1.,1.,1.,1.,1.,1.,0,0])

# read in the data points to be overplotted
data_points_read = "/Users/bandari/Documents/git.repos/rrlfe/ew_products/all_data_input_mcmc_20220130_run_1.csv"
df_choice = pd.read_csv(data_points_read)

# set the file name to write
file_name_write = "junk.pdf"

def expanded_layden_all_coeffs(coeff_array,H,F):

    # definition of coefficients as of 2020 Mar 9:
    # K = a + bH + cF + dHF + f(H^{2}) + g(F^{2}) + hF(H^{2}) + kH(F^{2}) + m(H^{3}) + n(F^{3})

    a_coeff = coeff_array[0]
    b_coeff = coeff_array[1]
    c_coeff = coeff_array[2]
    d_coeff = coeff_array[3]
    f_coeff = coeff_array[4]
    g_coeff = coeff_array[5]
    h_coeff = coeff_array[6]
    k_coeff = coeff_array[7]
    m_coeff = coeff_array[8]
    n_coeff = coeff_array[9]

    K_calc = a_coeff + b_coeff*H + c_coeff*F + d_coeff*H*F + \
        f_coeff*np.power(H,2.) + g_coeff*np.power(F,2.) + \
        h_coeff*F*np.power(H,2.) + k_coeff*H*np.power(F,2.) + \
        m_coeff*np.power(H,3.) + n_coeff*np.power(F,3.)

    return K_calc


# make some isometallicity lines for the plot
isometal_balmer_abcissa = np.arange(2,14,0.2)
retrieved_K_isometal_neg3pt0 = expanded_layden_all_coeffs(coeff_array=coeff_array, H=isometal_balmer_abcissa, F=-3.0)
retrieved_K_isometal_neg2pt5 = expanded_layden_all_coeffs(coeff_array=coeff_array, H=isometal_balmer_abcissa, F=-2.5)
retrieved_K_isometal_neg2pt0 = expanded_layden_all_coeffs(coeff_array=coeff_array, H=isometal_balmer_abcissa, F=-2.)
retrieved_K_isometal_neg1pt5 = expanded_layden_all_coeffs(coeff_array=coeff_array, H=isometal_balmer_abcissa, F=-1.5)
retrieved_K_isometal_neg1pt0 = expanded_layden_all_coeffs(coeff_array=coeff_array, H=isometal_balmer_abcissa, F=-1.)
retrieved_K_isometal_neg0pt5 = expanded_layden_all_coeffs(coeff_array=coeff_array, H=isometal_balmer_abcissa, F=-0.5)
retrieved_K_isometal_pos0pt0 = expanded_layden_all_coeffs(coeff_array=coeff_array, H=isometal_balmer_abcissa, F=0.0)
retrieved_K_isometal_pos0pt2 = expanded_layden_all_coeffs(coeff_array=coeff_array, H=isometal_balmer_abcissa, F=0.2)


plt.clf()
# kludge to make legend come out right
df_choice["[Fe/H]"] = df_choice["feh"]
df_choice["log(g)"] = df_choice["logg"]

cmap = "Reds"
sns.set(font_scale=1.9)
sns.set_style(style='white')

g = sns.relplot(
    data=df_choice,
    x="EW_Balmer", y="EW_CaIIK",
    hue="[Fe/H]", size="log(g)", edgecolor="k", legend="full",
    palette=cmap, sizes=(30,160)
)

feh_nonred = np.array([-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.2])

p02 = sns.lineplot(x=isometal_balmer_abcissa, y=retrieved_K_isometal_pos0pt2, color="k", alpha=0.5)
plt.text(isometal_balmer_abcissa[0]-0.6,retrieved_K_isometal_pos0pt2[0],"+0.2")
p00 = sns.lineplot(x=isometal_balmer_abcissa, y=retrieved_K_isometal_pos0pt0, color="k", alpha=0.5)
plt.text(isometal_balmer_abcissa[0]-0.6,retrieved_K_isometal_pos0pt0[0],"+0.0")
m05 = sns.lineplot(x=isometal_balmer_abcissa, y=retrieved_K_isometal_neg0pt5, color="k", alpha=0.5)
plt.text(isometal_balmer_abcissa[0]-0.6,retrieved_K_isometal_neg0pt5[0],"-0.5")
m10 = sns.lineplot(x=isometal_balmer_abcissa, y=retrieved_K_isometal_neg1pt0, color="k", alpha=0.5)
plt.text(isometal_balmer_abcissa[0]-0.6,retrieved_K_isometal_neg1pt0[0],"-1.0")
m15 = sns.lineplot(x=isometal_balmer_abcissa, y=retrieved_K_isometal_neg1pt5, color="k", alpha=0.5)
plt.text(isometal_balmer_abcissa[0]-0.6,retrieved_K_isometal_neg1pt5[0],"-1.5")
m20 = sns.lineplot(x=isometal_balmer_abcissa, y=retrieved_K_isometal_neg2pt0, color="k", alpha=0.5)
plt.text(isometal_balmer_abcissa[0]-0.6,retrieved_K_isometal_neg2pt0[0],"-2.0")
m25 = sns.lineplot(x=isometal_balmer_abcissa, y=retrieved_K_isometal_neg2pt5, color="k", alpha=0.5)
plt.text(isometal_balmer_abcissa[0]-0.6,retrieved_K_isometal_neg2pt5[0],"-2.5")

plt.xlim([1,14])
plt.legend(loc="upper right", ncol=3)

g.fig.set_size_inches(28,8)

g.set_ylabels(r"$W_{K}$"+" ($\AA$)", fontsize=22)
g.set_xlabels(r"$W_{B}$"+" ($\AA$)", fontsize=22)
g.set(ylim=(0, 25))

g.savefig(file_name_write, bbox_inches='tight')
