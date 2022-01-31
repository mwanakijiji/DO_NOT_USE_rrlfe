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
coeff_array= np.array([25.62681273081429,-2.7529892780145238,12.51076634183292,-1.128153778868404,
                    0.08955467256229102,1.4456813531059578,0.02987747113827764,-0.07646288477131769,0,0])
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
isometal_balmer_abcissa = np.arange(2,16,0.2)
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

# plot points to show median error bars
# define the line along which they will lie
balmer_dummy = np.arange(7,16,1)
contour_4_points = np.add(8.,expanded_layden_all_coeffs(coeff_array=coeff_array, H=balmer_dummy, F=0.2))
bool_p02 = (df_choice["feh"] == 0.2)
#import ipdb; ipdb.set_trace()
plt.errorbar(balmer_dummy[0],[20],
                xerr=np.median(df_choice["err_EW_Balmer_based_Robo"][bool_p02]),
                yerr=np.median(df_choice["err_EW_CaIIK_from_robo"][bool_p02]), ecolor="k") # [Fe/H] = +0.2
plt.text(balmer_dummy[0]+0.1,20+0.5,"[Fe/H] = +0.2")
'''
bool_p00 = (df_choice["feh"] == 0.0)
plt.errorbar(balmer_dummy[1],contour_4_points[1],
                xerr=np.median(df_choice["err_EW_Balmer_based_Robo"][bool_p00]),
                yerr=np.median(df_choice["err_EW_CaIIK_from_robo"][bool_p00]), ecolor="k")  # [Fe/H] = +0.0

bool_m05 = (df_choice["feh"] == -0.5)
plt.errorbar(balmer_dummy[2],contour_4_points[2],
                xerr=np.median(df_choice["err_EW_Balmer_based_Robo"][bool_m05]),
                yerr=np.median(df_choice["err_EW_CaIIK_from_robo"][bool_m05]), ecolor="k") # [Fe/H] = -0.5
'''
bool_m10 = (df_choice["feh"] == -1.0)
plt.errorbar(balmer_dummy[3],[15],
                xerr=np.median(df_choice["err_EW_Balmer_based_Robo"][bool_m10]),
                yerr=np.median(df_choice["err_EW_CaIIK_from_robo"][bool_m10]), ecolor="k") # [Fe/H] = -1.0
plt.text(balmer_dummy[3]+0.1,15+0.5,"-1.0")
'''
bool_m15 = (df_choice["feh"] == -1.5)
plt.errorbar(balmer_dummy[4],contour_4_points[4],
                xerr=np.median(df_choice["err_EW_Balmer_based_Robo"][bool_m15]),
                yerr=np.median(df_choice["err_EW_CaIIK_from_robo"][bool_m15]), ecolor="k")  # [Fe/H] = -1.5

bool_m20 = (df_choice["feh"] == -2.0)
plt.errorbar(balmer_dummy[5],contour_4_points[5],
                xerr=np.median(df_choice["err_EW_Balmer_based_Robo"][bool_m20]),
                yerr=np.median(df_choice["err_EW_CaIIK_from_robo"][bool_m20]), ecolor="k")  # [Fe/H] = -2.0
'''
bool_m25 = (df_choice["feh"] == -2.5)
plt.errorbar(balmer_dummy[6],[10],
                xerr=np.median(df_choice["err_EW_Balmer_based_Robo"][bool_m25]),
                yerr=np.median(df_choice["err_EW_CaIIK_from_robo"][bool_m25]), ecolor="k")  # [Fe/H] = -2.5
plt.text(balmer_dummy[6]+0.1,10+0.5,"-2.5")

# plot isometallicity contours
feh_nonred = np.array([-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.2])
plt.text(isometal_balmer_abcissa[0]-0.65,retrieved_K_isometal_pos0pt2[0]+2,"[Fe/H]")
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

plt.xlim([1,16])
plt.legend(loc="upper right", ncol=3)

g.fig.set_size_inches(28,8)

g.set_ylabels(r"$W_{K}$"+" ($\AA$)", fontsize=22)
g.set_xlabels(r"$W_{B}$"+" ($\AA$)", fontsize=22)
g.set(ylim=(0, 25))

g.savefig(file_name_write, bbox_inches='tight')
