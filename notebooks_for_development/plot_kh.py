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
isometal_balmer_abcissa = np.arange(2,16,0.2)
retrieved_K_isometal_neg3pt0 = expanded_layden_all_coeffs(coeff_array=coeff_array, H=isometal_balmer_abcissa, F=-3.0)
retrieved_K_isometal_neg2pt5 = expanded_layden_all_coeffs(coeff_array=coeff_array, H=isometal_balmer_abcissa, F=-2.5)
retrieved_K_isometal_neg2pt0 = expanded_layden_all_coeffs(coeff_array=coeff_array, H=isometal_balmer_abcissa, F=-2.)
retrieved_K_isometal_neg1pt5 = expanded_layden_all_coeffs(coeff_array=coeff_array, H=isometal_balmer_abcissa, F=-1.5)
retrieved_K_isometal_neg1pt0 = expanded_layden_all_coeffs(coeff_array=coeff_array, H=isometal_balmer_abcissa, F=-1.)
retrieved_K_isometal_neg0pt5 = expanded_layden_all_coeffs(coeff_array=coeff_array, H=isometal_balmer_abcissa, F=-0.5)
retrieved_K_isometal_pos0pt0 = expanded_layden_all_coeffs(coeff_array=coeff_array, H=isometal_balmer_abcissa, F=0.0)
retrieved_K_isometal_pos0pt2 = expanded_layden_all_coeffs(coeff_array=coeff_array, H=isometal_balmer_abcissa, F=0.2)

# plot it
plt.clf()
plt.figure(figsize=(20,10))

# underplot isometallicity lines
color_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
plt.plot(isometal_balmer_abcissa, retrieved_K_isometal_pos0pt2, linestyle="-", label="Isometal, Fe/H=+0.2", linewidth=5, color=color_list[-2])
plt.plot(isometal_balmer_abcissa, retrieved_K_isometal_pos0pt0, linestyle="-", label="Isometal, Fe/H=+0.0", linewidth=5, color=color_list[-3])
plt.plot(isometal_balmer_abcissa, retrieved_K_isometal_neg0pt5, linestyle="-", label="Isometal, Fe/H=-0.5", linewidth=5, color=color_list[-4])
plt.plot(isometal_balmer_abcissa, retrieved_K_isometal_neg1pt0, linestyle="-", label="Isometal, Fe/H=-1.0", linewidth=5, color=color_list[-5])
plt.plot(isometal_balmer_abcissa, retrieved_K_isometal_neg1pt5, linestyle="-", label="Isometal, Fe/H=-1.5", linewidth=5, color=color_list[-6])
plt.plot(isometal_balmer_abcissa, retrieved_K_isometal_neg2pt0, linestyle="-", label="Isometal, Fe/H=-2.0", linewidth=5, color=color_list[-7])
plt.plot(isometal_balmer_abcissa, retrieved_K_isometal_neg2pt5, linestyle="-", label="Isometal, Fe/H=-2.5", linewidth=5, color=color_list[-8])
plt.plot(isometal_balmer_abcissa, retrieved_K_isometal_neg3pt0, linestyle="--", label="Isometal, Fe/H=-3.0", linewidth=5, color=color_list[-9])


# calculated predicted
retrieved_K = expanded_layden_all_coeffs(coeff_array, df_choice["EW_Balmer"], df_choice["feh"])

# data points
#print(len(df_choice["final_feh_center"]))
# non-redundant list of metallicities

###
plt.clf()
# kludge to make legend come out right
df_choice["[Fe/H]"] = df_choice["feh"]
df_choice["log(g)"] = df_choice["logg"]

#fig, ax = plt.subplots()
#ax.annotate('Fastest Car', xy=(5,5), xytext=(5,5), xycoords='data', fontsize=20)

cmap = "Reds"
sns.set(font_scale=1.5)
sns.set_style(style='white')
'''
sns.relplot(
    data=df_choice,
    x="EW_Balmer", y="EW_CaIIK",
    hue="[Fe/H]", size="log(g)", edgecolor="k",
    palette=cmap, sizes=(30,160))
'''
g = sns.relplot(
    data=df_choice,
    x="EW_Balmer", y="EW_CaIIK",
    hue="[Fe/H]", size="log(g)", edgecolor="k",
    palette=cmap, sizes=(30,160)
)

feh_nonred = np.array([-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.2])

#ax.text(isometal_balmer_abcissa[0],retrieved_K_isometal_pos0pt2[0],"hello",color="k")


sns.lineplot(x=isometal_balmer_abcissa, y=retrieved_K_isometal_pos0pt2)
sns.lineplot(x=isometal_balmer_abcissa, y=retrieved_K_isometal_pos0pt0)
sns.lineplot(x=isometal_balmer_abcissa, y=retrieved_K_isometal_neg0pt5)
sns.lineplot(x=isometal_balmer_abcissa, y=retrieved_K_isometal_neg1pt0)
sns.lineplot(x=isometal_balmer_abcissa, y=retrieved_K_isometal_neg1pt5)
sns.lineplot(x=isometal_balmer_abcissa, y=retrieved_K_isometal_neg2pt0)
sns.lineplot(x=isometal_balmer_abcissa, y=retrieved_K_isometal_neg2pt5)

#sns.lineplot(x=isometal_balmer_abcissa, y=retrieved_K_isometal_neg3pt0, style=True, dashes=[(2,2)])



#handles, labels = ax.get_legend_handles_labels()

#g.legend_.set_title(None)

'''
g = sns.relplot(
    data=df_choice,
    x="EW_Balmer", y="EW_CaIIK",
    hue="[Fe/H]", size="log(g)", edgecolor="k",
    palette=cmap, sizes=(30,160)
)
'''
'''
f = sns.lineplot(
    data=df_choice,
    x="EW_Balmer", y="EW_CaIIK"
)
'''


'''
for num_met in range(0,len(feh_nonred)):
    cond_metal = (df_choice["feh"] == feh_nonred[num_met])
    df_subset = df_choice.where(cond_metal)
    sns.lineplot(x=isometal_balmer_abcissa, y=retrieved_K_isometal_pos0pt2)

    plt.plot(isometal_balmer_abcissa, retrieved_K_isometal_pos0pt2, linestyle="-", label="Isometal, Fe/H=+0.2", linewidth=5, color=color_list[-2])
    plt.plot(isometal_balmer_abcissa, retrieved_K_isometal_pos0pt0, linestyle="-", label="Isometal, Fe/H=+0.0", linewidth=5, color=color_list[-3])
    plt.plot(isometal_balmer_abcissa, retrieved_K_isometal_neg0pt5, linestyle="-", label="Isometal, Fe/H=-0.5", linewidth=5, color=color_list[-4])
    plt.plot(isometal_balmer_abcissa, retrieved_K_isometal_neg1pt0, linestyle="-", label="Isometal, Fe/H=-1.0", linewidth=5, color=color_list[-5])
    plt.plot(isometal_balmer_abcissa, retrieved_K_isometal_neg1pt5, linestyle="-", label="Isometal, Fe/H=-1.5", linewidth=5, color=color_list[-6])
    plt.plot(isometal_balmer_abcissa, retrieved_K_isometal_neg2pt0, linestyle="-", label="Isometal, Fe/H=-2.0", linewidth=5, color=color_list[-7])
    plt.plot(isometal_balmer_abcissa, retrieved_K_isometal_neg2pt5, linestyle="-", label="Isometal, Fe/H=-2.5", linewidth=5, color=color_list[-8])
    plt.plot(isometal_balmer_abcissa, retrieved_K_isometal_neg3pt0, linestyle="--", label="Isometal, Fe/H=-3.0", linewidth=5, color=color_list[-9])

    #plt.errorbar(df_choice["EW_Balmer"].where(cond_metal), df_choice["EW_CaIIK"].where(cond_metal),
    #    xerr=df_choice["err_EW_Balmer_based_Robo"].where(cond_metal), yerr=df_choice["err_EW_CaIIK_from_robo"].where(cond_metal),
    #    marker="o", color=color_list[-9+num_met], markersize=12, ecolor="gray", ls="none", label=str(feh_nonred[num_met]))
'''

g.fig.set_size_inches(28,8)
#sns.set(rc={'figure.figsize':(40,8.27)}) # figure size
#sns.set(font_scale = 25)

g.set_ylabels(r"$W_{K}$"+" ($\AA$)", fontsize=20)
g.set_xlabels(r"$W_{B}$"+" ($\AA$)", fontsize=20)
g.set(ylim=(0, 25))

#g.set(xlim=(0, 20))
'''
ax.set_ylabel(r"$W_{K}$"+" ($\AA$)", fontsize=20)
ax.set_xlabel(r"$W_{B}$"+" ($\AA$)", fontsize=20)

plt.legend(loc='upper right')
plt.savefig("junk.pdf")
'''
#g.set_xticklabels([str(i) for i in sns.get_xticks()], fontsize = 20)
#import ipdb; ipdb.set_trace()
# replace labels
#new_labels = ['[Fe/H]', 'log(g)']
#for t, l in zip(g._legend.texts, new_labels):
#    t.set_text(l)
#plt.legend(labels=['[Fe/H]', 'log(g)'])
g.savefig(file_name_write, bbox_inches='tight')
###
'''
feh_nonred = np.array([-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.2])
for num_met in range(0,len(feh_nonred)):
    cond_metal = (df_choice["feh"] == feh_nonred[num_met])
    plt.errorbar(df_choice["EW_Balmer"].where(cond_metal), df_choice["EW_CaIIK"].where(cond_metal),
        xerr=df_choice["err_EW_Balmer_based_Robo"].where(cond_metal), yerr=df_choice["err_EW_CaIIK_from_robo"].where(cond_metal),
        marker="o", color=color_list[-9+num_met], markersize=12, ecolor="gray", ls="none", label=str(feh_nonred[num_met]))
#plt.scatter(df_choice["EW_Balmer"], retrieved_K,
#    label="Retrieved, Modified Layden eqn")
# connect the empirical-retrieved dots, using list comprehension
#[plt.plot([df_choice["EW_Balmer"][j],df_choice["EW_Balmer"][j]],
#          [df_choice["EW_CaIIK"][j],retrieved_K[j]], color="k") for j in range(len(df_choice["feh"]))]

plt.ylabel("CaIIK "+r"$W_{K}$"+" ($\AA$)", fontsize=25)
plt.xlabel("Balmer "+r"$W_{B}$"+" ($\AA$)", fontsize=25)
#plt.title(str(new_coeffs_array[t]) + "\nBIC = " + str(bic))
plt.legend(fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylim([0,25])
plt.savefig(file_name_write)
plt.close()

print("Wrote",file_name_write)
print("----------")
'''


# baseline check
'''
abcd_older = [12.513685,-0.78716521,3.8778512,-0.24297523]
abcd_now = [19.09,-1.48,6.28,-0.49] # get more precise values later

K_baseline = original_layden_abcd(coeff_array=abcd_now,H=df_choice["EW_Balmer"],F=df_choice["feh"])
'''
