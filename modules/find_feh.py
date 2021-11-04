import glob
import pandas as pd
import os
import numpy as np
import sys
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits
from rrlyrae_metallicity.modules2 import *
from . import *

class findFeH():
    '''
    Sample the a, b, c, d posteriors as found via the MCMC, to find Fe/H given the equivalent widths
    of the science spectra
    '''

    def __init__(self,
                good_ew_info_file = config_apply["data_dirs"]["DIR_EW_PRODS"]+config_apply["file_names"]["RESTACKED_EW_DATA_W_NET_BALMER_ERRORS"],
                calib_read_file = config_apply["data_dirs"]["DIR_SRC"] + config_apply["file_names"]["CALIB_SOLN"],
                 write_pickle_dir = config_apply["data_dirs"]["DIR_PICKLE"],
                 verbose = False):

        self.ew_file = good_ew_info_file
        self.calib_file = calib_read_file
        self.write_pickle_dir = write_pickle_dir

        # the table with the EW data
        self.ew_data = pd.read_csv(self.ew_file).copy(deep=True)
        # the file containing the MCMC posterior chains of coefficients
        # (note this does not distinguish yet between 4 or 8 coefficients)
        self.mcmc_chain = Table.read(self.calib_file, hdu=1)

        # retrieve Teff vs. Balmer EW coefficients
        # (note file is being read in a second time; not elegant)
        hdul = fits.open(self.calib_file)
        self.soln_header = hdul[1].header

    def __call__(self):

        pass


    def feh_layden(self,coeff_a,coeff_b,coeff_c,coeff_d,H,K):
        '''
        Finds Fe/H given equivalent widths (in angstroms), from
        K = a + b*H + c*[Fe/H] + d*H*[Fe/H]  (Layden 1994 Eqn. 7)
        '''

        feh = np.divide(K-coeff_a-coeff_b*H,coeff_c+coeff_d*H)

        return feh

    def feh_abcdfghk(self,coeff_a,coeff_b,coeff_c,coeff_d,coeff_f,coeff_g,coeff_h,coeff_k,H,K):
        '''
        Finds Fe/H given equivalent widths (in angstroms), from
        K = a + b*H + c*[Fe/H] + d*H*[Fe/H] + f*(H^2) + g*([Fe/H]^2) + h*(H^2)*[Fe/H] + k*H*([Fe/H]^2)
        '''

        A_cap = coeff_g + coeff_k*H
        B_cap = coeff_c + coeff_d*H + coeff_h*np.power(H,2)
        C_cap = coeff_a + coeff_b*H + coeff_f*np.power(H,2) - K

        # since this involves a quadratic, there are two roots
        F_pos = np.divide(-B_cap + np.sqrt(np.power(B_cap,2.)-4*A_cap*C_cap),2*A_cap)
        F_neg = np.divide(-B_cap - np.sqrt(np.power(B_cap,2.)-4*A_cap*C_cap),2*A_cap)

        return F_pos, F_neg


    def pickle_feh_retrieval(self, write_out_filename):
        '''
        Find a Fe/H value for a combination of coefficients
        from the MCMC chain, and sample from the Balmer and
        CaIIK EWs, given their errors

        INPUTS:
        write_out_filename: the file name of everything, incl. retrieved Teff and Fe/H

        OUTPUTS:
        (csv is written to disk)
        final_table: dataframe equivalent of the written csv file, for unit testing
        '''

        ## ## find/input EWs for a single spectrum here; use stand-in EWs for the moment
        # number of samples to take within the Gaussian errors around Balmer, CaIIK EWs
        N_EW_samples = 1

        # loop over samples in the MCMC chain
        N_MCMC_samples = len(self.mcmc_chain)

        # check if there is already something else in pickle directory
        preexisting_file_list = glob.glob(self.write_pickle_dir + "/*.{*}")
        if (len(preexisting_file_list) != 0):
            logging.info("------------------------------")
            logging.info("Directory to pickle Fe/H retrievals to is not empty!!")
            logging.info(self.write_pickle_dir)
            logging.info("------------------------------")
            input("Do what you want with those files, then hit [Enter]")

        # add columns to data table to include retrieved Fe/H values
        self.ew_data["feh_retrieved"] = np.nan
        self.ew_data["err_feh_retrieved"] = np.nan
        self.ew_data["teff_retrieved"] = np.nan

        # loop over the rows of the table of good EW data, with each row
        # corresponding to a spectrum
        for row_num in range(0,len(self.ew_data)):

            logging.info("Finding Fe/H for spectrum " + str(self.ew_data.iloc[row_num]["realization_spec_file_name"]))

            Balmer_EW = self.ew_data.iloc[row_num]["EW_Balmer"]
            CaIIK_EW = self.ew_data.iloc[row_num]["EW_CaIIK"]
            err_Balmer_EW = self.ew_data.iloc[row_num]["EW_Balmer"]
            err_CaIIK_EW = self.ew_data.iloc[row_num]["EW_CaIIK"]

            # initialize array
            feh_sample_array = np.nan*np.ones((N_MCMC_samples, N_EW_samples))
            for t in range(0,N_MCMC_samples):

                # loop over each sample within the Gaussian around the EW errors
                for integral_piece in range(0,N_EW_samples):

                    ## ## MIGHT BE OVERKILL: USE OF ROBO ERROR AS GAUSSIAN WIDTH LEADS TO RIDICULOUSLY WIDE DISTRIBUTIONS
                    # set the offset (note mu=0; this is a relative offset)
                    offset_H = 0 # np.random.normal(loc = 0.0, scale = err_Balmer_EW)
                    offset_K = 0 # np.random.normal(loc = 0.0, scale = err_CaIIK_EW)

                    # find one value of Fe/H given those samples in Balmer and CaIIK EWs
                    if (len(self.mcmc_chain.columns)==4):
                        feh_1sample = self.feh_layden(coeff_a = self.mcmc_chain["a"][t],
                                          coeff_b = self.mcmc_chain["b"][t],
                                          coeff_c = self.mcmc_chain["c"][t],
                                          coeff_d = self.mcmc_chain["d"][t],
                                          H = Balmer_EW + offset_H,
                                          K = CaIIK_EW + offset_K)
                    elif (len(self.mcmc_chain.columns)==8):
                        feh_1sample = self.feh_abcdfghk(coeff_a = self.mcmc_chain["a"][t],
                                          coeff_b = self.mcmc_chain["b"][t],
                                          coeff_c = self.mcmc_chain["c"][t],
                                          coeff_d = self.mcmc_chain["d"][t],
                                          coeff_f = self.mcmc_chain["f"][t],
                                          coeff_g = self.mcmc_chain["g"][t],
                                          coeff_h = self.mcmc_chain["h"][t],
                                          coeff_k = self.mcmc_chain["k"][t],
                                          H = Balmer_EW + offset_H,
                                          K = CaIIK_EW + offset_K)
                        feh_1sample = feh_1sample[0] # just want positive answer

                    feh_sample_array[t][integral_piece] = feh_1sample


                print("Spectrum number " + str(row_num) + " out of " + str(len(self.ew_data)))
                print("MCMC sample " + str(t) + " out of " + str(N_MCMC_samples))
                print("-----")

            # write the results (note this pickle file just corresponds to one spectrum)
            self.ew_data.at[row_num,"feh_retrieved"] = np.median(feh_sample_array)
            self.ew_data.at[row_num,"err_feh_retrieved"] = np.std(feh_sample_array)
            self.ew_data.at[row_num,"teff_retrieved"] = np.add(
                                                                np.multiply(self.ew_data.iloc[row_num]["EW_Balmer"],self.soln_header["SLOPE_M"]),
                                                                self.soln_header["YINT_B"]
                                                                )

        final_table = self.ew_data.copy()

        final_table.to_csv(write_out_filename, index=False)
        logging.info("Wrote out retrieved [Fe/H] and Teff to " + write_out_filename)

        return final_table


    def compare_feh_synthetic(self):
        '''
        Retrieves pickle files and plots injected and retrieved Fe/H
        '''

        # excavate all pickle files in the directory
        pickle_list = glob.glob(self.write_pickle_dir + "/*.p")

        # initialize data frame to hold values
        # cols:
        # 'inj_feh':        injected [Fe/H]
        # 'err_inj_feh'     plus/minus error in injected [Fe/H]
        # 'retr_med_feh'    retrieved [Fe/H]
        # 'lower_sig_feh'   1-sigma lower bound of [Fe/H]
        # 'upper_sig_feh'   1-sigma upper bound of [Fe/H]
        # 'logg'            injected logg
        # 'Teff'            injected effective temperature Teff

        df = pd.DataFrame(columns=["inj_feh", "err_inj_feh", "retr_med_feh",
                                    "lower_err_ret_feh", "upper_err_ret_feh", "logg", "teff", "pickle_file_name"]) #, index=range(len(pickle_list)))


        for file_name in pickle_list:

            # load each item in pickle file (maybe redundant, since it is one dictionary)
            with open(file_name, "rb") as f:
                data_all = pickle.load(f)

            # calculate errors (just stdev for now, probably want to improve this later)
            feh_retrieved = np.nanmedian(data_all["feh_sample_array"])
            err_feh_retrieved = np.nanstd(data_all["feh_sample_array"])

            values_to_add = {"inj_feh": data_all["injected_feh"],
                            "err_inj_feh": data_all["err_injected_feh"],
                            "logg": data_all["logg"],
                            "teff": data_all["Teff"],
                            "retr_med_feh": feh_retrieved,
                            "lower_err_ret_feh": err_feh_retrieved,
                            "upper_err_ret_feh": err_feh_retrieved,
                            "pickle_file_name": os.path.basename(file_name)}

            row_to_add = pd.Series(values_to_add, name="x")
            df = df.append(row_to_add)

        import ipdb; ipdb.set_trace()

        print(data_all)

        # plot retrieved and injected metallicities
        # matplotlib to show error bars

        fig, axes = plt.subplots(2, 1, figsize=(15, 24), sharex=True)
        fig.suptitle("Retrieval comparison, from MCMC file\n" + str(self.calib_file))

        # Fe/H difference
        df["feh_diff"] = np.subtract(df["retr_med_feh"],df["inj_feh"])

        # introduce scatter in x
        scatter_x = 0.1*np.random.rand(len(df["inj_feh"]))
        df["inj_feh_scatter"] = np.add(scatter_x,df["inj_feh"])

        cmap = sns.color_palette("YlOrBr", as_cmap=True)
        #cmap = sns.rocket_palette(rot=-.2, as_cmap=True)

        axes[0].plot([-2.5,0.5],[-2.5,0.5],linestyle="--",color="k",zorder=0)

        # underplot error bars
        axes[0].errorbar(x=df["inj_feh_scatter"],y=df["retr_med_feh"],xerr=df["err_inj_feh"],yerr=df["lower_err_ret_feh"],linestyle="",color="k",zorder=1)

        g_abs = sns.scatterplot(
            ax=axes[0],
            data=df,
            x="inj_feh_scatter", y="retr_med_feh",
            hue="teff", size="logg",
            edgecolor="k",
            palette=cmap, sizes=(50, 150),
            zorder=10
        )
        axes[0].set_ylabel("Retrieved: [Fe/H]$_{r}$")

        axes[1].plot([-2.5,0.5],[0,0],linestyle="--",color="k",zorder=0)
        g_diff = sns.scatterplot(
            ax=axes[1],
            data=df,
            x="inj_feh_scatter", y="feh_diff",
            hue="teff", size="logg",
            edgecolor="k",
            palette=cmap, sizes=(50, 150),
            legend=False,
            zorder=10
        )
        axes[1].set_ylabel("Residual: [Fe/H]$_{r}$-[Fe/H]$_{i}$")
        axes[1].set_xlabel("Injected: [Fe/H]$_{i}$")
        #axes[0].set_ylim([-3.,10.])
        #axes[1].set_ylim([-0.45,0.8])

        plt.savefig("/Users/bandari/Desktop/junk.pdf")

        #g.set(xscale="log", yscale="log")
        #g.ax.xaxis.grid(True, "minor", linewidth=.25)
        #g.ax.yaxis.grid(True, "minor", linewidth=.25)
        #g.despine(left=True, bottom=True)
