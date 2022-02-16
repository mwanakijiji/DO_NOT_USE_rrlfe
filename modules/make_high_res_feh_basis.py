'''
Reads in literature metallicities and makes new Fe/H basis
'''

import pickle
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from . import *

class LitMetallicities():
    '''
    Class to
    1.   read in Fe/H values from the literature
    2.   initialize data set cross-referencing functionality
    '''

    def __init__(self,
                 source_dir=config_red["data_dirs"]["DIR_LIT_HIGH_RES_FEH"]):

        # stand-in that consists of our program star names
        self.our_program_stars = pd.read_csv(source_dir + "our_program_stars_names_only.csv")

        # Fe/H from Layden+ 1994; this may serve as the common basis for RRabs
        self.layden_feh = pd.read_csv(source_dir + "layden_1994_abundances.dat")
        # RES: "rather low"

        # Fe/H Clementini+ 1995
        self.clementini_feh = pd.read_csv(source_dir + "clementini_1995_abundances.dat")

        # Fe/H Fernley+ 1996
        self.fernley96_feh = pd.read_csv(source_dir + "fernley_1996_abundances.dat")
        # RES: 60,000, FeI & FeII, 5900-8100 A

        # Fe/H from Fernley+ 1997
        self.fernley97_feh = pd.read_csv(source_dir + "fernley_1997_abundances.dat")
        # RES: 60,000, two FeII lines, 5900-8100 A

        # log(eps) from Lambert+ 1996
        self.lambert_logeps = pd.read_csv(source_dir + "lambert_1996_abundances.dat")
        # RES: ~23,000, FeII + photometric models, 3600-9000 A

        # Fe/H from Wallerstein and Huang 2010, arXiv 1004.2017
        self.wallerstein_feh = pd.read_csv(source_dir + "wallerstein_huang_2010_abundances.dat")
        # RES: ~30,000, FeII

        # Fe/H from Chadid+ 2017 ApJ 835.2:187 (FeI and II lines)
        self.chadid_feh = pd.read_csv(source_dir + "chadid_2017_abundances.dat")
        # RES: 38000, FeI & FeII, 3400-9900 A

        # Fe/H from Liu+ 2013 Res Ast Astroph 13:1307
        self.liu_feh = pd.read_csv(source_dir + "liu_2013_abundances.dat")
        # RES: ~60,000, FeI (& FeII?), 5100-6400 A

        # Fe/H from Nemec+ 2013
        self.nemec_feh = pd.read_csv(source_dir + "nemec_2013_abundances.dat")
        # RES: ~65,000 or 36,000, FeI & FeII, 5150-5200 A

        # Fe/H from Solano+ 1997
        self.solano_feh = pd.read_csv(source_dir + "solano_1997_abundances.dat")
        # RES: 22,000 & 19,000, strong FeI lines, 4160-4390 & 4070-4490 A

        # Fe/H from Pancino+ 2015 MNRAS 447:2404
        self.pancino_feh = pd.read_csv(source_dir + "pancino_2015_abundances.dat")
        # RES: >30,000, FeI (weighted average), 4000-8500 A

        # Fe/H from Sneden+ 2017
        self.sneden_feh = pd.read_csv(source_dir + "sneden_2017_abundances.dat")
        # RES: ~27,000 (at 5000 A), FeI & FeII, 3400-9000 A


        # Now modify some of the datasets on an individual basis

        # convert Lambert's values, which are in terms of log(eps)
        # FeH = log(epsFe) - log(epsFe,sol)
        #     = log(epsFe) - log(NFe,sol/NH,sol)
        #     = log(epsFe) - 7.51 # value of 7.51 from Anstee+ 1997, MNRAS
        self.lambert_logeps["feh"] = np.subtract(self.lambert_logeps["log_eps_fe_spec"], 7.51)

        # average the values in Chadid from FeI and FeII lines
        self.chadid_feh["feh"] = np.mean([self.chadid_feh["fehI"].values,
                                          self.chadid_feh["fehII"].values], axis=0)

        ## ## INCLUDE SINGLE DATA PT FROM KOLENBERG+ 2010? (SEE CHADID+ 2017, FIG. 7)

        # FYI: average Fe/H values in Liu+ 2013 which were taken at different phases
        # liu_feh.groupby(liu_feh["name"], axis=0, as_index=False).mean()

        # FYI: average Fe/H values in Sneden+ 1997 which were taken at different epochs
        # sneden_feh.groupby(sneden_feh["name"], axis=0, as_index=False).mean()

        # Fe/H from Kemper+ 1982; this might serve as the common basis for RRcs
        self.kemper_feh = pd.read_csv(source_dir + "kemper_1982_abundances.dat")

        # Fe/H from Govea+ 2014
        ## ## note: Govea+ has abundances for each phase value, and this
        ## ## includes NLTE phases; how to get single Fe/H?
        self.govea_feh = pd.read_csv(source_dir + "govea_2014_abundances.dat")

        # put all read literature metallicities into a dictionary for testing
        '''
        # this doesn't really work, putting dataframes into a dictionary
        self.collated_feh = {"our_program_stars": self.our_program_stars,
                         "layden_feh": self.layden_feh,
                         "clementini_feh": self.clementini_feh,
                         "lambert_logeps": self.lambert_logeps,
                         "wallerstein_feh": self.wallerstein_feh,
                         "chadid_feh": self.chadid_feh,
                         "liu_feh": self.liu_feh,
                         "nemec_feh": self.nemec_feh,
                         "fernley96_feh": self.fernley96_feh,
                         "fernley97_feh": self.fernley97_feh,
                         "solano_feh": self.solano_feh,
                         "pancino_feh": self.pancino_feh,
                         "sneden_feh": self.sneden_feh,
                         "kemper_feh": self.kemper_feh,
                         "govea_feh": self.govea_feh}
        '''

        #####################

        # initialize arrays: essential info
        empir_spec_name_array = []
        star_name_array = []
        H_data_array = []
        K_data_array = []
        err_H_data_array = []
        err_K_data_array = []

        # initialize arrays: other info
        Hbet_data_array = []
        err_Hbet_data_array = []
        Hgam_data_array = []
        err_Hgam_data_array = []
        rHgam_data_array = [] # rescaled Hgamma
        err_rHgam_data_array = []
        Hdel_data_array = []
        err_Hdel_data_array = []
        Heps_data_array = []
        err_Heps_data_array = []


    def __call__(self):

        # make a list of all unique EMPIRICAL spectrum names
        unique_spec_names = line_data.drop_duplicates(subset="empir_spec_name")["empir_spec_name"]


    def return_some_raw_data(self):
        '''
        Return some read-in data, for testing
        '''

        return self.clementini_feh, self.pancino_feh


    def matchmaker(self,
                   input_table,
                   basis_table,
                   basis_dataset_name,
                   highres_dataset_name):
        '''
        Find what stars are common to two input tables, and return array of FeHs from the first table

        INPUTS:
        input_table: table I'm interested in checking for overlapping stars
            (pandas dataframe with col ["name"]: star name; col ["feh"]: Fe/H)
        basis_table: table with the names for which I am looking for repeats in the other table
            (pandas dataframe with col ["name"]: star name; col ["feh"]: Fe/H)
        basis_dataset_name: string indicating the name of the basis dataset
        highres_dataset_name: string indicating the name of the highres dataset

        OUTPUTS:
        pandas dataframe with
        1. overlapping star names
        2. FeHs from the input_table
        3. FeHs from the basis_table
        4. residuals in FeH: FeH_input - FeH_basis
        5. string indicating the high-res dataset being matched
        '''

        self.input_table = input_table
        self.basis_table = basis_table

        input_feh = [] # Fe/H of high-res study
        basis_feh = [] # Fe/H of basis (ex. Layden 1994)
        star_name_array = [] # name of star

        # scan over each row in input table
        for row in range(0, len(input_table)):
            if (basis_table["name"] == input_table["name"][row]).any():
                input_feh = np.append(input_feh, input_table["feh"][row])
                basis_feh = np.append(basis_feh,
                                      basis_table.loc[basis_table["name"] == input_table["name"][row]]["feh"])
                star_name_array = np.append(star_name_array, input_table["name"][row])

        d = dict()
        d["name_star"] = star_name_array
        d["feh_highres"] = input_feh
        d["feh_basis"] = basis_feh
        d["name_highres_dataset"] = np.repeat(highres_dataset_name, len(star_name_array))
        d["name_basis_dataset"] = np.repeat(basis_dataset_name, len(star_name_array))

        df = pd.DataFrame(data=d)

        return df


    def match_highres_w_basis(self, star_type):
        '''
        Find what stars overlap with basis data set, and return star name, data set names,
            FeH values, residuals

        The functionality of LitMetallicities is inherited
        N.b. There are no offsets applied yet (as are applied in Chadid+ 2017 plots)

        INPUTS:
        input_table: table of high-res-derived Fe/H values, which I want to cross-ref with
            a basis (like Layden 94 or Kemper 82)
        basis_table: the table which serves as the basis set

        OUTPUTS:
        df: concatenated frames containing
            ["name_star"]: star_name_array
            ["feh_highres"]: Fe/H from high-res study
            ["feh_basis"]: Fe/H from basis set
            ["name_highres_dataset"]: string indicating the high-res dataset
            ["name_basis_dataset"]: string indicating the basis set
        '''

        # define the basis data set (Layden 1994 for RRabs, or Kemper+ 1982 for RRcs)
        if star_type == "RRab":
            type_string = "ab"
            basis_set = self.layden_feh
            basis_string = "Layden RRab basis set" # string for plots
            basis_name = "layden_1994"
        elif star_type == "RRc":
            type_string = "c"
            basis_set = self.kemper_feh
            basis_string = "Kemper RRc basis set"
            basis_name = "kemper_1982"
        else:
            sys.exit("Error! No RR Lyrae subtype chosen.")


        ## match ALL available high-res studies with the basis set
        ## ## are these all the high-res studies I want? are more? any fewer?
        ## ## also indicate what tables, etc. in the papers these numbers are from
        pd_Lambert_1996 = self.matchmaker(input_table=self.lambert_logeps,
                                          basis_table=basis_set,
                                          basis_dataset_name=basis_name,
                                          highres_dataset_name="lambert_1996")
        pd_Nemec_2013 = self.matchmaker(input_table=self.nemec_feh,
                                        basis_table=basis_set,
                                        basis_dataset_name=basis_name,
                                        highres_dataset_name="nemec_2013")
        pd_Pancino_2015 = self.matchmaker(input_table=self.pancino_feh,
                                        basis_table=basis_set,
                                        basis_dataset_name=basis_name,
                                        highres_dataset_name="pancino_2015")
        pd_Chadid_2017 = self.matchmaker(input_table=self.chadid_feh,
                                         basis_table=basis_set,
                                         basis_dataset_name=basis_name,
                                         highres_dataset_name="chadid_2017")
        pd_Fernley_1997 = self.matchmaker(input_table=self.fernley97_feh,
                                          basis_table=basis_set,
                                          basis_dataset_name=basis_name,
                                          highres_dataset_name="fernley_1997")
        pd_Solano_1997 = self.matchmaker(input_table=self.solano_feh,
                                         basis_table=basis_set,
                                         basis_dataset_name=basis_name,
                                         highres_dataset_name="solano_1997")
        pd_Wallerstein_2010 = self.matchmaker(input_table=self.wallerstein_feh,
                                              basis_table=basis_set,
                                              basis_dataset_name=basis_name,
                                              highres_dataset_name="wallerstein_2010")

        # for Liu+ 2013, we need to group multiple Fe/H values by star name
        # (the grouping is done here rather than further up because a bug causes
        # the grouped column to disappear)
        self.liu_feh_grouped = self.liu_feh.groupby(self.liu_feh["name"],
                                                    axis=0,
                                                    as_index=False).mean()
        pd_Liu_2013 = self.matchmaker(input_table=self.liu_feh_grouped,
                                      basis_table=basis_set,
                                      basis_dataset_name=basis_name,
                                      highres_dataset_name="liu_2013") # Liu+ 2013

        # for Govea+ 2014, we need to group multiple Fe/H_I and Fe/H_II values by star name,
        # and get single Fe/H_i values by averaging them for each star
        # (the grouping is done here rather than further upstream because otherwise a bug causes
        # the grouped column to disappear)
        ## ## note that the phases of these spectra are also averaged;
        ## ## maybe I should winnow them according to our own phase boundaries
        self.govea_feh_grouped = self.govea_feh.groupby(self.govea_feh["name"],
                                                        axis=0,
                                                        as_index=False).mean()
        # now, average the Fe/H_I and Fe/H_II values to get single Fe/H values
        self.govea_feh_grouped["feh"] = self.govea_feh_grouped[["feIh", "feIIh"]].mean(axis=1)
        pd_Govea_2014 = self.matchmaker(input_table=self.govea_feh_grouped,
                                        basis_table=basis_set,
                                        basis_dataset_name=basis_name,
                                        highres_dataset_name="govea_2014") # Govea+ 2014

        # merge dataframes
        pd_collected = [pd_Lambert_1996,
                        pd_Nemec_2013,
                        pd_Pancino_2015,
                        pd_Liu_2013,
                        pd_Chadid_2017,
                        pd_Fernley_1997,
                        pd_Solano_1997,
                        pd_Wallerstein_2010,
                        pd_Govea_2014]
        df = pd.concat(pd_collected).reset_index()

        return df


def plot_mapping(input_matches,
                 mapped,
                 title_string,
                 plot_file_name,
                 write_plot_subdir=config_red["data_dirs"]["DIR_FYI_INFO"],
                 write_plot=True):

    # if plot is not to be written
    if write_plot == False:
        return

    plt.clf()
    limits = [-3.0, 0.5] # limits of Fe/H to plot
    fig, axs = plt.subplots(1, 1, figsize=(10, 10))
    axs.plot([limits[0], limits[1]], [limits[0], limits[1]], linestyle="--") # 1-to-1 line
    axs.scatter(input_matches["feh_basis"], mapped) # input_matches vs. basis
    axs.set_xlim(limits[0], limits[1])
    axs.set_ylabel("Fe/H, high-res")
    axs.set_xlabel("Fe/H, basis")
    # plot the as-is Fe/H values as reported in all the high-res basis sets
    axs.scatter(input_matches["feh_basis"], input_matches["feh_highres"],
                        alpha=0.5,
                        label="feh_highres_all_preshift")
    # add star names
    for point in range(0, len(input_matches["name_star"])):
        axs.annotate(
                    input_matches["name_star"][point],
                    xy=(input_matches["feh_basis"][point],
                        mapped[point]),
                    xytext=(input_matches["feh_basis"][point]+0.1,
                        mapped[point]+0.06),
                    textcoords="data", ha="right", va="bottom",
                    fontsize=10,
                    arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=0"))
    fig.suptitle(title_string)
    plt.savefig(write_plot_subdir + plot_file_name, overwrite=True)


def return_offsets(data_postmatch,
                   chadid_offset=True):
    '''
    Fit linear regression to input high-res Fe/H and the basis set,
    and find offset for each dataset to match some constant

    INPUTS:
    data_postmatch: output from match_w_highres_basis(), containing
        Fe/H data from high-res studies and basis set

    OUTPUTS:
    df: dataframe containing
    ["name_highres_dataset"]: name indicating high-res study
    ["offset_highres_residuals"]: offset for the entire high-res
        dataset, which needs to be applied to the vector (feh_highres-feh_basis)
    '''

    # initialize dataframe
    col_names = ["name_highres", "offset_highres_residuals"]
    df_offsets = pd.DataFrame(columns=col_names)

    # loop over each high-res dataset
    highres_names = data_postmatch["name_highres_dataset"].unique() # names of high-res datasets
    index_counter = 0 # initialize this for appending dataframes
    for dataset_num in range(0, len(highres_names)):

        this_dataset_name = highres_names[dataset_num] # name of this dataset
        this_dataset = data_postmatch[data_postmatch["name_highres_dataset"].str.match(this_dataset_name)]

        # need at least N+1 data points
        N = 2
        if (len(this_dataset["feh_highres"]) > N):

            # find linear regression of residuals
            coeff = np.polyfit(this_dataset["feh_basis"],
                               np.subtract(this_dataset["feh_highres"],
                                           this_dataset["feh_basis"]), 1)
            print("Data to fit residuals to:")
            print(this_dataset)
            print("Lin reg fit to residuals:")
            print(coeff)
            limits = [-3.0, 0.5] # Fe/H limits to display
            line = np.multiply(coeff[0], limits)+coeff[1] # points to plot linear regression

            # find offset between residuals and Chadid+ 2017 at Fe/H=-1.25 (see their Fig. 6)
            if chadid_offset:
                chadid_y_125 = -0.10583621694962 # from Chadid line at Fe/H=-1.25
                feh_resid_to_peg = chadid_y_125 # peg the feh residuals to this y- point
                feh_basis_loc = -1.25 # corresponding x- value (Fe/H in the basis dataset)
            else:
                feh_resid_to_peg = 0. # peg the data at this y- point
                feh_basis_loc = 0. # location in the basis dataset
                print("Offset corresponds to (0, 0)!")

            # y-value of the unshifted linear regression line at Fe/H=-1.25
            this_y_125 = np.multiply(coeff[0], feh_basis_loc)+coeff[1]

            # offset to shift the linear regression line of the residuals to feh_resid_to_peg
            net_offset = chadid_y_125 - this_y_125

            # line_offset = np.add(line, net_offset)

            print("-----------------")
            print("Calculating offsets for dataset "+this_dataset_name)
            print("Y_offset to add to residuals in order to overlap with Chadid+ 2017 at Fe/H=-1.25:")
            print(net_offset)
            print("Number of overlapping stars:")
            print(len(this_dataset["feh_basis"]))

            # append name of dataset and its offset to array
            dict_this = {"name_highres":[this_dataset_name],
                         "offset_highres_residuals":[net_offset]}
            #print(pd.DataFrame(dict_this, index=[dataset_num]))

            # append info from this dataset
            #print(dict_this)
            #print(list(dict_this.items()))
            #print(pd.DataFrame(dict_this))
            df_offsets = df_offsets.append(pd.DataFrame(dict_this))
            index_counter += 1 # increase the counter

        elif (len(this_dataset["feh_highres"]) <= N):
            pass

    return df_offsets


def make_basis_via_offsets(df_to_offset,
                           df_offsets,
                           plot_string,
                           csv_string,
                           feh_basis_csv_dir=config_red["data_dirs"]["DIR_FYI_INFO"],
                           make_plot=True):
    '''
    Apply offsets (which may be from RRabs, RRcs, combo, etc.) to data to make a basis

    INPUTS:
    df_to_offset: dataframe containing Fe/H values which we want to find residuals for, offset, and map
        ["name_star"]: star name
        ["feh_highres"]: Fe/H from high-res study
        ["feh_basis"]: Fe/H from basis set
        ["name_highres_dataset"]: string indicating the high-res dataset
        ["name_basis_dataset"]: string indicating the basis set

    df_offsets: dataframe containing the offset values and high-res dataset names
        ["name_highres_dataset"]: name indicating high-res study
        ["offset_highres_dataset_residuals"]: offset value to add to Fe/H values

    plot_string: string of plot filename (if make_plot=True)
    feh_basis_csv_dir: directory to contain the FYI csv of literature Fe/H values for making a basis
    make_plot: write plot or not

    OUTPUTS:
    d: dictionary with
        "m_merged_highres": slope of high_res_feh vs. feh_basis
        "b_merged_highres": y-intercept of " " " "
        "m_merged_shifted_resid": slope of high_res_feh residuals (i.e., high_res_feh minus feh_basis) vs. feh_basis
        "b_merged_shifted_resid": y-intercept of " " " "
    '''

    dict_merged_this_basis = {} # initialize dictionaries
    dict_not_merged_this_basis = {}

    # names of high-res datasets (which will also be the keys to dict_merged_this_basis)
    highres_names = df_to_offset["name_highres_dataset"].unique()

    # for each high-res dataset name, apply the offsets to Fe/H residuals
    for this_dataset_name in highres_names:

        this_dataset = df_to_offset[
            df_to_offset["name_highres_dataset"].str.match(this_dataset_name)
            ]

        # retrieve the required offset
        this_resid_offset = (df_offsets["offset_highres_residuals"].
                             loc[df_offsets["name_highres"] == this_dataset_name].values)

        # if there is an offset that could be found for this dataset
        if this_resid_offset:
            this_dataset.loc[:, "residuals_no_shift"] = (this_dataset["feh_highres"] -
                                                         this_dataset["feh_basis"])
            this_dataset.loc[:, "residuals_no_shift"] = (this_dataset["feh_highres"] -
                                                         this_dataset["feh_basis"])

            # add the residuals as found above
            this_dataset.loc[:, "residuals_shifted"] = (this_dataset["residuals_no_shift"] +
                                                        this_resid_offset)

            # fold this correction into the non-residual Fe/H values
            this_dataset.loc[:, "feh_highres_post_shift"] = (this_dataset["feh_highres"] +
                                                             this_resid_offset)

            # add dataframe to dictionary;
            # each key corresponds to a high-res dataset, and each value is
            # a dataframe (this is good for plotting)
            dict_not_merged_this_basis[this_dataset_name] = this_dataset

        else:
            continue

    # merge ["residuals_shifted"] and ["feh_highres"] across all dataframes
    # (i.e., across all high-res studies) in the dictionary (for finding net Fe/H mapping)
    vals_of_interest = dict_not_merged_this_basis.copy()
    pd_merged = pd.concat(vals_of_interest.values(), ignore_index=True)

    # A value is trying to be set on a copy of a slice from a DataFrame.
    # Try using .loc[row_indexer, col_indexer] = value instead

    # find best-fit line to Fe/H plot of high_res vs. basis
    # (note that user may have used a flag to make Fe/H values be offset)
    limits = [-3.0, 0.5] # Fe/H limits to display

    # regression line for high-res Fe/H
    # (note this is NOT before adding in the shifts)
    m_merged_highres, b_merged_highres = np.polyfit(pd_merged["feh_basis"],
                                                    pd_merged["feh_highres"], 1)
    line_highres = np.multiply(m_merged_highres, limits) + b_merged_highres

    # regression line for high-res Fe/H AFTER adding the shifts back in
    m_merged_highres_postshift, b_merged_highres_postshift = np.polyfit(pd_merged["feh_basis"],
                                                                        pd_merged["feh_highres_post_shift"], 1)

    # regression line for merged residuals
    m_merged_shifted_resid, b_merged_shifted_resid = np.polyfit(pd_merged["feh_basis"],
                                                                pd_merged["residuals_shifted"], 1)

    # make best-fit line for residuals
    line_shifted_resid = (np.multiply(m_merged_shifted_resid, limits) +
                          b_merged_shifted_resid)

    # save the data to a csv
    pd_merged.to_csv(feh_basis_csv_dir +
                     csv_string +
                     config_red["file_names"]["MERGED_LIT_FEH_CSV"])

    # save a plot (raw high_res vs. basis on top; shifted residuals vs. basis in middle;
    # corrected high_res vs. basis on bottom)
    plt.clf()
    limits = [-3.0, 0.5]
    fig, axs = plt.subplots(3, 1, figsize=(10, 15), sharex=True)

    ## plot 0: raw high-res Fe/H

    # 1-to-1 line
    axs[0].plot([limits[0], limits[1]], [limits[0], limits[1]], linestyle="--")
    # best-fit line
    axs[0].plot([limits[0], limits[1]], np.add(np.multiply(m_merged_highres,
                                                         [limits[0], limits[1]]),
                                             b_merged_highres),
                                             linestyle="--")
    # input vs. basis
    axs[0].scatter(pd_merged["feh_basis"], pd_merged["feh_highres"])
    # for keeping data points separated by datasets, plot with a list comprehension
    [axs[0].scatter(dict_not_merged_this_basis[key]["feh_basis"],
                    dict_not_merged_this_basis[key]["feh_highres"],
                    label=key) for key in dict_not_merged_this_basis]
    # add star names
    for point in range(0, len(pd_merged["feh_basis"])):
        axs[0].annotate(
            pd_merged["name_star"][point],
            xy=(pd_merged["feh_basis"][point],
                pd_merged["feh_highres"][point]),
            xytext=(pd_merged["feh_basis"][point]+0.1,
                    pd_merged["feh_highres"][point]+0.06),
            textcoords="data", ha="right", va="bottom",
            fontsize=10,
            arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=0"))

    axs[0].set_xlim(limits[0], limits[1])
    axs[0].set_ylabel("Fe/H, high-res; no offsets applied")
    axs[0].set_title("m = "+str(m_merged_highres)+
                     ", b = "+
                     str(b_merged_highres)+
                     "; (blue line: 1-to-1; orange line: best fit)")
    axs[0].legend() # indicates high-res dataset names

    ## plot 1: shifted residuals
    axs[1].axhline(y=0, linestyle="--") # dashed line at y=0
    # input vs. basis
    axs[1].scatter(df_to_offset["feh_basis"], np.subtract(df_to_offset["feh_highres"],
                                                          df_to_offset["feh_basis"]))
    # input vs. basis
    axs[1].scatter(pd_merged["feh_basis"], pd_merged["residuals_shifted"], alpha=0.5)
    # for keeping data points separated by datasets, plot with a list comprehension
    [axs[1].scatter(dict_not_merged_this_basis[key]["feh_basis"],
                    dict_not_merged_this_basis[key]["residuals_shifted"],
                    label=key) for key in dict_not_merged_this_basis]
    # add star names
    for point in range(0, len(pd_merged["feh_basis"])):
        axs[1].annotate(
            pd_merged["name_star"][point],
            xy=(pd_merged["feh_basis"][point],
                pd_merged["residuals_shifted"][point]),
            xytext=(pd_merged["feh_basis"][point]+0.1,
                    pd_merged["residuals_shifted"][point]+0.06),
            textcoords="data", ha="right", va="bottom",
            fontsize=10,
            arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=0"))
    axs[1].set_xlabel("Fe/H, Basis")
    axs[1].set_ylabel("Fe/H Residuals: high-res minus basis set")
    axs[1].set_title("m = "+str(m_merged_shifted_resid)+
                     ", b = "+
                     str(b_merged_shifted_resid)+
                     "; (blue line: zero)")
    axs[1].set_ylim([-1.5, 1.5])
    axs[1].legend() # indicates high-res dataset names


    ## plot 2: shifted high-res Fe/H (the basis for the calibration!)

    # 1-to-1 line
    axs[2].plot([limits[0], limits[1]], [limits[0], limits[1]], linestyle="--")
    # input vs. basis
    axs[2].plot([limits[0], limits[1]],
                np.add(
                    np.multiply(m_merged_highres_postshift, [limits[0], limits[1]]),
                    b_merged_highres_postshift
                    ),
                    linestyle="--") # best-fit line
    axs[2].scatter(pd_merged["feh_basis"], pd_merged["feh_highres_post_shift"])
    # for keeping data points separated by datasets, plot with a list comprehension
    [axs[2].scatter(dict_not_merged_this_basis[key]["feh_basis"],
                    dict_not_merged_this_basis[key]["feh_highres_post_shift"],
                    label=key) for key in dict_not_merged_this_basis]
    # add star names
    for point in range(0, len(pd_merged["feh_basis"])):
        axs[2].annotate(
            pd_merged["name_star"][point],
            xy=(pd_merged["feh_basis"][point],
                pd_merged["feh_highres_post_shift"][point]),
            xytext=(pd_merged["feh_basis"][point]+0.1,
                    pd_merged["feh_highres_post_shift"][point]+0.06),
            textcoords="data", ha="right", va="bottom",
            fontsize=10,
            arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=0"))

    axs[2].set_xlim(limits[0], limits[1])
    axs[2].set_ylabel("Fe/H, high-res; WITH offsets applied")
    axs[2].set_title("m = "+
                     str(m_merged_highres_postshift)+
                     ", b = "+
                     str(b_merged_highres_postshift)+
                     "; (blue line: 1-to-1; orange line: best fit)")
    axs[2].legend() # indicates high-res dataset names

    if make_plot:
        plt.savefig(plot_string, overwrite=True)

    plt.clf()

    print("Scatter in residuals before offsets:")
    print(np.std(np.subtract(df_to_offset["feh_highres"],
                             df_to_offset["feh_basis"])))
    print("Number of data points:")
    print(len(np.subtract(df_to_offset["feh_highres"],
                          df_to_offset["feh_basis"])))
    print("-----")
    print("Scatter in residuals after offset shifts:")
    print(np.std(pd_merged["residuals_shifted"]))
    print(len(pd_merged["residuals_shifted"]))

    # best-fit line coeffs for high-res vs. basis
    #d["coeff_merged_highres"] = [m_merged_highres, b_merged_highres]
    # best-fit line coeffs for (residuals: high-res minus basis) vs. basis
    #d["coeff_merged_resid"] = [m_merged_resid, b_merged_resid]

    d = {
        "m_merged_highres": m_merged_highres,
        "b_merged_highres": b_merged_highres,
        "m_merged_shifted_resid": m_merged_shifted_resid,
        "b_merged_shifted_resid": b_merged_shifted_resid,
        "m_merged_highres_postshift": m_merged_highres_postshift,
        "b_merged_highres_postshift": b_merged_highres_postshift
    }

    return d


def calc_feh_program_stars(pickle_subdir=config_red["data_dirs"]["DIR_PICKLE"]):
    '''
    Calculate metallicities for the program stars which form
    the basis of the metallicity calibration, by using the
    remapping relationships found in make_basis_via_offsets()

    INPUTS:
    basis_set: basis set used for either RRab (such as Layden 1994)
        or RRc (such as Kemper+ 1982)
    '''

    # instantiate literature metallicity values
    lit_metallicity_basis = LitMetallicities()

    # RRab, RRc matches between high-res studies and basis sets
    rrab_matches = lit_metallicity_basis.match_highres_w_basis("RRab")
    rrc_matches = lit_metallicity_basis.match_highres_w_basis("RRc")

    # find necessary offsets to feh_highres-feh_basis, for each high-res study
    print("======= STEP 1a: CALCULATE RRAB OFFSETS ========")
    rrab_offsets = return_offsets(data_postmatch=rrab_matches)
    print("======= STEP 1a: PRINT RRAB OFFSETS ========")
    print(rrab_offsets)
    print("======= STEP 2a: CALCULATE RRC OFFSETS ========")
    rrc_offsets = return_offsets(data_postmatch=rrc_matches)
    print("======= STEP 2b: PRINT RRC OFFSETS ========")
    print(rrc_offsets)

    # get best-fit line info for the mappings, based on RRL sub-type and origin of offsets
    print("======= STEP 3: MAKE RRAB BASIS W RRAB OFFSETS ========")
    rrab_basis_w_rrab_offsets = make_basis_via_offsets(df_to_offset=rrab_matches,
                                                           df_offsets=rrab_offsets,
                                                           plot_string=(config_red["data_dirs"]["DIR_FYI_INFO"]+
                                                                        "rrab_w_rrab_offsets.png"),
                                                           csv_string="rrab_w_rrab_offsets",
                                                           make_plot=True)
    print("======= STEP 4: MAKE RRAB BASIS W RRC OFFSETS ========")
    rrab_basis_w_rrc_offsets = make_basis_via_offsets(df_to_offset=rrab_matches,
                                                          df_offsets=rrc_offsets,
                                                          plot_string=(config_red["data_dirs"]["DIR_FYI_INFO"]+
                                                                       "rrab_w_rrc_offsets.png"),
                                                          csv_string="rrab_w_rrc_offsets",
                                                          make_plot=True)
    print("======= STEP 5: MAKE RRC BASIS W RRC OFFSETS ========")
    rrc_basis_w_rrc_offsets = make_basis_via_offsets(df_to_offset=rrc_matches,
                                                         df_offsets=rrc_offsets,
                                                         plot_string=(config_red["data_dirs"]["DIR_FYI_INFO"]+
                                                                      "rrc_w_rrc_offsets.png"),
                                                         csv_string="rrc_w_rrc_offsets",
                                                         make_plot=True)
    print("======= STEP 6: MAKE RRC BASIS W RRAB OFFSETS ========")
    rrc_basis_w_rrab_offsets = make_basis_via_offsets(df_to_offset=rrc_matches,
                                                          df_offsets=rrab_offsets,
                                                          plot_string=(config_red["data_dirs"]["DIR_FYI_INFO"]+
                                                                       "rrc_w_rrab_offsets.png"),
                                                          csv_string="rrc_w_rrab_offsets",
                                                          make_plot=True)

    # use the bases to put Fe/H values on a common scale
    # i.e., to have ONE high-res-spectroscopically determined
    # Fe/H value for making the metallicity calibration

    # mapping is
    # [Fe/H]_highres_convention = m*[Fe/H]_basis_set + b

    # RRabs, using ab-based offsets
    rrab_feh_highres_ab_offsets = np.add(np.multiply(rrab_basis_w_rrab_offsets["m_merged_highres"],
                                                     rrab_matches["feh_basis"]),
                                         rrab_basis_w_rrab_offsets["b_merged_highres"])
    print(rrab_feh_highres_ab_offsets)
    # RRabs, using c-based offsets
    rrab_feh_highres_c_offsets = np.add(np.multiply(rrab_basis_w_rrc_offsets["m_merged_highres"],
                                                    rrab_matches["feh_basis"]),
                                         rrab_basis_w_rrc_offsets["b_merged_highres"])
    print(rrab_feh_highres_c_offsets)
    # RRcs, using ab-based offsets
    rrc_feh_highres_ab_offsets = np.add(np.multiply(rrc_basis_w_rrab_offsets["m_merged_highres"],
                                                    rrc_matches["feh_basis"]),
                                        rrc_basis_w_rrab_offsets["b_merged_highres"])
    print(rrc_feh_highres_ab_offsets)
    # RRcs, using c-based offsets
    rrc_feh_highres_c_offsets = np.add(np.multiply(rrc_basis_w_rrc_offsets["m_merged_highres"],
                                                   rrc_matches["feh_basis"]),
                                       rrc_basis_w_rrc_offsets["b_merged_highres"])
    print(rrc_feh_highres_c_offsets)

    # RRabs with RRab offsets
    plot_mapping(input_matches=rrab_matches,
                     mapped=rrab_feh_highres_ab_offsets,
                     title_string="RRab Fe/H mapping, w/ ab-based offsets",
                     plot_file_name="rrab_w_ab_offsets_basis.png")

    # RRabs with RRc offsets
    plot_mapping(input_matches=rrab_matches,
                     mapped=rrab_feh_highres_c_offsets,
                     title_string="RRab Fe/H mapping, w/ c-based offsets",
                     plot_file_name="rrab_w_c_offsets_basis.png")

    # RRcs with RRab offsets
    plot_mapping(input_matches=rrc_matches,
                     mapped=rrc_feh_highres_ab_offsets,
                     title_string="RRc Fe/H mapping, w/ ab-based offsets",
                     plot_file_name="rrc_w_ab_offsets_basis.png")

    # RRabs with RRab offsets
    plot_mapping(input_matches=rrc_matches,
                     mapped=rrc_feh_highres_c_offsets,
                     title_string="RRc Fe/H mapping, w/ c-based offsets",
                     plot_file_name="rrc_w_c_offsets_basis.png")

    pickle.dump( [rrab_matches, rrab_feh_highres_ab_offsets],
                     open( pickle_subdir + config_red["file_names"]["RRAB_RRAB_OFFSETS"], "wb" ) )
    pickle.dump( [rrab_matches, rrab_feh_highres_c_offsets],
                     open( pickle_subdir + config_red["file_names"]["RRAB_RRC_OFFSETS"], "wb" ) )
    pickle.dump( [rrc_matches, rrc_feh_highres_ab_offsets],
                     open( pickle_subdir + config_red["file_names"]["RRC_RRAB_OFFSETS"], "wb" ) )
    pickle.dump( [rrc_matches, rrc_feh_highres_c_offsets],
                     open( pickle_subdir + config_red["file_names"]["RRC_RRC_OFFSETS"], "wb" ) )

    # PRINT m (slope) and b (y-intercept) info for the mapping ## ## FYI ONLY, AT THE MOMENT
    print("rrab_basis_w_rrab_offsets:")
    print(rrab_basis_w_rrab_offsets)
    print("rrab_basis_w_rrc_offsets:")
    print(rrab_basis_w_rrc_offsets)
    print("rrc_basis_w_rrc_offsets:")
    print(rrc_basis_w_rrc_offsets)
    print("rrc_basis_w_rrab_offsets:")
    print(rrc_basis_w_rrab_offsets)

    return


# entry point
if __name__ == '__main__':
    sys.exit(main())
