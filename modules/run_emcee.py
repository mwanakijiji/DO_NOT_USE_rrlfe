'''
This is an emcee wrapper for fitting the Layden '94 metallicity
calibration to equivalent widths of RR Lyrae spectra
'''

import sys
import git
import time
import numpy as np
import pandas as pd
import emcee
import corner
import logging
from scipy.optimize import least_squares
from astropy.io import fits
from . import *


def corner_plot(model,
                mcmc_text_output_file_name = config_red["data_dirs"]["DIR_BIN"] + config_red["file_names"]["MCMC_OUTPUT"],
                corner_plot_putput_file_name = config_red["data_dirs"]["DIR_BIN"] + config_red["file_names"]["MCMC_CORNER"]):
    '''
    Reads in MCMC output and writes out a corner plot
    '''

    if (model == "abcd"):

        # corner plot (requires 'storechain=True' in enumerate above)
        test_samples = pd.read_csv(mcmc_text_output_file_name, delim_whitespace=True, nrows=5) # read in first rows to check column number
        samples = pd.read_csv(mcmc_text_output_file_name, delim_whitespace=True, usecols=(1,2,3,4), names=["a", "b", "c", "d"])
        fig = corner.corner(samples, labels=["$a$", "$b$", "$c$", "$d$"],
                            quantiles=[0.16, 0.5, 0.84],
                            title_fmt='.2f',
                            show_titles=True,
                            verbose=True,
                            title_kwargs={"fontsize": 12})
        fig.savefig(corner_plot_putput_file_name)
        logging.info("--------------------------")
        logging.info("Corner plot of MCMC posteriors written out to")
        print(corner_plot_putput_file_name)

        # if its necessary to read in MCMC output again
        #data = np.loadtxt(self.mcmc_text_output, usecols=range(1,5))

        # This code snippet from Foreman-Mackey's emcee documentation, v2.2.1 of
        # https://emcee.readthedocs.io/en/stable/user/line.html#results
        a_mcmc, b_mcmc, c_mcmc, d_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                             zip(*np.percentile(samples, [16, 50, 84], axis=0)))

        print("--------------------------")
        print("Coefficients a, b, c, d, and errors (see corner plot):")
        print("coeff a: " + " ".join(map(str,a_mcmc)))
        print("coeff b: " + " ".join(map(str,b_mcmc)))
        print("coeff c: " + " ".join(map(str,c_mcmc)))
        print("coeff d: " + " ".join(map(str,d_mcmc)))
        logging.info("--------------------------")
        logging.info("Coefficients a, b, c, d, and errors (see corner plot):")
        logging.info("coeff a: " + " ".join(map(str,a_mcmc)))
        logging.info("coeff b: " + " ".join(map(str,b_mcmc)))
        logging.info("coeff c: " + " ".join(map(str,c_mcmc)))
        logging.info("coeff d: " + " ".join(map(str,d_mcmc)))

        #logging.info("--------------------------")
        #logging.info("MCMC data written to ")
        #logging.info(self.mcmc_text_output)


    elif (model == "abcdfghk"):
        # corner plot (requires 'storechain=True' in enumerate above)
        # just first few lines to test
        test_samples = pd.read_csv(mcmc_text_output_file_name, delim_whitespace=True, nrows=5) # read in first rows to check column number
        samples = pd.read_csv(mcmc_text_output_file_name, delim_whitespace=True, usecols=(1,2,3,4,5,6,7,8), names=["a", "b", "c", "d", "f", "g", "h", "k"])
        fig = corner.corner(samples, labels=["$a$", "$b$", "$c$", "$d$", "$f$", "$g$", "$h$", "$k$"],
                            quantiles=[0.16, 0.5, 0.84],
                            title_fmt='.2f',
                            show_titles=True,
                            verbose=True,
                            title_kwargs={"fontsize": 12})
        fig.savefig(corner_plot_putput_file_name)
        logging.info("--------------------------")
        logging.info("Corner plot of MCMC posteriors written out to")
        print(str(corner_plot_putput_file_name))

        # if its necessary to read in MCMC output again
        #data = np.loadtxt(self.mcmc_text_output, usecols=range(1,5))

        # This code snippet from Foreman-Mackey's emcee documentation, v2.2.1 of
        # https://emcee.readthedocs.io/en/stable/user/line.html#results
        a_mcmc, b_mcmc, c_mcmc, d_mcmc, f_mcmc, g_mcmc, h_mcmc, k_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                             zip(*np.percentile(samples, [16, 50, 84], axis=0)))

        print("--------------------------")
        print("Coefficients a, b, c, d, f, g, h, k, and errors (see corner plot):")
        print("coeff a: " + " ".join(map(str,a_mcmc)))
        print("coeff b: " + " ".join(map(str,b_mcmc)))
        print("coeff c: " + " ".join(map(str,c_mcmc)))
        print("coeff d: " + " ".join(map(str,d_mcmc)))
        print("coeff f: " + " ".join(map(str,f_mcmc)))
        print("coeff g: " + " ".join(map(str,g_mcmc)))
        print("coeff h: " + " ".join(map(str,h_mcmc)))
        print("coeff k: " + " ".join(map(str,k_mcmc)))
        logging.info("--------------------------")
        logging.info("Coefficients a, b, c, d, f, g, h, k, and errors (see corner plot):")
        logging.info("coeff a: " + " ".join(map(str,a_mcmc)))
        logging.info("coeff b: " + " ".join(map(str,b_mcmc)))
        logging.info("coeff c: " + " ".join(map(str,c_mcmc)))
        logging.info("coeff d: " + " ".join(map(str,d_mcmc)))
        logging.info("coeff f: " + " ".join(map(str,f_mcmc)))
        logging.info("coeff g: " + " ".join(map(str,g_mcmc)))
        logging.info("coeff h: " + " ".join(map(str,h_mcmc)))
        logging.info("coeff k: " + " ".join(map(str,k_mcmc)))


    else: # pragma: no cover

        logging.error("Error! No calibration model chosen for the MCMC posteriors!")

    return test_samples


def lnprob(walker_pos_array,
           Teff_pass,
           measured_H_pass,
           measured_F_pass,
           measured_K_pass,
           err_measured_H_pass,
           err_measured_F_pass,
           err_measured_K_pass):
    '''
    Nat log of probability density

    INPUTS:
    walker_pos_array: array of coefficients (regardless of model)
    Teff_pass: Teff (a vestigial MCMC constant; this is NOT astrophysical Teff)
    measured_H_pass: Balmer EW
    measured_F_pass: [Fe/H]
    measured_K_pass: CaIIK EW
    err_measured_H_pass: error in Balmer EW
    err_measured_F_pass: error in [Fe/H]
    err_measured_K_pass: error in CaIIK EW


    OUTPUTS:
    ln(prior*like)
    '''

    # walker_pos is the proposed walker position in N-D (likely 4-D) space
    # (i.e., these are the inputs to the model)
    lp = lnprior(walker_pos_array) # prior
    if not np.isfinite(lp): # afoul of prior
      return -np.inf
    result = -np.divide(1, 2*Teff_pass)*chi_sqd_fcn(measured_H_pass,
                                                    measured_F_pass,
                                                    measured_K_pass,
                                                    err_measured_H_pass,
                                                    err_measured_F_pass,
                                                    err_measured_K_pass,
                                                    walker_pos_array)

    return lp + result # ln(prior*like)



def lnprior(theta):
    '''
    Prior

    INPUTS:
    theta: array of parameter values

    OUTPUTS: 0 or -inf (top-hat priors only)
    '''

    if (len(theta) == 4):
        # Layden '94 relation
        a_test, b_test, c_test, d_test = theta
    elif (len(theta) == 8):
        # updated relation
        a_test, b_test, c_test, d_test, f_test, g_test, h_test, k_test = theta

    # top-hat priors
    if ((np.abs(a_test) < 40) and
        (np.abs(b_test) < 5) and
        (np.abs(c_test) < 20) and
        (np.abs(d_test) < 10)):
        return 0.0

    return -np.inf


def function_K(coeffs_pass,Bal_pass,F_pass):
    '''
    Function which gives CaIIK EW as function of Balmer, [Fe/H]

    INPUTS:
    coeffs_pass: array of coefficients
    Bal_pass: Balmer EWs
    F_pass: [Fe/H]

    OUTPUTS:
    K_pass: CaIIK EW
    '''

    if (len(coeffs_pass) == 4):
        # Layden '94 relation

        K_pass = coeffs_pass[0] \
                    + coeffs_pass[1]*Bal_pass \
                    + coeffs_pass[2]*F_pass \
                    + coeffs_pass[3]*Bal_pass*F_pass

    elif (len(coeffs_pass) == 8):
        # updated relation

        K_pass = coeffs_pass[0] \
                    + coeffs_pass[1]*Bal_pass \
                    + coeffs_pass[2]*F_pass \
                    + coeffs_pass[3]*Bal_pass*F_pass \
                    + coeffs_pass[4]*np.power(Bal_pass,2.) \
                    + coeffs_pass[5]*np.power(F_pass,2.) \
                    + coeffs_pass[6]*np.power(Bal_pass,2.)*F_pass \
                    + coeffs_pass[7]*Bal_pass*np.power(F_pass,2.)

    return K_pass


def K_residual(coeffs_pass,Bal_pass,F_pass,y):
    # residual function, which we want to minimize for the Levenberg-Marquardt fit
    # ('y' here is the measured CaIIK EW)

    return y - function_K(coeffs_pass,Bal_pass,F_pass)


def sigma_Km_sqd(coeffs_pass,Bal_pass,err_Bal_pass,Feh_pass,err_Feh_pass):
    # def of model CaIIK error squared (this is general, regardless of number of coeffs):
    # sigma_Km^2 = (del_K/del_H)^2 * sig_H^2 + (del_K/del_F)^2 * sig_F^2

    # case of 4 coefficients
    if (len(coeffs_pass) == 4):

        dKdH = coeffs_pass[1] + coeffs_pass[3]*Feh_pass
        dKdF = coeffs_pass[2] + coeffs_pass[3]*Bal_pass

    # case of 8 coefficients
    elif (len(coeffs_pass) == 8):

        dKdH = coeffs_pass[1] + coeffs_pass[3]*Feh_pass + 2.*coeffs_pass[4]*Bal_pass + \
                    2.*coeffs_pass[6]*Feh_pass*Bal_pass + coeffs_pass[7]*np.power(Feh_pass,2.)
        dKdF = coeffs_pass[2] + coeffs_pass[3]*Bal_pass + 2.*coeffs_pass[5]*Feh_pass + \
                    coeffs_pass[6]*np.power(Bal_pass,2.) + 2.*coeffs_pass[7]*Bal_pass*Feh_pass

    else: # pragma: no cover
        logging.error("Number of coefficients does not make sense!")

    # 'sigma_Km squared'
    sigma_Km_2 = np.power(dKdH*err_Bal_pass,2.) + np.power(dKdF*err_Feh_pass,2.)

    return sigma_Km_2


def chi_sqd_fcn(Bal_pass,
                Feh_pass,
                Caiik_pass,
                sig_Bal_pass,
                sig_Feh_pass,
                sig_Caiik_pass,
                coeffs_pass):
    '''
    Chi-squared

    INPUTS:
    Bal_pass: Balmer EW (angstroms)
    Feh_pass: [Fe/H]
    Caiik_pass: CaIIK EW (angstroms)
    err_Bal_pass: error in Balmer EW (angstroms)
    err_Feh_pass: error in [Fe/H]
    err_Caiik_pass: error in CaIIK EW (angstroms)
    coeffs_pass: array of coefficients

    OUTPUTS:
    val: chi^2
    '''

    # def. of chi-squared for collection of datapoints which each have
    # subscript i:
    # X2 = Sigma_i [(K0,i - Km,i)^2/(sigma_K0,i^2 + sigma_Km,i^2)]
    # K0: measured CaIIK EW (error sigma_K0)
    # Km: model CaIIK EW (error sigma_Km)

    # note the functions below accept  arrays of either 4 or 8 coefficients
    Caiik_model_array = function_K(coeffs_pass=coeffs_pass,
                                    Bal_pass=Bal_pass,
                                    F_pass=Feh_pass)

    sigma_Km_sqd_array = sigma_Km_sqd(coeffs_pass=coeffs_pass,
                                    Bal_pass=Bal_pass,
                                    err_Bal_pass=sig_Bal_pass,
                                    Feh_pass=Feh_pass,
                                    err_Feh_pass=sig_Feh_pass)

    numerator_array = np.power(np.subtract(Caiik_pass,Caiik_model_array),2.)

    denominator_array = np.add(sigma_Km_sqd_array,np.power(sig_Caiik_pass,2.))

    i_element = np.divide(numerator_array,denominator_array)

    val = np.sum(i_element)

    return val


def write_soln_to_fits(model,
                        mcmc_text_output_file_name = config_red["data_dirs"]["DIR_BIN"] + config_red["file_names"]["MCMC_OUTPUT"],
                        teff_data_retrieve_file_name = config_red["data_dirs"]["DIR_BIN"] + config_red["file_names"]["TREND_TEFF_VS_BALMER"],
                        soln_write_name = config_red["data_dirs"]["DIR_BIN"] + config_red["file_names"]["CALIB_SOLN"],
                        test_flag=False):
    '''
    Takes the full reduction solution and writes it to a FITS file with
    1) The MCMC posteriors in tabular form
    2) Meta-data in FITS header

    test_flag: if True, then terminal prompts are suppressed to enable continuous integration
    '''

    # initialize FITS header and append keys
    hdr = fits.Header()

    # retrieve git hash
    repo = git.Repo(search_parent_directories=True)
    sha = repo.head.object.hexsha
    hdr["HASH"] = (sha, "Hash of rrlfe pipeline")

    # get Teff vs Balmer line info
    # set compound datatype
    dtype=np.rec.fromrecords([['string_key', 189.6752158]]).dtype
    # load data, skipping header and hash corresponding to that file
    teff_data = np.loadtxt(teff_data_retrieve_file_name, skiprows=2, usecols=(0,1), delimiter=':', dtype=dtype)

    dict_teff_data = {}
    for key, val in teff_data:
        dict_teff_data.update({key: val})

    hdr["MODEL"] = (model, "Calibration type")
    hdr["TEFF_MIN"] = (dict_teff_data["Teff_min"], "Minimum Teff for linear Teff vs. Balmer EW fit")
    hdr["TEFF_MAX"] = (dict_teff_data["Teff_max"], "Maximum Teff for linear Teff vs. Balmer EW fit")
    hdr["SLOPE_M"] = (dict_teff_data["m"], "Slope of Teff vs. Balmer EW")
    hdr["ESLOPE_M"] = (dict_teff_data["err_m"], "Error in slope of Teff vs. Balmer EW")
    hdr["YINT_B"] = (dict_teff_data["b"], "Y-intercept of Teff vs. Balmer EW")
    hdr["EYINT_B"] = (dict_teff_data["err_b"], "Error in y-intercept of Teff vs. Balmer EW")
    # comment explaining the solution
    hdr["COMMENT"] = "Coefficients are defined as "
    hdr["COMMENT"] = "K = a + bH + cF + dHF + f(H^2) + g(F^2) + h(H^2)F + kH(F^2)"

    # read in posterior in csv form
    if (model == "abcd"):

        # corner plot (requires 'storechain=True' in enumerate above)
        samples = pd.read_csv(mcmc_text_output_file_name, delim_whitespace=True, usecols=(1,2,3,4), names=["a", "b", "c", "d"])
        c1 = fits.Column(name="a", array=np.array(samples.iloc[:,0].values), format="D")
        c2 = fits.Column(name="b", array=np.array(samples.iloc[:,1].values), format="D")
        c3 = fits.Column(name="c", array=np.array(samples.iloc[:,2].values), format="D")
        c4 = fits.Column(name="d", array=np.array(samples.iloc[:,3].values), format="D")
        table_hdu = fits.BinTableHDU.from_columns([c1, c2, c3, c4], header=hdr)


    elif (model == "abcdfghk"):
        # corner plot (requires 'storechain=True' in enumerate above)
        # just first few lines to test
        samples = pd.read_csv(mcmc_text_output_file_name, delim_whitespace=True, usecols=(1,2,3,4,5,6,7,8), names=["a", "b", "c", "d", "f", "g", "h", "k"])
        c1 = fits.Column(name="a", array=np.array(samples.iloc[:,0].values), format="D")
        c2 = fits.Column(name="b", array=np.array(samples.iloc[:,1].values), format="D")
        c3 = fits.Column(name="c", array=np.array(samples.iloc[:,2].values), format="D")
        c4 = fits.Column(name="d", array=np.array(samples.iloc[:,3].values), format="D")
        c5 = fits.Column(name="f", array=np.array(samples.iloc[:,4].values), format="D")
        c6 = fits.Column(name="g", array=np.array(samples.iloc[:,5].values), format="D")
        c7 = fits.Column(name="h", array=np.array(samples.iloc[:,6].values), format="D")
        c8 = fits.Column(name="k", array=np.array(samples.iloc[:,7].values), format="D")
        table_hdu = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7, c8], header=hdr)


    # write out
    if not test_flag: # pragma: no cover
        if os.path.exists(soln_write_name):

            input("A calibration solution file already exists! \nDo " +\
                    "what you want with that file and hit [ENTER] (will overwrite)")

        table_hdu.writeto(soln_write_name, overwrite=True)
        logging.info("Full calibration MCMC posterior written to " + soln_write_name)

    # return FITS table for testing
    return table_hdu



class RunEmcee():
    '''
    Run the emcee MCMC to obtain coefficients a, b, c, d (+ f, g, h, k)
    '''

    def __init__(self,
                 scraped_ews_good_only_file_name = config_red["data_dirs"]["DIR_EW_PRODS"] + config_red["file_names"]["RESTACKED_EW_DATA_GOOD_ONLY_TEFFFIT"],
                 mcmc_text_output_file_name = config_red["data_dirs"]["DIR_BIN"] + config_red["file_names"]["MCMC_OUTPUT"]
                 ):

        # name of file with final K, H, FeH, and error values (and not the others from the noise-churned spectra)
        self.scraped_ew_filename = scraped_ews_good_only_file_name

        # name of file of the MCMC output
        self.mcmc_text_output = mcmc_text_output_file_name


    def __call__(self, model, post_burn_in_links = 3e5):
        '''
        INPUTS

        model: list of coefficients to use as the model
            'abcd':     corresponds to Layden '94
            'abcdfghk': corresponds to K = a + b*H + c*F + d*H*F + f*(H^2) + g*(F^2) + h*(H^2)*F + k*H*(F^2)
        post_burn_in_links: chain links following burn-in
        '''

        # read in EWs, Fe/Hs, phases, errors, etc.
        logging.info("--------------------------")
        logging.info("Reading in data from " + self.scraped_ew_filename)

        ## ## make df_choice.Spectrum -> df_choice["Spectrum etc.
        df_choice = pd.read_csv(self.scraped_ew_filename,delim_whitespace=False)

        # EWs in table are in angstroms and are mislabeled as mA (2020 Jan 12)
        name = df_choice['orig_spec_file_name']
        #caii = np.divide(df_choice['K'], 1000.)
        caii = df_choice['EW_CaIIK']
        #ecaii = np.divide(df_choice['err_K'], 1000.)
        ecaii = df_choice['err_EW_CaIIK_from_robo'] # might try other error sources later
        #ecaii = df_choice['err_EW_CaIIK']
        #ave = np.divide(df_choice['balmer'], 1000.)
        ave = df_choice['EW_Balmer']
        eave = df_choice['err_EW_Balmer_from_Robo']
        #eave = np.divide(df_choice['err_balmer'], 1000.)
        ## ## THE BELOW FEH VALUES NEED TO BE CHECKED/FIXED
        feh = df_choice['feh']
        efeh = df_choice['err_feh']

        # non-zero starting points for coefficients fghk (i.e., those beyond the Layden
        # model)
        f_init = 0.1
        g_init = 0.1
        h_init = 0.1
        k_init = 0.1
        m_init = 0. # 4th-order term; vestigial
        n_init = 0. # 4th-order term; vestigial
        #sigma_a_layden = 0.416
        #sigma_b_layden = 0.076
        #sigma_c_layden = 0.285
        #sigma_d_layden = 0.052

        # starting position, before adding a perturbation

        if model == 'abcd':
            # coeffs_pass = [a,b,c,d]

            nwalkers = 8 # number of MCMC chains (at least 2x number of parameters)

            param_array_0 = [float(a_layden),
                            float(b_layden),
                            float(c_layden),
                            float(d_layden)]

        elif model == 'abcdfghk':
            # coeffs_pass = [a,b,c,d,f,g,h,k]

            nwalkers = 16 # number of MCMC chains (at least 2x number of parameters)

            param_array_0 = [float(a_layden),
                            float(b_layden),
                            float(c_layden),
                            float(d_layden),
                            float(f_init),
                            float(g_init),
                            float(h_init),
                            float(k_init)]

        ################# do a Levenberg-Marquardt fit to check consistency #################
        # this is for abcd or abcdfghk calibrations, whichever is being done by the MCMC
        lstsq_soln = least_squares(K_residual, param_array_0, args=(ave, feh, caii), method="lm")
        logging.info("--------------------------")
        logging.info("Levenberg-Marquardt soln made. Vector array: ")
        logging.info(lstsq_soln.x)

        ################# MCMC setup #################

        logging.info("--------------------------")
        logging.info("Setting up MCMC ...")

        ndim = len(param_array_0) # dimensions of space to explore

        # convert the one starting point into a nwalkers*ndim array with gaussian-offset starting points
        p0 = [np.add(param_array_0,
                     np.multiply(param_array_0, 1e-4*np.random.randn(ndim))) for i in range(nwalkers)]

        # set up sampler
        sampler = emcee.EnsembleSampler(nwalkers,
                                        ndim,
                                        lnprob,
                                        args=[Teff, ave, feh, caii, eave, efeh, ecaii])

        # burn-in
        burn_in = 3e4
        posAfterBurn, prob, state = sampler.run_mcmc(p0, burn_in)

        # post-burn-in
        start_time = time.time()

        ################# SAVE PROGRESSIVELY TO TEXT FILE #################
        ## ## refer to these code snippets from Foreman-Mackey's website
        # IMPORTANT: sampler will only have memory of the last iteration if
        # storechain flag is set to False

        logging.info("--------------------------")
        logging.info("Saving MCMC chains to text file ...")

        # post-burn-in calculate and save iteratively
        f = open(self.mcmc_text_output, "w")
        f.close()
        progBarWidth = 30
        start_time = time.time()
        for i, result in enumerate(sampler.sample(posAfterBurn,
                                                  iterations=post_burn_in_links,
                                                  storechain=True)):
            position = result[0]
            f = open(self.mcmc_text_output, "a") # append
            for k in range(position.shape[0]): # loop over number of chains
                position_string = str(position[k]).strip("[]") # convert to string
                f.write("{0:4d} {1:s}\n".format(k, " ".join(str(p) for p in position[k])))
            n = int((progBarWidth+1) * float(i) / post_burn_in_links) # update progress bar
            sys.stdout.write("\r[{0}{1}]".format("#" * n, " " * (progBarWidth - n)))
            f.close()
        elapsed_time = time.time() - start_time
        sys.stdout.write(" Done!\n")
        sys.stdout.write("{0:s} {1:10d} {2:s}\n".format("Elapsed time: ",
                                                        int(elapsed_time), "sec"))
        logging.info("--------------------------")
        sys.stdout.write("MCMC chain data written out to\n")
        logging.info(self.mcmc_text_output)
