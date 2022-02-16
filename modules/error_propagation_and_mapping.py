'''
Finds errors for program star metallicities
'''

import pickle
import time
import multiprocessing
import scipy
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
from multiprocessing import Pool
#from rrlfe.modules2 import *

class feh_plotter():
    '''
    Class containing a bunch of the functions we will use to map metallicities
    '''

    def __init__(self):

        pass

    def __call__(self):

        pass

    def cdf_fcn(self, array_input):
        '''
        Return CDF of an unsorted input array of values
        '''

        number_cum_norm = np.divide(np.arange(len(array_input)),
                                    len(array_input))
        array_input_sort = np.sort(array_input)
        #array_cdf = np.divide(np.cumsum(array_input_sort),
        #                      np.cumsum(array_input_sort)[-1])

        return array_input_sort, number_cum_norm


    def cdf_gauss(self, x_range, mu, sig):
        '''
        Returns the CDF of a general Gaussian function for any mu and sig
        '''

        # rescale x -> x* = (x-mu)/sig
        x_range_adjust = np.divide(np.subtract(x_range, mu), sig)

        # erf(x*/sqrt(2))
        erf_return = scipy.special.erf(np.divide(x_range_adjust, np.sqrt(2)))

        # return (1/2)*(1 + erf(x*/sqrt(2)))
        return np.multiply(0.5, np.add(1., erf_return))


    def pickle_plot_info(self,
                         name_star,
                         feh_mapped_array,
                         write_pickle_subdir=config["data_dirs"]["DIR_PICKLE"]):
        '''
        Find sigmas and pickle the info

        INPUTS:
        name_star: string defining the star
        feh_mapped_array: list of Fe/H values for this star, post-mapping
        write_pickle_subdir: directory to write the pickled Fe/H info to
        '''

        x_vals, y_vals = self.cdf_fcn(np.ravel(feh_mapped_array))

        # fit a Gaussian
        popt, pcov = optimize.curve_fit(self.cdf_gauss, x_vals, y_vals)

        print("Line parameters")
        print(popt)

        xvals_interp = np.arange(x_vals[0], x_vals[-1], 0.001)
        yvals_interp = np.interp(xvals_interp, x_vals, y_vals)

        # find element of value closest to 0.5
        percent_bar = 0.5
        idx = np.abs(yvals_interp - percent_bar).argmin()


        # ---------------------------------------------------------------------------

        # SIGMA DEFINITION 1: FIND MEDIAN AND SIGMA BRACKETS AROUND IT
        # find element of value closest to 1-sigma limit (on low side)
        one_sigma_perc = 0.682689492
        percent_bar_1sig_low = 0.5 - 0.5*one_sigma_perc
        idx_1sig_low = np.abs(yvals_interp - percent_bar_1sig_low).argmin()
        # find element of value closest to 1-sigma limit (on high side)
        percent_bar_1sig_high = 0.5 + 0.5*one_sigma_perc
        idx_1sig_high = np.abs(yvals_interp - percent_bar_1sig_high).argmin()

        # SIGMA DEFINITION 2: FIND NARROWEST REGION CONTAINING 1-SIGMA WORTH OF POINTS
        shortest_xrange = xvals_interp[-1] - xvals_interp[0]
        shortest_xrange_lower = xvals_interp[0]

        for t in range(0, len(xvals_interp)):

            lower_bar_y = yvals_interp[t]
            upper_bar_y = yvals_interp[t] + one_sigma_perc

            # break if the range will go beyond data points
            if (upper_bar_y > 0.99):
                break

            idx_1sig_here = t
            idx_1sig_above = np.abs(yvals_interp - upper_bar_y).argmin()

            if (np.subtract(xvals_interp[idx_1sig_above],
                        xvals_interp[idx_1sig_here]) < shortest_xrange):

                shortest_xrange = xvals_interp[idx_1sig_above] - xvals_interp[idx_1sig_here]
                shortest_xrange_lower = xvals_interp[idx_1sig_here]
                shortest_xrange_upper = xvals_interp[idx_1sig_above]
                shortest_xrange_halfway = (0.5*np.subtract(shortest_xrange_upper,
                                                           shortest_xrange_lower) +
                                           shortest_xrange_lower)

        # ---------------------------------------------------------------------------

        print("Fe/H at 50 percentile")
        feh_50_perc = xvals_interp[idx]
        print(feh_50_perc)

        print("1-sigma interval")
        feh_1sig_low = xvals_interp[idx_1sig_low]
        feh_1sig_high = xvals_interp[idx_1sig_high]
        print(feh_1sig_low)
        print(feh_1sig_high)

        # pickle the data for this one star, to avoid choking the machine
        # with too much plot-making all at once

        # replace space with underscore
        name_star_underscore = str(name_star).replace(" ", "_")
        pickle_write_name = (write_pickle_subdir + "plot_info_" +
                             name_star_underscore + ".pkl")
        cdf_gauss_info = self.cdf_gauss(x_vals, *popt)
        with open(pickle_write_name, "wb") as f:
                pickle.dump((name_star,
                     feh_mapped_array,
                     x_vals,
                     y_vals,
                     xvals_interp,
                     cdf_gauss_info,
                     idx,
                     idx_1sig_low,
                     idx_1sig_high,
                     shortest_xrange_lower,
                     shortest_xrange_upper,
                     shortest_xrange_halfway), f)

        # return FeH based on definition 1 (median and sigma brackets)
        # and definition 2 (narrowest region containing 1-sigma worth of points)
        return (feh_1sig_low,
                feh_50_perc,
                feh_1sig_high,
                shortest_xrange_lower,
                shortest_xrange_halfway,
                shortest_xrange_upper)


    def write_cdf_hist_plot(self,
                            name_star,
                            read_pickle_subdir=config["data_dirs"]["DIR_PICKLE"],
                            write_plot_subdir=config["data_dirs"]["DIR_FYI_INFO"],
                            write_plot=True):
        '''
        Takes the pickled plot info and saves CDF and histogram plots

        INPUTS:
        name_star: string ID of the star
        read_pickle_subdir: directory to read the pickled Fe/H info from
        write_plot_subdir: directory to write the FYI plots to
        write_plot: write plot or not
        '''

        print("Making CDF and histogram plots of FeH for " + name_star + "...")

        # replace space with underscore
        name_star_underscore = str(name_star).replace(" ", "_")
        pickle_read_name = (read_pickle_subdir + "plot_info_" +
                            name_star_underscore + ".pkl")

        # open the pickle file
        with open(pickle_read_name, 'rb') as f:
            name_star,feh_mapped_array,x_vals,y_vals,xvals_interp,cdf_gauss_info,\
              idx,idx_1sig_low,idx_1sig_high,shortest_xrange_lower,\
              shortest_xrange_upper,shortest_xrange_halfway = pickle.load(f)

        # if no plot is to be written
        if write_plot == False:
            return

        plt.clf()
        plt.plot(x_vals, y_vals)
        plt.plot(x_vals, cdf_gauss_info, linestyle=":", color="k")
        plt.axvline(xvals_interp[idx], color='blue')
        plt.axvline(xvals_interp[idx_1sig_low], color='blue')
        plt.axvline(xvals_interp[idx_1sig_high], color='blue')
        plt.axvline(shortest_xrange_lower, color='orange')
        plt.axvline(shortest_xrange_upper, color='orange')
        plt.axvline(shortest_xrange_halfway, color='orange')
        plt.xlabel("Fe/H")
        plt.ylabel("CDF")
        plt.title(name_star + "\n" +
                  "Fe/H based on median (blue): " +
                  "{:.{}f}".format( xvals_interp[idx], 3) +
                  ", +" +
                  "{:.{}f}".format( np.subtract(xvals_interp[idx_1sig_high],
                                                xvals_interp[idx]), 3) +
                  ", -" + "{:.{}f}".format( np.subtract(xvals_interp[idx],
                                                        xvals_interp[idx_1sig_low]), 3) +
                  "\n"+
                  "Fe/H based on shortest range (orange): " +
                  "{:.{}f}".format( shortest_xrange_halfway, 3) +
                  ", +" +
                  "{:.{}f}".format( np.subtract(shortest_xrange_upper,
                                                shortest_xrange_halfway), 3) +
                  ", -" +
                  "{:.{}f}".format( np.subtract(shortest_xrange_halfway,
                                                shortest_xrange_lower), 3))
        plt.tight_layout()
        plt.savefig(write_plot_subdir + name_star + "_cdf.pdf")
        plt.close()

        plt.clf()
        plt.hist(np.ravel(feh_mapped_array), bins=100)
        plt.title(name_star + "\n" + "std = "+str(np.std(np.ravel(feh_mapped_array))))
        plt.xlabel("Fe/H")
        plt.tight_layout()
        plt.savefig(write_plot_subdir + name_star + "_hist.pdf")
        plt.close()



    def do_bootstrap(self,
                     read_pickle_subdir=config["data_dirs"]["DIR_PICKLE"]):
        '''
        Do bootstrap on high-res Fe/H values to find the mapping

        INPUTS:
        read_pickle_subdir: directory containing the pickle file of Fe/H data

        OUTPUTS:
        m_array: array of slopes for each bootstrap step
        b_array: array of y-intercepts for each bootstrap step
        params_list_star_feh: star names and basis set Fe/H ## ## for what?
        data_1: original data which is fed into the bootstrap
        '''

        # read in actual data
        ## ## N.b. this is just the RRabs with RRab offsets for now
        real_data_1 = pickle.load( open( read_pickle_subdir
                                         + config["file_names"]["RRAB_RRAB_OFFSETS"], "rb" ) )

        # arrange the data in a way we can use
        # N.b. This is NOT fake data; I'm just appropriating the old variable name
        ## ## Note the ersatz Layden errors for now; need to revisit this with values from his paper
        data_1 = { "star_name": real_data_1[0]["name_star"],
                "feh_lit": real_data_1[0]["feh_highres"],
                "feh_layden": real_data_1[0]["feh_basis"],
                "err_feh_lit": np.zeros(len(real_data_1[0]["feh_basis"])),
                "err_feh_layden": 0.07*np.ones(len(real_data_1[0]["feh_basis"]))}
        #dataset_1 = pd.DataFrame(data=data_1)


        # # Find the linear regression line to high res literature Fe/H vs. basis set Fe/H values

        # Put Fe/H values into a useable form
        feh_sample = np.transpose([data_1["feh_layden"], data_1["feh_lit"]])

        # Bootstrap
        N_samples = int(1e4)
        # set RNG for reproducibility of the bootstrap
        with NumpyRNGContext(1):
            bootresult = bootstrap(feh_sample, N_samples)

        # populate the arrays with bootstrap results
        m_array = np.nan*np.ones(len(bootresult)) # initialize
        b_array = np.nan*np.ones(len(bootresult))
        for boot_n in range(0, len(bootresult)):
            test_fit = np.polyfit(bootresult[boot_n, :, 0], bootresult[boot_n, :, 1], 1)
            m_array[boot_n] = test_fit[0]
            b_array[boot_n] = test_fit[1]

        # consolidate info, remove extra dimension
        name_star = data_1["star_name"][:].values
        feh_test = data_1["feh_layden"][:].values
        params_array = np.squeeze([[name_star], [feh_test]])

        # arrange into a list for parallel processing
        params_list_star_feh = list(np.transpose(params_array))

        return m_array, b_array, params_list_star_feh, data_1


class feh_mapper(feh_plotter):

    def __init__(self,
                 write_pickle_subdir=config["data_dirs"]["DIR_PICKLE"]):

        self.write_pickle_subdir = write_pickle_subdir

    def __call__(self):
        pass

    def map_feh_one_star(self, params_element):
        '''
        Maps the Fe/H values for one star (in a single function for parallel processing)
        '''

        time_start = time.time()
        feh_mapped_array = np.nan*np.ones(len(m_array)) # initialize array

        # get name and Layden Fe/H of star
        name_star = params_element[:][0]
        feh_test = params_element[:][1]
        print("Star:")
        print(name_star)
        print("Layden Fe/H:")
        print(feh_test)

        for sample_num in range(0, len(m_array)):

            feh_mapped_1sample = m_array[sample_num]*feh_test + b_array[sample_num]
            feh_mapped_array[sample_num] = feh_mapped_1sample

        ## take one basis set Fe/H (integral over a Gaussian)
        ## and find what the mapped value should be

        # number of samples to take within the Gaussian error around Layden's Fe/H value
        N = 100
        gaussian_spread = 0.07 ## ## change this in future
        layden_feh = feh_test # this is the discrete value

        # N_m_samples x N_Layden_samples
        feh_mapped_array = np.nan*np.ones((len(m_array), N))

        # loop over each sample within the Gaussian around Layden's Fe/H
        for integal_piece in range(0, N):

            # set the offset (note mu=0; this is a relative offset)
            offset = np.random.normal(loc=0.0, scale=gaussian_spread)

            # loop over all (m,b) combinations found further above
            for sample_num in range(0, len(m_array)):

                ## ## is it layden_feh*(1. + offset) or (layden_feh + offset) ?
                feh_mapped_1sample = (m_array[sample_num]*layden_feh*(1. + offset) +
                                      b_array[sample_num])
                feh_mapped_array[sample_num][integal_piece] = feh_mapped_1sample

        # pickle plot info
        self.pickle_plot_info(name_star,
                              feh_mapped_array,
                              write_pickle_subdir=self.write_pickle_subdir)

        print("Elapsed time:")
        print(str(time.time() - time_start))
        print("--------------------------")


    def do(self):

        # do the bootstrap
        global m_array # have to make this global for multiprocessing to work
        global b_array
        m_array, b_array, params_list_star_feh, data_1 = self.do_bootstrap()

        # parallel process the Fe/H info
        #ncpu = multiprocessing.cpu_count()
        pool = Pool(ncpu)
        outdat = pool.map(self.map_feh_one_star, params_list_star_feh) # FeH info is pickled here
        pool.close()

        # now make and save the plots of Fe/H values
        # (this is done in series to avoid memory chokes)
        for t in range(0, len(data_1["star_name"])):
            this_star = data_1["star_name"][t]
            print("Writing Fe/H CDF for star " + this_star)
            self.write_cdf_hist_plot(this_star)
