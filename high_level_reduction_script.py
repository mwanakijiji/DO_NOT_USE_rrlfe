'''
This is the high-level script which runs all the pieces of the pipeline to
obtain updated Layden coefficients [a, b, c, d]
'''

import sys
from conf import *
from modules import *
from modules import (compile_normalization,
                      create_spec_realizations,
                      run_robo,
                      scrape_ew_and_errew,
                      ca_correction,
                      consolidate_pre_mcmc,
                      run_emcee,
                      teff_retrieval)

def main():

    model_choice = "abcdfghk"
    '''
    # make all the directories
    make_dirs(objective = "find_calib") ## find_calib as opposed to apply_calib

    # compile the C spectral normalization script

    compile_normalization.compile_bkgrnd()

    # Take list of unnormalized empirical spectra and noise-churned the
    # spectra, normalize them, and write out normalizations
    ## ## just 1 or 2 realizations for testing (default is 100)

    create_spec_realizations.create_spec_realizations_main(num = 100, noise_level=0.07, spec_file_type="ascii.no_header")

    # run_robospect on normalized synthetic spectra
    run_robo.main()

    # scrape_ew_from_robo and calculate EWs + err_EW
    scraper_instance = scrape_ew_and_errew.Scraper()
    scraper_instance() # call instance
    '''
    data_checked = scrape_ew_and_errew.quality_check()

    # put the good EW data into a table with
    # rows corresponding to files and cols for the lines

    data_stacked = scrape_ew_and_errew.stack_spectra()

    data_net_balmer = scrape_ew_and_errew.generate_net_balmer()

    data_errors = scrape_ew_and_errew.generate_addl_ew_errors(groupby_parent = True)

    data_add_metadata = scrape_ew_and_errew.add_synthetic_meta_data()

    # finds the Teff calibration
    temp = teff_retrieval.temp_vs_balmer()
    '''
    # run_emcee with input data_table_winnowed
    # coeff defs: K = a + bH + cF + dHF + f(H^2) + g(F^2) + h(H^2)F + kH(F^2) + m(H^3) + n(F^3)
    # where K is CaII K EW; H is Balmer EW; F is [Fe/H]

    emcee_instance = run_emcee.RunEmcee()
    #emcee_instance(model = 'abcd') # call instance
    emcee_instance(model = model_choice)

    posterior_write = run_emcee.write_soln_to_fits(model = model_choice)

    posterior_sample = run_emcee.corner_plot(model = model_choice)
    '''

# entry point
if __name__ == '__main__':
    sys.exit(main())
