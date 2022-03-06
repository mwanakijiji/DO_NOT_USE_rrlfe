'''
This is the reduction pipeline for applying the updated Layden coefficients
[a, b, c, d] to low-resolution spectra
'''

import sys
from conf import *
from modules import *
from modules import (compile_normalization,
                      create_spec_realizations,
                      run_robo,
                      scrape_ew_and_errew,
                      find_feh)

def main():

    # set the objective
    # "apply_calib": apply pre-determined coefficients to spectra to find Fe/H
    # "find_abcd": determine coefficients [a,b,c,d] in the first place
    objective_choice = "apply_calib"

    # Make all the directories
    make_dirs(objective = objective_choice)

    # Compile the C spectral normalization script
    compile_normalization.compile_bkgrnd()

    # Take list of unnormalized empirical science spectra and normalize them
    # (N.b. for application of solution, only make 1 realization, with 0 noise;
    # note also that application of the calibration uses non-default directories)
    create_spec_realizations.create_spec_realizations_main(num = 1,
                                                            spec_file_type="ascii.no_header",
                                                            noise_level="None",
                                                            input_spec_list_dir = config_apply["data_dirs"]["DIR_SRC"],
                                                            input_list = config_apply["data_dirs"]["DIR_SRC"] + config_apply["file_names"]["LIST_SPEC_APPLY"],
                                                            unnorm_empirical_spectra_dir = config_apply["data_dirs"]["DIR_SCI_SPECTRA"],
                                                            unnorm_noise_churned_spectra_dir = config_apply["data_dirs"]["DIR_REZNS_SPEC"],
                                                            bkgrnd_output_dir = config_apply["data_dirs"]["DIR_REZNS_SPEC_NORM"],
                                                            final_dir = config_apply["data_dirs"]["DIR_REZNS_SPEC_NORM_FINAL"])

    # run_robospect on normalized spectra
    # (N.b. Robospect directory is with 'config' because that's where the program was installed)
    run_robo.main(
                normzed_spec_source_dir = config_apply["data_dirs"]["DIR_REZNS_SPEC_NORM_FINAL"],
                write_dir=config_apply["data_dirs"]["DIR_ROBO_OUTPUT"],
                robo_dir=config_apply["sys_dirs"]["DIR_ROBO"]
                )

    # scrape_ew_from_robo and calculate EWs + err_EW

    scraper_instance = scrape_ew_and_errew.Scraper(subdir = config_apply["data_dirs"]["DIR_ROBO_OUTPUT"],
                                                   file_scraped_info = config_apply["data_dirs"]["DIR_EW_PRODS"]+config_apply["file_names"]["SCRAPED_SCIENCE_SPECTRA_FILE_NAME"])
    scraper_instance() # call instance

    # follow-up functions
    data_checked = scrape_ew_and_errew.quality_check(
                    read_in_filename = config_apply["data_dirs"]["DIR_EW_PRODS"]+config_apply["file_names"]["SCRAPED_SCIENCE_SPECTRA_FILE_NAME"],
                    write_out_filename = config_apply["data_dirs"]["DIR_EW_PRODS"]+config_apply["file_names"]["SCRAPED_EW_DATA_GOOD_ONLY"]
                    )

    # put the good EW data into a table with
    # rows corresponding to files and cols for the lines
    data_stacked = scrape_ew_and_errew.stack_spectra(
        read_in_filename = config_apply["data_dirs"]["DIR_EW_PRODS"]+config_apply["file_names"]["SCRAPED_EW_DATA_GOOD_ONLY"],
        write_out_filename = config_apply["data_dirs"]["DIR_EW_PRODS"]+config_apply["file_names"]["RESTACKED_EW_DATA_GOOD_ONLY"],
        input_list = config_apply["data_dirs"]["DIR_SRC"]+config_apply["file_names"]["LIST_SPEC_APPLY"])


    data_net_balmer = scrape_ew_and_errew.generate_net_balmer(
        read_in_filename = config_apply["data_dirs"]["DIR_EW_PRODS"]+config_apply["file_names"]["RESTACKED_EW_DATA_GOOD_ONLY"],
        write_out_filename = config_apply["data_dirs"]["DIR_EW_PRODS"]+config_apply["file_names"]["RESTACKED_EW_DATA_W_NET_BALMER"])

    data_errors = scrape_ew_and_errew.generate_addl_ew_errors(
        read_in_filename = config_apply["data_dirs"]["DIR_EW_PRODS"]+config_apply["file_names"]["RESTACKED_EW_DATA_W_NET_BALMER"],
        write_out_filename = config_apply["data_dirs"]["DIR_EW_PRODS"]+config_apply["file_names"]["RESTACKED_EW_DATA_W_NET_BALMER_ERRORS"],
        groupby_parent = False)

    # find Fe/H values, sampling from the a, b, c, d posteriors and while
    # incorporating equivalent width errors
    find_feh_instance = find_feh.findFeH()

    # find Fe/H and write
    find_feh_instance.pickle_feh_retrieval(write_out_filename = config_apply["data_dirs"]["DIR_BIN"]+config_apply["file_names"]["RETRIEVED_VALS"])
    '''
    # retrieve pickle files and compare values (only for case of synthetic values with injected and retrieved Fe/H)
    find_feh_instance.compare_feh_synthetic()
    '''

# entry point
if __name__ == '__main__':
    sys.exit(main())
