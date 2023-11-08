from bagle import model_fitter
from bagle import model
#from bagle import munge
#from mad 
import mad_munge
from bagle import multinest_utils
from mpi4py import MPI
import argparse
import time
import os
import sys
import pylab as plt
import numpy as np
from datetime import datetime
import json

def run_bagle(target, phot_data, modstr):
    # Load needed data from munge
    data = mad_munge.getdata2(target,
                          phot_data=phot_data,
                          ast_data = [])
    priors = mad_munge.getpriors(target)
    datestr = datetime.today().strftime('%Y_%m_%d')

    # Set up directory structure if needed
    if not os.path.isdir('bagle_fits/'+target+'/fits_'+datestr+'/'+modstr):
        os.makedirs('bagle_fits/'+target+'/fits_'+datestr+'/'+modstr)

    outdir = 'bagle_fits/'+target+'/fits_'+datestr+'/'+modstr+'/'
    outbase = 'b0_'
    os.makedirs(outdir, exist_ok=True)

    fitter = model_fitter.PSPL_Solver(data,
                                      model.PSPL_Phot_Par_Param1,
                                      importance_nested_sampling = False,
                                      n_live_points = 400,
                                      evidence_tolerance = 0.5,
                                      sampling_efficiency = 0.8,
                                      outputfiles_basename=outdir + outbase)

    # Adjust the priors to encompass both possible solutions
    fitter.priors['t0'] = model_fitter.make_gen(*priors['t0'])
    fitter.priors['u0_amp'] = model_fitter.make_gen(-1.0, 1.0)
    fitter.priors['tE'] = model_fitter.make_gen(*priors['tE'])
    fitter.priors['piE_E'] = model_fitter.make_gen(-1, 1)
    fitter.priors['piE_N'] = model_fitter.make_gen(-1, 1)
    fitter.priors['b_sff1'] = model_fitter.make_gen(*priors['srcfrac'])
    #fitter.priors['b_sff2'] = model_fitter.make_gen(0, 1.01)
    fitter.priors['mag_src1'] = model_fitter.make_gen(*priors['Isrc'])
    #fitter.priors['mag_src2'] = model_fitter.make_gen(18, 21)

    ##########
    # BELOW: Do not change.
    #    Standard code for all runs.
    ##########

    # Process arguments to see if we are solving or not.
    parser = argparse.ArgumentParser()
    parser.add_argument('-N', '--nosolve', help='run plotting only, no solving', action='store_true')

    args = parser.parse_args()

    if (args.nosolve == False):
        t0 = time.time()
     
        fitter.solve()
     
        t1 = time.time()
        print('Runtime: {0:.0f} sec'.format(t1 - t0))

        
    # Catch the case where calling with MPI. Only use a single thread here.
    comm = MPI.COMM_WORLD
    comm.Barrier()

    if comm.Get_rank() == 0:
        print('Making plots')
        fitter.summarize_results()
        fitter.separate_modes()
        fitter.summarize_results_modes()
        fitter.write_summary_maxL()
        fitter.plot_dynesty_style(fit_vals=None)
        fitter.plot_dynesty_style_hel_to_geo_phot(60142, fit_vals='maxL')
        plt.close('all')
        best_mod = fitter.get_best_fit_model(def_best='maxL')
        fitter.plot_model_and_data(best_mod, suffix='_maxL',
                                   zoomx=[[59900, 60300], None, None],
                                   zoomy=[[16.8, 15], None, None],
                                   zoomy_res=[[-0.1, 0.1], None, None])
        plt.close('all')

def run_all():
    # Get required inputs for fit
    modstr = 'pspl_phot_par'
    query_output = json.load(open('query_output.json'))
    target_list = list(query_output['ra'].keys())
    for target in target_list:
        run_bagle(target, list(query_output['data_sets'][target].keys()), modstr)

run_all()
