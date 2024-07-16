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
import pandas as pd
from datetime import date, datetime
import json
from glob import glob

def run_bagle(target, phot_data, modstr):
    # Load needed data from munge
    data = mad_munge.getdata2(target,
                          phot_data=phot_data,
                          ast_data = [])
    priors = mad_munge.getpriors(target)
    print(priors)
    datestr = datetime.today().strftime('%Y_%m_%d')

    # Set up directory structure if needed
    if not os.path.isdir('bagle_fits/'+target+'/fits_'+datestr+'/'+modstr):
        os.makedirs('bagle_fits/'+target+'/fits_'+datestr+'/'+modstr)

    outdir = 'bagle_fits/'+target+'/fits_'+datestr+'/'+modstr+'/'
    outbase = 'b0_'
    os.makedirs(outdir, exist_ok=True)
    
    if modstr=='pspl_phot_par':
        model_selection = model.PSPL_Phot_Par_Param2
    elif modstr=='bspl_phot_par':
        model_selection = model.BSPL_Phot_Par_Param1
    elif modstr=='bspl_phot_nopar':
        model_selection = model.BSPL_Phot_noPar_Param1

    fitter = model_fitter.PSPL_Solver(data,
                                      model_selection,
                                      importance_nested_sampling = False,
                                      n_live_points = 400,
                                      evidence_tolerance = 0.5,
                                      sampling_efficiency = 0.8,
                                      outputfiles_basename=outdir + outbase)

    # Adjust the priors to encompass both possible solutions
    fitter.priors['t0'] = model_fitter.make_gen(*priors['t0'])
    fitter.priors['u0_amp'] = model_fitter.make_gen(-1.0, 1.0)
    fitter.priors['tE'] = model_fitter.make_gen(*priors['tE'])
    fitter.priors['piE_E'] = model_fitter.make_norm_gen(-0.02, 0.12)
    fitter.priors['piE_N'] = model_fitter.make_norm_gen(-0.03, 0.13)
    fitter.priors['b_sff1'] = model_fitter.make_gen(0.0,1.05)
    #fitter.priors['b_sff2'] = model_fitter.make_gen(0, 1.01)
    if modstr[:4]=='pspl':
        fitter.priors['mag_base1'] = model_fitter.make_gen(*priors['Ibase'])
    #fitter.priors['mag_src2'] = model_fitter.make_gen(18, 21)
    
    if modstr[:4]=='bspl':
        fitter.priors['sep'] = model_fitter.make_gen(0.0,2.0)
        fitter.priors['phi'] = model_fitter.make_gen(0.0,360.0)
        fitter.priors['mag_src_pri1'] = model_fitter.make_gen(priors['Ibase'][0], priors['Ibase'][1]+3)
        fitter.priors['mag_src_sec1'] = model_fitter.make_gen(priors['Ibase'][0], priors['Ibase'][1]+3)
        
    if target[0]=='K' or target[0]=='M':
        fitter.priors['mag_base1'] = model_fitter.make_gen(data['med_mag1']-1.0,data['med_mag1']+1.0)
        print(data['med_mag1']-1.0,data['med_mag1']+1.0)

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
                                   zoomx=[[60310, 60780], None, None],
                                   zoomy=[[np.max(data['mag']), np.min(data['mag'])], None, None],
                                   zoomy_res=[[-0.5, 0.5], None, None])
        plt.close('all')

def run_all(modstr='pspl_phot_par'):
    # Get required inputs for fit
    query_output = json.load(open(sorted(glob('query_output*'))[-1]))
    target_list = list(query_output['ra'].keys())
    #file = open('ignore_events_list.txt','r')
    ignored_events = []#file.readlines()
    #file.close()
    for target in target_list:
        #if (target not in ignored_events) and (target[0]!='O')and (target[0]!='M'):
        try:
            print(target)
            run_bagle(target, list(query_output['data_sets'][target].keys()), modstr)
        #else:
        except:
            print('skip', target)

def run_one(target, modstr='pspl_phot_par'):
    query_output = json.load(open(sorted(glob('query_output*'))[-1]))
    run_bagle(target, list(query_output['data_sets'][target].keys()), modstr)

run_all() #one('MB24028')  #('KB241458')
