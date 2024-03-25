import numpy as np
#import pylab as plt
#from microlens.jlu import model
#from microlens.jlu import model_fitter
#import dynesty
#from dynesty import utils as dyutil
#from dynesty import plotting as dyplot
from astropy.table import Table
from astropy.time import Time
from astropy import units 
from astropy.coordinates import SkyCoord
from astropy import time as atime, coordinates as coord, units as u
#from multiprocessing import Pool, cpu_count
#import time
#import pickle
#import pdb
import os
from datetime import date
import json

mad_dir = os.getcwd()+'/' #'/u/mhuston/code/MAD/'
data = json.load(open(mad_dir+'query_output_' + str(date.today()) + '.json'))

ra = data['ra']
dec = data['dec']
photom_ogle = data['photom_ogle']
for photom in photom_ogle:
	photom_ogle.update({photom:mad_dir+photom_ogle[photom]})
photom_moa = data['photom_moa']
for photom in photom_moa:
	photom_moa.update({photom:mad_dir+photom_moa[photom]})
photom_kmt = data['photom_kmt']
for photom in photom_kmt:
	photom_kmt.update({photom:mad_dir+photom_kmt[photom]})
data_set_info = data['data_sets']
data_sets = {}
for event_id in data_set_info:
    data_sets[event_id] = {}
    for data_set in data_set_info[event_id]:
        data_sets[event_id][data_set] = data[data_set_info[event_id][data_set]][event_id]
print(data_sets)

def getpriors(target):
    # Get alert fit info as starting points to set priors for BAGLE fit
    #alertkeys = ['t0', 't0_err', 'tE', 'tE_err', 'Ibase', 'Ibase_err', 'Isrc', 'Isrc_err', 'srcfrac', 'srcfrac_err']
    alertkeys = ['t0', 'tE', 'Ibase']
    alertfit = {}
    for key in alertkeys:
        alertfit[key] = data[key][target]

    # Calculate reasonable priors based on alert fit
    #priorkeys = ['t0', 'tE', 'Isrc', 'srcfrac']
    priors = {}

    priors['t0'] = [alertfit['t0'] - alertfit['tE']/2, alertfit['t0'] + alertfit['tE']/2]
    priors['tE'] = [min(alertfit['tE']/2, 50), alertfit['tE'] + alertfit['tE']/2]
    priors['Ibase'] = [alertfit['Ibase'] - 0.2, alertfit['Ibase'] + 0.2]

    # Don't put real limits on blending
    #priors['srcfrac'] = [0.001,1.05]

    print(priors)
    return priors

def getdata2(target, phot_data=['I_OGLE'], ast_data=['Kp_Keck'],
             time_format='mjd', verbose=False):
    """
    Get the photometric and astrometric data for the specified target. 
    Specify the types of data through the phot_data and ast_data lists. 

    Inputs
    ----------
    target : str
        Target name (lower case)

    Optional Inputs
    --------------------
    phot_data : list
        List of strings specifying the data sets. Options include:
        I_OGLE, Kp_Keck, Ch1_Spitzer, MOA

    ast_data : list
        List of strings specifying the data sets. Options include:
        Kp_Keck

    time_format : string
        The time format (default = 'mjd') such as mjd, year, jd.

    verbose : bool
        Print out extra information. 

    Returns
    ----------
    data : dict
        A ditionary containing the data. For each photometric data set, the dictionary
        contains:
        
            t_phot1
            mag1
            mag_err1

        where the '1' at the end is incremented for additional data sets. Note the
        index is assigned according to the order int he list. Note that if only a single
        photometry data set is requested, the returned keys are t_phot, mag, mag_err
        with a 1 on the end.

        For each astrometric data set, the dictionary contains:
        
            t_ast1
            xpos1
            ypos1
            xpos_err1
            ypos_err1

        where the index is incremented for additional astrometric data sets (useful for the future).

        There are two additional entries in the dictionary which contain the R.A. and Dec.
        of the lensing target... this is the photocenter of the joint lens/source system and 
        is hard-coded in module tables. Note that if only a single
        astrometry data set is requested, the returned keys are t_ast, xpos, ypos, xpos_err, ypos_err
        with no index on the end.

        data['raL']
        data['decL']
    """
    data = {}

    # Load the RA and Dec
    target_coords = SkyCoord(ra[target], dec[target], 
                             unit = (units.hourangle, units.deg), frame = 'icrs')
    data['target'] = target
    data['raL'] = target_coords.ra.degree
    data['decL'] = target_coords.dec.degree

    # Keep track of the data files we used.
    phot_files = []
    ast_files = []
    
    # Load up the photometric data.
    for pp in range(len(phot_data)):
        filt = phot_data[pp]

        if filt not in data_sets[target].keys():
            raise RuntimeError('Failed to find photometric data set {0:s} for {1:s}'.format(filt, target))
        
        phot_files.append(data_sets[target][filt])
                              
        if filt == 'I_OGLE':
            # Read in photometry table.
            pho = Table.read(data_sets[target][filt], format = 'ascii')
            t = Time(pho['mjd'], format='mjd', scale='utc')
            m = pho['mag']
            me = pho['mag_err']

        if filt == 'Kp_Keck':
            pho = Table.read(data_sets[target][filt])
            tdx = np.where(pho['name'] == target)[0][0]
            t = Time(pho['t'][tdx, :], format='jyear', scale='utc')
            m = pho['m'][tdx, :]

            # Add empirical photometric errors and radii
            pho['m0e'] = np.nanstd(pho['m'], axis=1)
            pho['dr'] = np.hypot(pho['x0'] - pho['x0'][tdx],
                                 pho['y0'] - pho['y0'][tdx])
            pho['dm'] = pho['m0'] - pho['m0'][tdx]
            
            # We can't use the normal me because it captures the source variability.
            # Calculate from surrounding stars. Don't use the target itself.
            dr = pho['dr']
            dm = np.abs(pho['dm'])

            # Iterate until we get some stars.
            n_neigh = 0
            nn = 1
            dr_start = 2.5
            dm_start = 1.5
            
            while n_neigh < 3:
                r_factor = 1.0 + (nn / 10.)  # Grow search radius by 10% each round.
                m_factor = 1.0 + (nn / 5.)   # Grow mag search by 20% each round.
                rdx = np.where((dr < dr_start*r_factor) & (dm < dm_start*m_factor) & (pho['m0e'] != 0))[0]
                rdx = rdx[rdx != tdx]  # Drop the target.
                n_neigh = len(rdx)
                nn += 1

            # For all the magniudes of the surrounding stars (in individual epochs),
            # mask out the invalid values.
            me_neigh = np.nanmean( pho['m0e'][rdx] )
            if verbose:
                print('Found {0:d} neighbors within:'.format(n_neigh))
                print('  dr = {0:0.2f} arcsec'.format(dr_start * r_factor))
                print('  dm = {0:0.2f} mag'.format(dm_start * m_factor))
                print(pho['name','m0','m0e','dr','dm'][rdx])

            if np.isnan(me_neigh):
                me_neigh = 0.025
                if verbose:
                    print('Using hard-coded me_neigh')

            if verbose: 
                print('me_neigh = {0:.3f} mag'.format(me_neigh))
                
            me = np.ones(len(t), dtype=float) * me_neigh

        #if filt == 'Ch1_Spitzer':
        #    pho = Table.read(data_sets[target][filt], format='ascii')
        #    t = Time(pho['col1']  + 2450000.0, format='jd', scale='utc')
        #    f = pho['col2']
        #    fe = pho['col3']
        #    m = 25.0 - 2.5 * np.log10(f)
        #    me = 1.086 * fe / f

        if filt == 'MOA':
            pho = Table.read(data_sets[target][filt], format='ascii')
            # Convert HJD provided by MOA into JD.
            # https://geohack.toolforge.org/geohack.php?pagename=Mount_John_University_Observatory&params=43_59.2_S_170_27.9_E_region:NZ-CAN_type:landmark
            moa = coord.EarthLocation(lat=-43.986667 * u.deg,lon=170.465*u.deg, height=1029*u.meter)
            #print(type(atime))
            t_hjd = atime.Time(pho['col1'], format='jd', scale = 'utc')
            ltt = t_hjd.light_travel_time(target_coords, 'heliocentric', location=moa)

            t = t_hjd - ltt
            m = pho['col5']
            me = pho['col6']

        if filt[0:3] == 'HST':
            pho = Table.read(data_sets[target][filt])
            if 'ob110462' in target:
                tdx = np.where(pho['name'] == 'OB110462')[0][0]
            elif target == 'ob110462_op_bc':
                tdx = np.where(pho['name'] == 'OB110462')[0][0]
            elif target == 'ob110462_new2':
                tdx = np.where(pho['name'] == 'OB110462')[0][0]
            else:
                tdx = np.where(pho['name'] == target.upper())[0][0]
            good_idx = np.where(~np.isnan(pho[tdx]['t']))[0] # get rid of nans
            t = Time(pho['t'][tdx, good_idx], format='jyear', scale='utc')
            m = pho['m'][tdx, good_idx]
            me = pho['me'][tdx, good_idx]
            # Make sure t is increasing
            if t[0] > t[-1]:
                t = t[::-1]
                m = m[::-1]
                me = me[::-1]

        if filt == 'KMT':
            pho = Table.read(data_sets[target][filt], format = 'ascii')
            t = Time(pho['mjd'], format='mjd', scale='utc')
            m = pho['mag']
            me = pho['mag_err']

        if filt == 'KMT_DIA':
            pho = Table.read(data_sets[target][filt], format='ascii')
            t = Time(pho['col1'] + 2450000.0, format='jd', scale='utc')
            m = 27.68-2.5*np.log10(pho['col2']+27300)
            me = -1.08 * pho['col3']/(pho['col2'] + 27300)

        # Set time to proper format
        if time_format == 'mjd':
            t = t.mjd
        if time_format == 'jyear':
            t = t.j_year
        if time_format == 'jd':
            t = t.jd

        # Insert the data into the dictionary.
        suffix = '{0:d}'.format(pp + 1)
        if len(phot_data) == 1:
            suffix = '1'

        data['t_phot' + suffix] = t
        data['mag' + suffix] = m
        data['mag_err' + suffix] = me

    for aa in range(len(ast_data)):
        filt = ast_data[aa]

        if filt not in data_sets[target].keys():
            raise RuntimeError('Failed to find astrometric data set {0:s} for {1:s}'.format(filt, target))

        ast_files.append(data_sets[target][filt])
        
        if filt == 'Kp_Keck':
            ast = Table.read(data_sets[target][filt])
            tdx = np.where(ast['name'] == target)[0][0]
            t = Time(ast['t'][tdx, :], format='jyear', scale='utc')
            x = ast['x'][tdx, :] * -1.0   # East in +x direction
            y = ast['y'][tdx, :]
            xe = ast['xe'][tdx, :]
            ye = ast['ye'][tdx, :]

        if filt[0:3] == 'HST':
            ast = Table.read(data_sets[target][filt])
            if 'ob110462' in target:
                tdx = np.where(ast['name'] == 'OB110462')[0][0]
            elif target == 'ob110462_op_bc':
                tdx = np.where(ast['name'] == 'OB110462')[0][0]
            elif target == 'ob110462_new2':
                tdx = np.where(ast['name'] == 'OB110462')[0][0]
            else:
                tdx = np.where(ast['name'] == target.upper())[0][0]
            good_idx = np.where(~np.isnan(ast[tdx]['t']))[0] # get rid of nans 
            t = Time(ast['t'][tdx, good_idx], format='jyear', scale='utc')
            x = ast['x'][tdx, good_idx] * -1.0   # East in +x direction
            y = ast['y'][tdx, good_idx]
            xe = ast['xe'][tdx, good_idx]
            ye = ast['ye'][tdx, good_idx]

            # Make sure t is increasing
            if t[0] > t[-1]:
                t = t[::-1]
                x = x[::-1]
                y = y[::-1]
                xe = xe[::-1]
                ye = ye[::-1]
            
        # Set time to proper format
        if time_format == 'mjd':
            t = t.mjd
        if time_format == 'jyear':
            t = t.j_year
        if time_format == 'jd':
            t = t.jd

        # Insert the data into the dictionary.
        suffix = '{0:d}'.format(aa + 1)
            
        data['t_ast' + suffix] = t
        data['xpos' + suffix] = x
        data['ypos' + suffix] = y
        data['xpos_err' + suffix] = xe
        data['ypos_err' + suffix] = ye

    # Keep a record of the types of data.
    data['phot_data'] = phot_data
    data['ast_data'] = ast_data

    data['phot_files'] = phot_files
    data['ast_files'] = ast_files
            
    return data

def getdata_OGLEphot(phot_file, ra, dec, time_format='mjd', add_tconst = True):
    """
    Set up OGLE photometry to be fit, given some file.
    (This is mostly for OGLE EWS stuff where we don't
    have any astrometry available.)
    """
    # Read in photometry table.
    pho = Table.read(phot_file, format = 'ascii')
    pho.rename_column('col1', 't')
    pho.rename_column('col2', 'm')
    pho.rename_column('col3', 'me')

    # Put all the times in MJD and fix tables.
    if add_tconst is True:
        # For OGLE-III EWS
        p_t = Time(pho['t'] + 2450000, format = 'jd', scale = 'utc')

    if add_tconst is False:
        # For OGLE-IV EWS
        p_t = Time(pho['t'], format = 'jd', scale = 'utc')
    
    if time_format=='mjd':
        pho['t'] = p_t.mjd

    # Prepare data to be fit
    target_coords = SkyCoord(ra, dec, 
                             unit = (units.hourangle, units.deg), frame = 'icrs')
    
    data = {}
    data['raL'] = target_coords.ra.degree
    data['decL'] = target_coords.dec.degree
    data['t_phot1'] = pho['t']
    data['mag1'] = pho['m']
    data['mag_err1'] = pho['me']
        
    return data 

