
import ftplib
from sqlalchemy import create_engine, text
from sqlalchemy.sql import select
from bs4 import BeautifulSoup
import numpy as np 
from io import BytesIO
import pandas as pd
import time 
from urllib.request import urlopen
import numexpr as ne
import requests
import time
import multiprocessing as mp
from itertools import repeat

# Setting up database stuff with SQLAlchemy.
engine = create_engine('sqlite:///microlensing.db')

def get_moa_lightcurves(year):
    """
    Function that grabs MOA lightcurves from the alert pages 
    and writes them to a table in the database.
    
    Note that we don't care about the delta flux photometry...
    what we really care about is the photometry in magnitudes.
    This function takes the reported calibration values from
    the MOA webpage, and converts the delta flux and flux error
    measurements into magnitude and magnitude errors.

    Parameters
    ----------
    year : int
        Year of the MOA alerts you want. 
        Valid choices are 2016 - 2022, inclusive.
        
    Outputs
    -------
    sqlite table called photometry in microlensing.db
    Columns are hjd (HJD - 245000), mag, mag_err, alert_name, and telescope.
    """
    # The delta flux measurements sometimes yield negative fluxes
    # after calibration. Ignore warnings so we don't have to deal
    # with the log10 complaining during the magnitude conversion.
    import warnings
    warnings.filterwarnings("ignore")

    # Go to the MOA alerts site and scrape the page.
    year = str(year)
    url = "http://www.massey.ac.nz/~iabond/moa/alert" + year + "/alert.php"
    response = urlopen(url)
    html = response.read()
    response.close()
    soup = BeautifulSoup(html,"html.parser")

    # Get a list of all the bulge microlensing alert directories.
    links = soup.find_all('a', href=True)
    alert_dirs = []
    for ii, link in enumerate(links):
        if 'BLG' in link.text:
            alert_dirs.append(links[ii]['href'])

    t0 = time.time() 
    # Go to the page for each bulge microlensing alert.
    for nn, alert_dir in enumerate(alert_dirs):
        # Scrape the page.
        url = "http://www.massey.ac.nz/~iabond/moa/alert" + year + "/" + alert_dir
        response = urlopen(url)
        html = response.read()
        response.close()
        soup = BeautifulSoup(html,"html.parser")

        # Get the magnitude and flux offsets, so we can convert
        # from delta flux to a magnitude.
        foo = soup.find('b').next_sibling
        moff = foo.split('=')[1].split('-')[0].strip(' ')
        bah = soup.find('sub').next_sibling
        foff = bah.split('+')[1].split(')')[0].strip(' ')

        # Convert those offsets from strings into floats
        m = ne.evaluate(moff)
        f = ne.evaluate(foff)

        # Grab the .dat file containing the photometry data into a pandas dataframe.
        url = "https://www.massey.ac.nz/~iabond/moa/alert" + year + "/fetchtxt.php?path=moa/ephot/phot-" + \
                alert_dir.strip('display.php?id=') + ".dat"
        bytes_data = requests.get(url).content
        df = pd.read_csv(BytesIO(bytes_data), 
                         delim_whitespace=True, skiprows=11, skipfooter=1, header=None, engine='python', 
                         names=['hjd', 'delta_flux', 'flux_err', 'foo1', 'foo2', 'foo3', 'foo4', 'foo5'])

        # Add columns for magnitude and magnitude error, using the conversion
        # values we just figured out.
        df['mag'] = m - 2.5*np.log10(df['delta_flux'] + f)
        df['mag_err'] = 1.09 * df['flux_err']/(df['delta_flux'] + f)
        
        # Add a column for the alert name (of the form MBYYNNN, YY=year, NNN=alert number)
        # and telescope (MOA)
        df['alert_name'] = 'MB' + year[2:] + str(nn + 1).zfill(3)  # need to make sure this always works.
        df['telescope'] = 'MOA'
        
        # Write HJD as HJD - 2450000 to match OGLE and KMTNet (less cumbersome digits)
        df['hjd'] -= 2450000

        # Get rid of all the nans which crop up during the conversion from delta flux to magnitude.
        df.dropna(axis='index', how='any', inplace=True)

        # Write out the HJD, mag, mag_err, telescope, and alert_name data into the table.
        cols = ['hjd', 'mag', 'mag_err', 'telescope', 'alert_name']
        df[cols].to_sql(con=engine, schema=None, name="photometry", if_exists="append", index=False)
    t1 = time.time() 
    
    print('Took {0:.2f} seconds'.format(t1-t0))
    
def get_ogle_lightcurves(year):
    """
    Function that grabs OGLE lightcurves from the alert website 
    and writes them to a table in the database.

    Parameters
    ----------
    year : int
        Year of the MOA alerts you want. 
        Valid choices are 2011 - 2029, inclusive.
        
    Outputs
    -------
    sqlite table called photometry in microlensing.db
    Columns are hjd (HJD - 245000), mag, mag_err, alert_name, and telescope.
    """
    # Go to the OGLE alert site and get the data with FTP.
    year = str(year)
    ftp = ftplib.FTP("ftp.astrouw.edu.pl")
    ftp.login()
    ftp.cwd("ogle/ogle4/ews/" + year + "/")
    
    # Figure out how many objects there are by counting the number of .tar.gz files 
    # in the directory (each .tar.gz file corresponds to an alert)
    nobj = sum('.tar.gz' in x for x in ftp.nlst())
 
    t0 = time.time() 
    for nn in np.arange(start=1, stop=nobj+1, step=1):
        # Grab the photometry for each alert.
        ftp.cwd("blg-" + str(nn).zfill(4))
        
        flo = BytesIO()
        ftp.retrbinary('RETR phot.dat', flo.write)
        flo.seek(0)
        df = pd.read_fwf(flo, header=0, 
                         names=['hjd', 'mag', 'mag_err', 'see', 'sky'], 
                         widths=[14, 7, 6, 5, 8])

        # Add a column for the alert name (of the form OBYYNNNN, YY=year, NNN=alert number)
        # and telescope (OGLE)
        df['alert_name'] = 'OB' + year[2:] + str(nn).zfill(4) 
        df['telescope'] = 'OGLE'
        
        # Write HJD as HJD - 2450000 (less cumbersome digits)
        df['hjd'] -= 2450000
        
        # Write out the HJD, mag, mag_err, telescope, and alert_name data into the table.
        cols = ['hjd', 'mag', 'mag_err', 'telescope', 'alert_name']
        df[cols].to_sql(con=engine, schema=None, name="photometry", if_exists="append", index=False)

        ftp.cwd("../")
    t1 = time.time() 
    ftp.close()
    
    print('Took {0:.2f} seconds'.format(t1-t0))
    
def get_kmtnet_lightcurves(year):
    """
    Function that grabs KMTNet lightcurves from the alert pages 
    and writes them to a table in the database.

    Parameters
    ----------
    year : int
        Year of the KMTNet alerts you want. 
        Valid choices are 2016 - 2022, inclusive.
        
    Outputs
    -------
    sqlite table called photometry in microlensing.db
    Columns are hjd (HJD - 245000), mag, mag_err, alert_name, and telescope (the pysis name).
    """
    # Figure out how many objects there are by counting how many columns
    # there are on the alert page.
    year = str(year)
    url = "https://kmtnet.kasi.re.kr/~ulens/event/" + year + "/"
    response = urlopen(url)
    html = response.read()
    response.close()
    soup = BeautifulSoup(html,"html.parser")
    nobj = len(soup.find_all('td')[0::15][1:])
    
    t0 = time.time()
    # Go to the KMTNet alerts site and scrape the page for each alert.
    for nn in np.arange(start=1, stop=nobj+1, step=1):
        url = "https://kmtnet.kasi.re.kr/~ulens/event/" + year + \
                "/view.php?event=KMT-" + year + "-BLG-" + str(nn).zfill(4)
        response = urlopen(url)
        html = response.read()
        response.close()
        soup = BeautifulSoup(html,"html.parser")

        # Get the names of all the different lightcurve files (pysis names).
        links = soup.find_all('a', href=True)
        pysis_names = links[3].get_text(separator=',').split(',')[:-2]
        
        # Note, we are only keeping I-band lightcurves (V-band ones are not useful). 
        for pysis_name in pysis_names:
            if '_I.pysis' in pysis_name:
                # Grab the photometry for each alert's I-band lightcurve data into a pands dataframe.
                url = "https://kmtnet.kasi.re.kr/~ulens/event/" + year + "/data/KB" + \
                        year[2:] + str(nn).zfill(4) + "/pysis/" + pysis_name
                bytes_data = requests.get(url).content
                df = pd.read_csv(BytesIO(bytes_data), 
                                 delim_whitespace=True, skiprows=1, header=None, 
                                 names=['hjd', 'Delta_flux', 'flux_err', 'mag', 
                                        'mag_err', 'fwhm', 'sky', 'secz'])

                # Add columns for the alert name (of the form KBYYNNNN, YY=year, NNNN=alert number)
                # and telescope (lightcurve's pysis file.)
                df['alert_name'] = 'KB' + year[2:] + str(nn).zfill(4) 
                df['telescope'] = pysis_name

                # Write out the HJD, mag, mag_err, telescope, and alert_name data into the table.
                cols = ['hjd', 'mag', 'mag_err', 'telescope', 'alert_name']
                df[cols].to_sql(con=engine, schema=None, name="photometry", 
                                if_exists="append", index=False)
    t1 = time.time()             
    
    print('Took {0:.2f} seconds'.format(t1-t0))

def get_moa_params(alert_dir, year, nn):  
    """
    Get all the different MOA alert parameters (along with their uncertainties)
    from the individual web pages. The uncertainties are not listed on the 
    front summary page unfortunately.
    
    The "_e" values are the uncertainties.
    """
    # Add a column for the alert name (of the form MBYYNNN, YY=year, NNN=alert number)
    alert_name = 'MB' + year[2:] + str(nn + 1).zfill(3)  # need to make sure this always works.
    
    # Go to the alert page and scrape the data.
    url = "http://www.massey.ac.nz/~iabond/moa/alert" + year + "/" + alert_dir
    response = urlopen(url)
    html = response.read()
    response.close()
    soup = BeautifulSoup(html,"html.parser")

    # Parse the scraped data.
    meta = soup.find('div', id="metadata").text
    RA = meta.split('RA:')[1].split('Dec:')[0]
    Dec = meta.split('RA:')[1].split('Dec:')[1].split('Current')[0]

    tmax_str = soup.find('div', id="lastphot").text.split('<td>=<td align=right>')[1]
    tmax = moa_str_to_float(tmax_str.split()[1])
    tmax_e = moa_str_to_float(tmax_str.split('<td>')[2].split()[0])

    tE_str = soup.find('div', id="lastphot").text.split('<td>=<td align=right>')[2]
    tE = moa_str_to_float(tE_str.split()[0])
    tE_e = moa_str_to_float(tE_str.split('<td>')[2].split()[0])

    u0_str = soup.find('div', id="lastphot").text.split('<td>=<td align=right>')[3]
    u0 = moa_str_to_float(u0_str.split()[0])
    u0_e = moa_str_to_float(u0_str.split('<td>')[2].split()[0].split('<')[0])

    Ibase_str = soup.find('div', id="lastphot").text.split('<td>=<td align=right>')[4]
    Ibase = moa_str_to_float(Ibase_str.split()[0])
    Ibase_e = moa_str_to_float(Ibase_str.split('<td>')[2].split()[0].split('<')[0])

    assessment = soup.find('div', id="metadata").find_all('td', align='right')[4].text
    
    return alert_name, RA, Dec, tmax, tmax_e, tE, tE_e, \
            u0, u0_e, Ibase, Ibase_e, assessment, url
    
def get_moa_alerts(year):
    """
    Function that grabs all the different MOA alert parameters 
    (along with their uncertainties) for any given alert year, 
    and write them into a database.
    
    Parameters
    ----------
    year : int
        Year of the OGLE alerts you want.
        Valid choices are 2001 - 2019, inclusive.
        
    Outputs
    -------
    sqlite table called moa_alerts_<YYYY> in microlensing.db
    """    
    # Go to the MOA alerts site and scrape the page.
    year = str(year)
    url = "http://www.massey.ac.nz/~iabond/moa/alert" + year + "/alert.php"
    response = urlopen(url)
    html = response.read()
    response.close()
    soup = BeautifulSoup(html,"html.parser")

    # Get a list of all the bulge microlensing alert directories.
    links = soup.find_all('a', href=True)
    alert_dirs = []
    for ii, link in enumerate(links):
        if 'BLG' in link.text:
            alert_dirs.append(links[ii]['href'])
        
    # Figure out how many alerts there are in total
    npages = len(alert_dirs)

    # Go to the page for each bulge microlensing alert and scrape the parameters.
    # This process is parallelized using Pool (it's really slow to have to loop
    # over all pages, and this is an embarassingly parallel process.)
    _t0 = time.time()     
    num_workers = mp.cpu_count()  -2
    pool = mp.Pool(processes=num_workers)
    parallel_results = pool.starmap(get_moa_params, 
                                    zip(alert_dirs, repeat(year), range(npages)))
    _t1 = time.time()
    
    # Put it all into a dataframe and write out to the database.
    df = pd.DataFrame(parallel_results,
                     columns = ['alert_name', 'RA', 'Dec', 't0', 't0_err', 'tE', 'tE_err', 
                                'u0', 'u0_err', 'Ibase', 'Ibase_err', 'class', 'alert_url'])
    
    # Write HJD as HJD - 2450000 (less cumbersome digits)
    df['t0'] -= 2450000
    
    # Fill in the other columns
    df['Isrc'] = np.nan
    df['Isrc_err'] = np.nan
    df['srcfrac'] = np.nan
    df['srcfrac_err'] = np.nan
    
    df.to_sql(con=engine, schema=None, name="alerts", if_exists="append", index=False)
    
    _t1 = time.time()
    print('Read MOA alerts in {0:.2f} seconds'.format(_t1-_t0))


def calculate_srcfrac(mag_src, mag_base):
    """
    Calculate the source flux fraction srcfrac given
    source magnitude mag_src and baseline magnitude mag_base.
    """
    exp = -0.4 * (mag_src - mag_base)
    srcfrac = 10**exp
    
    return srcfrac

    
def get_ogle_params(year, nn, reg):  
    """
    Get all the different OGLE alert parameters (along with their uncertainties)
    from the individual web pages. The uncertainties are not listed on the 
    front summary page or lenses.par file in the data download unfortunately.
    
    The "_e" values are the uncertainties.
    """
    
    # Add a column for the alert name (of the form OBYYNNNN, YY=year, NNN=alert number)
    # Macy modification - use OB for BLG, OD for DG, and OG for GD
    alert_name = 'O' + reg[0] + year[2:] + str(nn + 1).zfill(4) 
    
    # Go to the alert page and scrape the data.
    url = "https://ogle.astrouw.edu.pl/ogle4/ews/" + year + \
            "/" + reg.lower() + "-" + str(nn+1).zfill(4) + ".html"
    response = urlopen(url)
    html = response.read()
    response.close()
    soup = BeautifulSoup(html,"html.parser")
    header_list = soup.find_all('table')[1].find('td').text.split()
    param_list = soup.find_all('table')[2].find('td').text.split()

    # Parse the scraped data.
    RA = header_list[7]
    Dec = header_list[10]
    Tmax = ogle_str_to_float(param_list, 1)
    Tmax_e = ogle_str_to_float(param_list, 3)
    tau =  ogle_str_to_float(param_list, 7)
    tau_e =  ogle_str_to_float(param_list, 9)
    Umin =  ogle_str_to_float(param_list, 11)
    Umin_e =  ogle_str_to_float(param_list, 13)
    fbl =  ogle_str_to_float(param_list, 23)
    fbl_e =  ogle_str_to_float(param_list, 25)
    Ibl =  ogle_str_to_float(param_list, 27)
    Ibl_e =  ogle_str_to_float(param_list, 29)
    I0 = ogle_str_to_float(param_list, 31)
    I0_e =  ogle_str_to_float(param_list, 33)

    return alert_name, RA, Dec, Tmax, Tmax_e, tau, tau_e, Umin, Umin_e, \
            fbl, fbl_e, Ibl, Ibl_e, I0, I0_e, url

def ogle_str_to_float(list_in, idx):
    """
    Little helper function to turn strings into floats.
    """
    try:
        return float(ne.evaluate(list_in[idx]))
    except:
        return np.nan
    
def moa_str_to_float(str_in):
    """
    Little helper function to turn strings into floats.
    """
    try:
        return float(ne.evaluate(str_in))
    except:
        return np.nan
        
def get_ogle_alerts(year):
    """
    Function that grabs all the different OGLE alert parameters 
    (along with their uncertainties) for any given alert year, 
    and write them into a database.
    
    Parameters
    ----------
    year : int
        Year of the OGLE alerts you want.
        Valid choices are 2001 - 2019, inclusive.
        
    Outputs
    -------
    sqlite table called ogle_alerts_<YYYY> in microlensing.db
    """
    # Go to the OGLE alerts site and scrape the page.
    year = str(year)
    url = "https://ogle.astrouw.edu.pl/ogle4/ews/" + year + "/ews.html"
    response = urlopen(url)
    html = response.read()
    response.close()
    soup = BeautifulSoup(html,"html.parser")

    # Figure out how many alert pages there are.
    tds = soup.find_all('td')[0::15]
    npages = len(tds)
    nblg = sum(list(map(lambda elem: int('BLG' in elem.get_text()), tds)))
    ndg = sum(list(map(lambda elem: int('DG' in elem.get_text()), tds)))
    ngd = sum(list(map(lambda elem: int('GD' in elem.get_text()), tds)))
    
    # Grab all the parameters from the OGLE pages.
    # Use pool to parallelize (it is very slow otherwise,
    # since we have to loop through each page individually.)
    _t0 = time.time()     
    num_workers = mp.cpu_count()  -2
    pool = mp.Pool(processes=num_workers)
    parallel_results_blg = pool.starmap(get_ogle_params, zip(repeat(year), range(nblg), repeat('BLG')))
    parallel_results_dg = pool.starmap(get_ogle_params, zip(repeat(year), range(ndg), repeat('DG')))
    parallel_results_gd = pool.starmap(get_ogle_params, zip(repeat(year), range(ngd), repeat('GD')))
    parallel_results = parallel_results_blg + parallel_results_dg + parallel_results_gd
    #print(parallel_results)
    _t1 = time.time() 

    # Put it all into a dataframe and write out to the database.
    df = pd.DataFrame(parallel_results,
                     columns =['alert_name', 'RA', 'Dec', 't0', 't0_err', 'tE', 'tE_err', 'u0', 'u0_err', 
                               'srcfrac', 'srcfrac_err', 'Ibase', 'Ibase_err', 'Isrc', 'Isrc_err', 'alert_url'])

    # Add in missing columns
    df['class'] = 'microlensing'
    
    # Write HJD as HJD - 2450000 (less cumbersome digits)
    df['t0'] -= 2450000

    df.to_sql(con=engine, schema=None, name="alerts", if_exists="append", index=False)
    
    print('Read OGLE alerts in {0:.2f} seconds'.format(_t1-_t0))
    
def get_kmtnet_alerts(year):
    """
    Function that grabs KMTNet alerts and writes the fit
    tE and Ibase parameters, as well as each alert's 
    classification, to a table in the database.
    
    Parameters
    ----------
    year : int
        Year of the KMTNet alerts you want.
        Valid choices are 2016 - 2022, inclusive.
        
    Outputs
    -------
    sqlite table called kmtnet_alerts_<YYYY> in microlensing.db
    Columns are alert_name, class, tE, Ibase, alert_url.
    """
    def kmtnet_str_to_float(item):
        try:
            return float(ne.evaluate(item.get_text().replace(u'\xa0', u'')))
        except:
            return
        
    # Go to the KMTNet alerts site and scrape the page.
    year = str(year)
    url = "https://kmtnet.kasi.re.kr/~ulens/event/" + year + "/"
    response = urlopen(url)
    html = response.read()
    response.close()
    soup = BeautifulSoup(html,"html.parser")

    # For some annoying reason, the KMTNet alerts system changes
    # across years randomly. Some years they report a single
    # classification for alerts, other years there are two 
    # classifications ("EF" and "AL", I don't know what it means).
    # For years where there are two classifications, I've picked 
    # AL classification arbitrarily.
    if year in ['2023','2022', '2020', '2017', '2016']:
        class_ = soup.find_all('td')[3::15][1:]
        RA = soup.find_all('td')[4::15][1:]
        Dec = soup.find_all('td')[5::15][1:]
        t_0 = soup.find_all('td')[6::15][1:]
        t_E = soup.find_all('td')[7::15][1:]
        u_0 = soup.find_all('td')[8::15][1:]
        Isource = soup.find_all('td')[9::15][1:]
        Ibase = soup.find_all('td')[10::15][1:]
    elif year in ['2021', '2019', '2018']:
        classEF = soup.find_all('td')[3::16][1:]
        classAL = soup.find_all('td')[4::16][1:]
        RA = soup.find_all('td')[5::16][1:]
        Dec = soup.find_all('td')[6::16][1:]
        t_0 = soup.find_all('td')[7::16][1:]
        t_E = soup.find_all('td')[8::16][1:]
        u_0 = soup.find_all('td')[9::16][1:]
        Isource = soup.find_all('td')[10::16][1:]
        Ibase = soup.find_all('td')[11::16][1:]
    else:
        raise Exception('Not a valid year')

    # Process output to get strings/floats as appropriate.
    RA_list = [item.get_text().replace(u'\xa0', u'') for item in RA]
    Dec_list = [item.get_text().replace(u'\xa0', u'') for item in Dec]
    t_0_list = [kmtnet_str_to_float(item) for item in t_0]
    t_E_list = [kmtnet_str_to_float(item) for item in t_E]
    u_0_list = [kmtnet_str_to_float(item) for item in u_0]
    Isource_list = [kmtnet_str_to_float(item) for item in Isource]
    Ibase_list = [kmtnet_str_to_float(item) for item in Ibase]
    if year in ['2023','2022', '2020', '2017', '2016']:
        class_list = [item.get_text().replace(u'\xa0', u'') for item in class_]
    elif year in ['2021', '2019', '2018']:
        classEF_list = [item.get_text().replace(u'\xa0', u'') for item in classEF]
        classAL_list = [item.get_text().replace(u'\xa0', u'') for item in classAL]

    # Get link to the alert page.
    if year in ['2023','2022', '2020', '2017', '2016']:
        alert_url = soup.find_all('td')[0::15][1:]
    elif year in ['2021', '2019', '2018']:
        alert_url = soup.find_all('td')[0::16][1:]
    else:
        raise Exception('Not a valid year')
    kmt_alert_url = 'https://kmtnet.kasi.re.kr/~ulens/event/' + year + '/'
    alert_url_list = [kmt_alert_url + item.find_all('a', href=True)[0]['href'] for item in alert_url]

    # Get alert name
    nn = len(t_E_list)
    alert_name = []
    for ii in np.arange(nn):
        alert_name.append('KB' + year[2:] + str(ii+1).zfill(4))

    if year in ['2023','2022', '2020', '2017', '2016']:
        # Put it all into a dataframe and write out to the database.
        df = pd.DataFrame(list(zip(alert_name, RA_list, Dec_list, t_0_list, t_E_list, u_0_list, 
                                   Isource_list, Ibase_list, class_list, alert_url_list)),
                         columns =['alert_name', 'RA', 'Dec', 't0', 'tE', 'u0',
                                   'Isrc', 'Ibase', 'class', 'alert_url'])
    elif year in ['2021', '2019', '2018']:
        df = pd.DataFrame(list(zip(alert_name, RA_list, Dec_list, t_0_list, t_E_list, u_0_list, 
                           Isource_list, Ibase_list, classEF_list, alert_url_list)),
                 columns =['alert_name', 'RA', 'Dec', 't0', 'tE', 'u0',
                           'Isrc', 'Ibase', 'class', 'alert_url'])
     
    df['t0_err'] = np.nan
    df['tE_err'] = np.nan
    df['u0_err'] = np.nan
    df['Ibase_err'] = np.nan
    df['Isrc_err'] = np.nan
    df['srcfrac'] = calculate_srcfrac(df['Isrc'], df['Ibase'])
    df['srcfrac_err'] = np.nan
    
    df.to_sql(con=engine, schema=None, name="alerts", if_exists="append", index=False)
