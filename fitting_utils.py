import json
from bs4 import BeautifulSoup
from urllib.request import urlopen
import time
import requests
import numexpr as ne
import pandas as pd
from io import BytesIO
import numpy as np 
import ftplib
from pathlib import Path

def moa_lightcurves_from_list(moa_list):
    #Download lightcurves from MOA website only for events in list, using code from query_alerts
    year = "20" + str(moa_list[0][2:4])
    import warnings
    warnings.filterwarnings("ignore")
    url = "http://www.massey.ac.nz/~iabond/moa/alert" + year + "/alert.php"
    response = urlopen(url)
    html = response.read()
    response.close()
    soup = BeautifulSoup(html,"html.parser")
    links = soup.find_all('a', href=True)
    alert_dirs = []
    for ii, link in enumerate(links):
        if 'BLG' in link.text:
            alert_dirs.append(links[ii]['href'])
    file_dirs = {}
    for nn, alert_dir in enumerate(alert_dirs):
        # Scrape the page.
        alert_name = 'MB' + year[2:] + str(nn + 1).zfill(3)
        if alert_name in moa_list:
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
                             names=['mjd', 'delta_flux', 'flux_err', 'foo1', 'foo2', 'foo3', 'foo4', 'foo5'])

            # Add columns for magnitude and magnitude error, using the conversion
            # values we just figured out.
            df['mag'] = m - 2.5*np.log10(df['delta_flux'] + f)
            df['mag_err'] = 1.09 * df['flux_err']/(df['delta_flux'] + f)
            df['mjd'] -= 2400000.5

            # Get rid of all the nans which crop up during the conversion from delta flux to magnitude.
            df.dropna(axis='index', how='any', inplace=True)
            cols = ['mjd', 'mag', 'mag_err']
            
            #Download dataframe object as a csv file to the MOA specific folder within MAD
            path = 'lightcurves/MOA/' + alert_name + '.csv'
            file_dirs[alert_name] = path
            file_path = Path(path)
            file_path.parent.mkdir(parents=True, exist_ok=True)
            df[cols].to_csv(file_path, index=False)

    return file_dirs
    \




def kmt_lightcurves_from_list(kmt_list):
    #Download lightcurves from KMTNet website only for events in list, using code from query_alerts
    year = "20" + str(kmt_list[0][2:4])  
    url = "https://kmtnet.kasi.re.kr/~ulens/event/" + year + "/"
    response = urlopen(url)
    html = response.read()
    response.close()
    soup = BeautifulSoup(html,"html.parser")
    nobj = len(soup.find_all('td')[0::15][1:])
    file_dirs = {}
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
                alert_name = 'KB' + year[2:] + str(nn).zfill(4)
                if alert_name in kmt_list:
                    # Grab the photometry for each alert's I-band lightcurve data into a pands dataframe.
                    url = "https://kmtnet.kasi.re.kr/~ulens/event/" + year + "/data/KB" + \
                            year[2:] + str(nn).zfill(4) + "/pysis/" + pysis_name

                    bytes_data = requests.get(url).content
                    df = pd.read_csv(BytesIO(bytes_data), 
                                     delim_whitespace=True, skiprows=1, header=None, 
                                     names=['mjd', 'Delta_flux', 'flux_err', 'mag', 
                                            'mag_err', 'fwhm', 'sky', 'secz'])
                    df['mjd'] -= 2450000

                    # Write out the HJD, mag, mag_err, telescope, and alert_name data into the table.
                    cols = ['mjd', 'mag', 'mag_err']

                    #Download dataframe object as a csv file to the MOA specific folder within MAD
                    path = 'lightcurves/MOA/' + alert_name + '.csv'
                    file_dirs[alert_name] = path
                    file_path = Path(path)
                    file_path.parent.mkdir(parents=True, exist_ok=True)
                    df[cols].to_csv(file_path, index=False)
    return file_dirs

def ogle_lightcurves_from_list(ogle_list):
    #Download lightcurves from OGLE website only for events in list, using code from query_alerts
    year = "20" + str(ogle_list[0][2:4])
    ftp = ftplib.FTP("ftp.astrouw.edu.pl")
    ftp.login()
    ftp.cwd("ogle/ogle4/ews/" + year + "/")
    prefs = ['blg','dg','gd']
    nobjs = [int(sum(pref in x for x in ftp.nlst())/2) for pref in prefs]
    file_dirs = {}
    for i_pref, pref in enumerate(prefs):
        for nn in np.arange(start=1, stop=nobjs[i_pref]+1, step=1):
            alert_name = 'O' + pref[0].upper() + year[2:] + str(nn).zfill(4)
            if alert_name in ogle_list:
                # Grab the photometry for each alert.
                ftp.cwd(pref+"-" + str(nn).zfill(4))
                
                flo = BytesIO()
                ftp.retrbinary('RETR phot.dat', flo.write)
                flo.seek(0)
                df = pd.read_fwf(flo, header=0, 
                                 names=['mjd', 'mag', 'mag_err', 'see', 'sky'], 
                                 widths=[14, 7, 6, 5, 8])
                df['mjd'] -= 2400000.5
                cols = ['mjd', 'mag', 'mag_err']

                #Download dataframe object as a csv file to the MOA specific folder within MAD
                path = 'lightcurves/OGLE/' + alert_name + '.csv'
                file_dirs[alert_name] = path
                file_path = Path(path)
                file_path.parent.mkdir(parents=True, exist_ok=True)
                df[cols].to_csv(file_path, index=False)

                ftp.cwd("../")
    return file_dirs
