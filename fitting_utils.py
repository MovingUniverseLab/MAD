import json
from bs4 import BeautifulSoup
from urllib.request import urlopen
import time

def moa_lightcurves_from_list(list):
    #Download lightcurves from MOA website only for events in list, using code from query_alerts
    year = "20" + str(list[0][2:4])
    import warnings
    warnings.filterwarnings("ignore")
    url = "http://www.massey.ac.nz/~iabond/moa/alert" + year + "/alert.php"
    response = urlopen(url)
    html = response.read()
    response.close()
    soup = BeautifulSoup(html,"html.parser")

def kmt_lightcurves_from_list(list):
    #Download lightcurves from KMTNet website only for events in list, using code from query_alerts
    year = "20" + str(list[0][2:4])  
    url = "https://kmtnet.kasi.re.kr/~ulens/event/" + year + "/"
    response = urlopen(url)
    html = response.read()
    response.close()
    soup = BeautifulSoup(html,"html.parser")
    nobj = len(soup.find_all('td')[0::15][1:])
    t0 = time.time()

def ogle_lightcurves_from_list(list):
    #Download lightcurves from OGLE website only for events in list, using code from query_alerts
    year = "20" + str(list[0][2:4])
    ftp = ftplib.FTP("ftp.astrouw.edu.pl")
    ftp.login()
    ftp.cwd("ogle/ogle4/ews/" + year + "/")