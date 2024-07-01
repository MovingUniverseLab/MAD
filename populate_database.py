import query_alerts

def run(year=2024,sources=['OGLE','KMTNet', 'MOA']):
    if 'MOA' in sources:
        query_alerts.get_moa_alerts(year)
        print('Downloaded MOA alerts from '+str(year)+' to database.')
    if 'OGLE' in sources:
        df1 = query_alerts.get_ogle_alerts(year)
        print('Downloaded OGLE alerts from '+str(year)+' to database.')
    if 'KMTNet' in sources:
        query_alerts.get_kmtnet_alerts(year)
        print('Downloaded KMTNet alerts from '+str(year)+' to database.')

    # Get light curves
    #query_alerts.get_moa_lightcurves(2023)
    #print('Downloaded MOA photometry from 2023 to database.')
    #query_alerts.get_ogle_lightcurves(2023)
    #print('Downloaded OGLE photometry from 2023 to database.')
    #query_alerts.get_kmtnet_lightcurves(2023)
    #print('Downloaded KMTNet photometry from 2023 to database.')

if __name__ == '__main__':
    run()
