import query_alerts

# Get alerts
query_alerts.get_moa_alerts(2023)
print('Downloaded MOA alerts from 2023 to database.')
query_alerts.get_ogle_alerts(2023)
print('Downloaded OGLE alerts from 2023 to database.')
query_alerts.get_kmtnet_alerts(2023)
print('Downloaded KMTNet alerts from 2023 to database.')

# Get light curves
query_alerts.get_moa_lightcurves(2023)
print('Downloaded MOA photometry from 2023 to database.')
query_alerts.get_ogle_lightcurves(2023)
print('Downloaded OGLE photometry from 2023 to database.')
query_alerts.get_kmtnet_lightcurves(2023)
print('Downloaded KMTNet photometry from 2023 to database.')

