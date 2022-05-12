import query_alerts

query_alerts.get_moa_lightcurves(2022)
print('Downloaded MOA photometry from 2022 to database.')
query_alerts.get_ogle_lightcurves(2019)
print('Downloaded OGLE photometry from 2019 to database.')
query_alerts.get_kmtnet_lightcurves(2022)
print('Downloaded KMTNet photometry from 2022 to database.')

query_alerts.get_moa_alerts(2022)
print('Downloaded MOA alerts from 2022 to database.')
query_alerts.get_ogle_alerts(2019)
print('Downloaded OGLE alerts from 2019 to database.')
query_alerts.get_kmtnet_alerts(2022)
print('Downloaded KMTNet alerts from 2022 to database.')
