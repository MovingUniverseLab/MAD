# Upload the alerts to skyportal via the API

import pandas as pd


def main():
    # This is your 'main' function where you can call other functions or initialize multiprocessing tasks
    print("Main function is starting.")
    year = 2023
    # WORKING:
    get_moa_alerts(year)                #
    # get_kmtnet_alerts(year)           # No 2024 data yet
    # get_ogle_alerts(year)
    # get_moa_lightcurves(2023)
#    get_kmtnet_lightcurves(2024)
#    get_ogle_lightcurves(2024)


# Import the sky_query_alert function
#
# Use this, from the skyportal documentation as a starting point.
# Remember: update the token, make sure the container is running.

token = ''


def api(method, endpoint, data=None):
    headers = {'Authorization': f'token {token}'}
    response = requests.request(method, endpoint, json=data, headers=headers)
    return response


response = api('GET', 'http://localhost:5000/api/sysinfo')

print(f'HTTP code: {response.status_code}, {response.reason}')
if response.status_code in (200, 400):
    print(f'JSON response: {response.json()}')

if __name__ == '__main__':
    # Only run the main function when this module is executed as the main script
    main()
