'''
   Upload the alerts to skyportal via the API
   Use sky_upload_alerts.ipynb to develop this script.
'''

import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
import requests
import urllib
import json
pd.options.display.max_columns = None

# In order to test out the API, you either need to have the container running on 9000 (access on 5000)
#   the non-container version running on 8000.


# Add token and host
# Add a token (either from .tokens.yaml or from profile on skyportal)
token = 'a3dfebed-59cf-4144-b112-58e0d7a9b21f'
host = "http://localhost:9000"  # http://localhost:8000/api/sources
headers = {'Authorization': f'token {token}'}


# Run this function to test if the API is READABLE
def api(method, endpoint, data=None):
    headers = {'Authorization': f'token {token}'}
    response = requests.request(method, endpoint, json=data, headers=headers)
    return response


def post_data_to_api(endpoint, data, api_token):
    """
    Posts data to the specified API endpoint.

    Parameters:
    - endpoint: str, the endpoint of the API (e.g., 'telescope' or 'sources')
    - data: dict, the data to be posted to the API

    Returns:
    - response: requests.Response object, the response from the API
    """
    # Add token and host
    # api_token = ''  # Add a token (either from .tokens.yaml or from profile on skyportal)
#    host = "http://localhost:9000"  # http://localhost:8000/api/sources

    print("Attempting endpoint: ", endpoint)
    headers = {'Authorization': f'token {api_token}'}          # Define header
    url = urllib.parse.urljoin(
        host, f'api/{endpoint}')        # Input end point

    # Check to see if the source already exists, source is not overwritten if same name is used.
    # response = requests.get(url,headers=headers,json=data)
    # print("Previous response: ",response.json())

    print("Endpoint:", {endpoint})
    if {endpoint} == {'photometry'}:
        # Make the POST request
        print("Uploading photometry")
        response = requests.post(url, headers=headers, json=data)
    if {endpoint} == {'sources'}:
        print("Uploading Sources")
        print("URL", url)
        print("Headers:", headers,)
        print("Data:", data)
        response = requests.post(url, headers=headers, json=data)
    else:
        # Make the POST request
        print("Making OTHER post")
        response = requests.post(url, headers=headers, json=data)

    # Check if the request was successful
    if response.status_code == 200:
        success = True
        print("Successful upload")
    else:
        success = False
        print("Failed  upload")

    print(response.json())                                    # Confirm output

    return response


# SOURCE UPLOAD TEST
data = {
    "id": "14gqr_v5",
    "ra": 355.36647,
    "dec": 35.646149,
    "group_ids": [1],
}
endp = 'sources'
print("token:", token)
post_data_to_api(endp, data, token)  # Already in skyportal, but a good test.
