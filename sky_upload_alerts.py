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


# USE THE NEXT 5 LINES TO TEST IF THE API IS REACHABLE.
# response = api('GET', 'http://localhost:9000/api/sysinfo')

# print(f'HTTP code: {response.status_code}, {response.reason}')
# if response.status_code in (200, 400):
#     print(f'API is Readable:  JSON response: {response.json()}')
# else:
#     print(f' API is not readable. JSON response: {response.json()}')


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


def put_data_to_api(endpoint, data, api_token):
    """
    PUT data, updating an exisitng endpoint

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
        response = requests.put(url, headers=headers, json=data)
    if {endpoint} == {'sources'}:
        print("Uploading Sources")
        print("URL", url)
        print("Headers:", headers,)
        print("Data:", data)
        response = requests.put(url, headers=headers, json=data)
    else:
        # Make the POST request
        response = requests.put(url, headers=headers, json=data)

    # Check if the request was successful
    if response.status_code == 200:
        success = True
        print("Successful upload")
    else:
        success = False
        print("Failed  upload")

    print(response.json())                                    # Confirm output

    return response


'''
   Uploading multiple sources at once, such as the MOA, KLT or OGLE targets
   requires looping over all of the targets and posting individually
    
    TODO: Add in default programs, start with MOA, OGLE, KMTNET
        Turn this into a function that can be called for the three surveys.

    Note: Queries to targets that already exist are not altered, but the
        query still returns a success.    
    EX:
    data={
            "id": "14gqr_v3",
            "ra": 353.36647,
            "dec": 33.646149,
            "group_ids": [1],
        }
'''

endp = 'sources'
origin = "KMTNET"

df = pd.read_csv('data_alerts/kmtnet_alerts.csv',
                 usecols=['alert_name', 'RA', 'Dec', 'alert_url'])
df = df.iloc[20:25]  # TESTING
#  We want to incorpoarate these columns later: t0,tE,u0, alert_url, related_event, others?
ra_tmp = []
dec_tmp = []
for i, item in enumerate(df['RA'].values):
    # print(item,df['Dec'].iloc[i])
    coords = SkyCoord(ra=item, dec=df['Dec'].iloc[i], unit=(
        u.hour, u.deg), frame='icrs')
    # Annotations is a list of dictionaries.
    annotations = [
        df['alert_url'].iloc[i]
    ]
    data = {
        "id": df['alert_name'].iloc[i],
        "ra": coords.ra.value,
        "dec": coords.dec.value,
        "origin": origin,
        "summary": "Great target",
        "group_ids": [1, 2, 4],
        # "annotations": [str(df['alert_url'].iloc[i])],
        # "comments": ["test comment","test2"],
        # "altdata": {
        #     "gaia": {
        #         "info": {
        #             "Teff": 5770
        #         }
        #     }
        # }
    }

    print("Posting data", i)
    print("Data: ", data)
    # data = json.dumps(data)
    post_data_to_api(endp, data, token)
    print("Type of annotation: ", type(annotations))
    print("Annotations:", annotations)
print("Data type:", type(data))


# KMNNET TELESCOPE
#  https://noirlab.edu/public/programs/ctio/kmtnet-16m-telescope/
endp = 'telescope'
data = {
    'name':	'Korea Microlensing Telescope Network',
    'nickname': 'KMTNET Telescope',
    'diameter': 1.6,
    'fixed_location': True,
    'lat': -30.16717777,
    'lon': 70.804788888,
    'elevation': 2167
}
# Coords of Mount John Observatory in New Zealand
print("Lat in Degrees: ", -1*(30 + 10/60 + 1.84/3600))
print("Lon in Degres: ", 70+48/60 + 17.24/3600)
# post_data_to_api(endp,data,token) # POST oritinal
post_data_to_api(endp, data, token)  # PUT an update.
