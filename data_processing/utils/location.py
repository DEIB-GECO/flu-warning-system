import requests
import urllib
from time import sleep
import random 
from tqdm import tqdm
import io
from flu_dev_refact.utils import network
import pandas as pd
import numpy as np

class Geoapify:

    headers= {
        "Content-Type": "application/json",
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }
    # For north american boundary boxes used
    # https://gist.github.com/frankydp/d3a7488582f303874fe6
    # -178.2,6.6,-49.0,83.3
    # http://bboxfinder.com/#10.487812,-181.054688,84.928321,-8.789063
    # -181.054688,10.487812,-8.789063,84.928321
    continent_bounding_box = {
        'North America': '-178.2,6.6,-49,83.3'
    }

    continent_country_whitelist = {
        'North America': ['Anguilla', 'Antigua and Barbuda', 'Bahamas', 'Barbados', 'Belize', 'Bermuda', 'Canada', 'Cayman Islands', 'Costa Rica', 'Cuba', 'Dominica', 'Dominican Republic', 'Greenland', 'Grenada', 'Guadeloupe', 'Guatemala', 'Haiti', 'Honduras', 'Jamaica', 'Martinique', 'Mexico', 'Montserrat', 'Navassa Island', 'Nicaragua', 'Panama', 'Puerto Rico', 'Saba', 'Saint-Barthélemy', 'Saint Eustatius', 'Saint Kitts and Nevis', 'Saint Lucia', 'Saint Martin', 'Saint Pierre and Miquelon', 'Saint Vincent and the Grenadines', 'El Salvador', 'Turks and Caicos Islands', 'United States', 'British Virgin Islands', 'United States Virgin Islands']
    }
    
    @staticmethod
    def check_request_status(job_url, time_between_attemps: int = 10):
        sleep(time_between_attemps)
        job_status = requests.get(f"{job_url}&format=csv")
        job_status.raise_for_status()       # a pending job has status code 202 (doesn't raise in raise_for_status)
        if job_status.status_code != 200:    
            raise TimeoutError(f"Job not yet completed. Status code: {job_status.status_code}")
        else:
            # instead of returning immediately, do one more call
            sleep(1)
            job_status = requests.get(f"{job_url}&format=csv")
            job_status.raise_for_status()       
            return job_status.text
            
    @staticmethod
    def check_request_fail_callback(exc, job_url):
        print(f"Fallback method for call to: {job_url}")
        print(f"Attempts exceeded. Last attempt failed with exception {str(exc)}:")
        raise exc
    
    @staticmethod
    def post_batch_request_with_csv_output(url, body_json_serilaizable, preexisting_job_url=None):
        if not preexisting_job_url:
            # send job
            job_response = requests.post(url, json=body_json_serilaizable, headers=Geoapify.headers)
            job_response.raise_for_status()

            job_url = job_response.json()['url']
        else:
            job_url = preexisting_job_url
        print(f"Waiting completion of job at URL {job_url}")

        # wait for job completion
        csv_text_response = network.retry_f(Geoapify.check_request_status, job_url, fail_callback=Geoapify.check_request_fail_callback, max_attempts=6)  # check job result for 5  minutes

        # parse job output
        return pd.read_csv(io.StringIO(initial_value=csv_text_response))
    
    @staticmethod
    def geocoding_batch_request(api_key, bounding_box_coordinates=None):
        # bias_proximity="&proximity:-102.85429547459711,43.23212984712197"
        filter_bbox = f"&filter=rect:{bounding_box_coordinates}" if bounding_box_coordinates else ""
        return f"https://api.geoapify.com/v1/batch/geocode/search?&apiKey={api_key}{filter_bbox}"
        
    @staticmethod
    def geocoding_batch_request_body_from_locations(unique_locations: list):
        assert len(unique_locations) == len(set(unique_locations)), "Input locations not unique"
        # remove continents
        body = pd.Series(unique_locations).str.removeprefix("North America / ").str.removeprefix("Europe / ")   # ...
        # invert order from most specific location to less specific
        body = body.str.strip("").str.split(" / ").str[::-1].str.join(", ").tolist()   # list of unique places, each described by a string with commas separating country, county, city in reverse order (decreasing specificity from left to right)
        unique_locations2geoapify_input = dict(zip(unique_locations,body))
        assert len(body) < 1000, "Exceeded body limit (1000 addresses)"
        return body, unique_locations2geoapify_input
    
    @staticmethod
    def parse_country_or_usa_state_from_query_output(output: pd.DataFrame):
        output.country = output.country.str.replace("United States of America", "United States")
        unidentified_locations_mask = (pd.isna(output['result_type'])) | (output['result_type'] == 'unknown')    # unknown, amenity, building, street, suburb, district, postcode, city, county, state, country

        # detect locations recognised as country or USA state
        il_mask = ~unidentified_locations_mask
        il = output.loc[il_mask]
        countries_mask = (il['result_type'] == 'country')    # may include USA
        usa_country_mask = (il['country_code'] == 'us') & (il['result_type'] == 'country')               # need treatment                
        countries_mask_no_usa = countries_mask & (~usa_country_mask)   # countries except USA
        # states of USA
        usa_states_mask = (il['country'] == 'United States') & (il['result_type'] == 'state')

        new_output = output[['query.text', 'country', 'state', 'result_type']].copy()

        new_output['country_or_state'] = pd.NA
        new_output['country_or_state'] = np.where(countries_mask_no_usa, new_output.country, new_output['country_or_state'])
        new_output['country_or_state'] = np.where(usa_states_mask, new_output.state, new_output['country_or_state'])

        # locations at resolution higher than country or USA state
        higher_res_mask = (il_mask) & (~countries_mask_no_usa) & (~usa_states_mask) & (~usa_country_mask)

        # replace USA counties/cities with state
        higher_res_mask_within_usa = higher_res_mask & (il['country'] == 'United States')
        new_output['country_or_state'] = np.where(higher_res_mask_within_usa, new_output.state, new_output['country_or_state'])

        # replace non-USA locations at higher resolution with country 
        higher_res_mask_outside_usa = higher_res_mask & (il['country'] != 'United States')
        new_output['country_or_state'] = np.where(higher_res_mask_outside_usa, new_output.country, new_output['country_or_state'])

        return new_output     
    

class GeoKB:
    def __init__(self):
        self.countries: pd.DataFrame = None
        self.usa_states: pd.DataFrame = None

    def coordinates_countries(self) -> pd.DataFrame:
        if self.countries is None:
            countries = pd.read_csv("/Users/tom/Developer/flu-dev-refact/data/temporal_spatial_analysis/country_coordinates.csv", sep=",", header=0, index_col=0, dtype={'Lat': float, 'Lon': float})[['Lat','Lon']]
            countries.index.rename('location_name', inplace=True)
            countries.rename(columns={'Lat': 'lat', 'Lon': 'lon'}, inplace=True)
            self.countries = countries
        return self.countries
    
    def coordinates_usa_states(self) -> pd.DataFrame:
        if self.usa_states is None:
            usa_states = pd.read_csv("/Users/tom/Developer/flu-dev-refact/data/temporal_spatial_analysis/usa_states_coordinates.csv", sep=",", header=0, index_col=0, dtype={'Lat': float, 'Lon': float})
            usa_states.index.rename('location_name', inplace=True)
            usa_states.rename(columns={'Lat': 'lat', 'Lon': 'lon'}, inplace=True)
            self.usa_states = usa_states
        return self.usa_states