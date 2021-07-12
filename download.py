#!/usr/bin/env python

import os
import time

import json
import requests


download_dir = './problems'


def download_prob(prob_id):
        api_token = os.environ['ICFP_2021_API_TOKEN']
        headers = { 'Authorization': 'Bearer ' + api_token }
        host = 'poses.live'
        file = f'/api/problems/{prob_id}'
        url = f'http://{host}{file}'

        try:
            response = requests.get(url, headers=headers)
        except HTTPError as http_err:
            print(f'HTTP error occurred" {http_err}')
        except Exception as err:
            print(f'Other error occurred" {err}')
        else:
            print(f'Problem #{prob_id}: Downloaded')

        json_content = response.json()
        with open(f'{download_dir}/{prob_id}.problem', 'w') as f:
            json.dump(json_content, f)  # indent=4
        print(f'Problem #{prob_id}: Saved')
        print('========================================')


def reformat(prob_id):
    fname = os.path.join('.', 'problems', f'{prob_id}.problem')

    with open(fname, 'r') as f:
        json_obj = json.load(f)

    with open(fname, 'w') as f:
        json.dump(json_obj, f)  # Rewrite with different format, as desired


if __name__ == '__main__':
    first = 1
    last = 132  # Previous values: 106, 88, 78, 59
    for prob_id in range(first, last + 1):
        download_prob(prob_id)
        time.sleep(0.5)
