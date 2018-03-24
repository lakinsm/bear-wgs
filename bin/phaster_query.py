#!/usr/bin/env python3

import sys
import os
import requests
import glob
import time
import argparse
from urllib.request import urlretrieve
from pathlib import Path


SLEEP_TIME = 30  # seconds


phaster_post = 'http://phaster.ca/phaster_api?contigs=1'
phaster_get = 'http://phaster.ca/phaster_api?acc='


def submit_job(filepath):
    with open(filepath, 'rb') as f:
        res = requests.post(url=phaster_post, data=f)
        if not res.ok or 'error' in res.json():
            print(res.json())
            print(res.url)
            raise ValueError('Bad POST request: {}, {}\n'.format(res.status_code, res.reason))
        else:
            sys.stdout.write('POST received by server, url={}...\n'.format(res.url))
            print(res.json())
    return res


def check_status_retrieve(accession):
    res = requests.get(phaster_get + accession)
    if not res.ok:
        raise ValueError('Bad GET request: {}, {}\n'.format(res.status_code, res.reason))
    return res


def job_scheduler(infolder, outfolder):
    if outfolder[-1] == '/':
        modified_outpath = outfolder[:-1]
    else:
        modified_outpath = outfolder
    check_out = Path(modified_outpath)
    if not check_out.exists() or not check_out.is_dir():
        os.mkdir(outfolder)
    zip_check = Path(modified_outpath + '/zipfiles')
    if not zip_check.exists() or not zip_check.is_dir():
        os.mkdir(modified_outpath + '/zipfiles')
    for file in glob.glob(infolder + '/*'):
        samplename = file.split('/')[-1].split('.')[0]
        out_path = modified_outpath + '/' + samplename + '.phaster.out'
        in_file = Path(file)
        out_file = Path(out_path)
        if out_file.is_file():
            sys.stdout.write('Cached output found for sample {}...\n'.format(samplename))
            continue
        if not in_file.is_file():
            raise ValueError('The input file is not a valid file: {}'.format(file))
        sys.stdout.write('Scheduling sample {}...\n'.format(samplename))
        submit_obj = submit_job(file)
        time.sleep(SLEEP_TIME)
        current_status_obj = check_status_retrieve(submit_obj.json()['job_id'])
        current_data = current_status_obj.json()
        while 'error' in current_data.keys():
            sys.stdout.write('Job is pending on remote host for sample {}...\n'.format(samplename))
            time.sleep(SLEEP_TIME)
            current_status_obj = check_status_retrieve(submit_obj.json()['job_id'])
            current_data = current_status_obj.json()
        while current_data['status'] == 'Running...' or 'submissions ahead of yours' in current_data['status']:
            sys.stdout.write('Job is currently running for sample {}...\n'.format(samplename))
            time.sleep(SLEEP_TIME)
            current_status_obj = check_status_retrieve(submit_obj.json()['job_id'])
            current_data = current_status_obj.json()
        if current_data['status'] == 'Complete':
            sys.stdout.write('Job complete for sample {}\n'.format(samplename))
            zipfile_path = current_status_obj.json()['zip']
            res_url = current_status_obj.json()['url']
            summary = current_status_obj.json()['summary']
            with open(out_path, 'w') as out:
                out.write(summary)
            with open(modified_outpath + '/all_response_urls.txt', 'a') as out:
                out.write(res_url + '\n')
            urlretrieve('http://' + zipfile_path, modified_outpath + '/zipfiles/' + samplename + '_phaster.zip')


parser = argparse.ArgumentParser('phaster_query.py')
parser.add_argument('input', type=str, default=None, help='Path to input folder containing only the files to be queried, in multi-fasta format')
parser.add_argument('output', type=str, default=None, help='Path to output folder')


if __name__ == '__main__':
    args = parser.parse_args()
    job_scheduler(args.input, args.output)

