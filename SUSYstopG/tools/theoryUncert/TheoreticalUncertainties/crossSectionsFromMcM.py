#!/usr/bin/env python
import requests, cookielib
import subprocess
import argparse
import os

requests.packages.urllib3.disable_warnings()

parser = argparse.ArgumentParser()
parser.add_argument("--dataset_name", "-d", type=str, required=True)
args = parser.parse_args()

temp_cookie_file = "mcmcookies.temp"
mcm_address = "https://cms-pdmv.cern.ch/mcm"
subprocess.call(["cern-get-sso-cookie", "--krb", "-r", "-u", mcm_address, 
    "-o", temp_cookie_file])
c = cookielib.MozillaCookieJar(temp_cookie_file)
c.load()
search_options = {"db_name" : "requests", 
    "page" : -1,
    "dataset_name" : args.dataset_name}
r = requests.get("/".join([mcm_address, "search"]), params=search_options, cookies=c, verify=False)
for sample in r.json()["results"]:
    gen_params = sample["generator_parameters"]
    if len(gen_params) and "GS" in sample["prepid"]:
            # Sometimes generator_parameters has ?empty? lists in it
            values = [i for i in gen_params if type(i) is dict]
            print "\nCross section for dataset %s" \
                  " from request %s" % (sample["dataset_name"], sample["prepid"])
            print "---->  sigma = %s pb" % values[0]["cross_section"]
os.remove(temp_cookie_file)
