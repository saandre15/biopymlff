#!/bin/bash

python3 -m venv venv
source venv/bin/activate
# Read TOML File to get AMBERHOME

pip install --user pipenv
poetry installecho $e