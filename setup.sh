#!/bin/bash

virtualenv venv
source venv/bin/activate
pip install tensorflow
pip install deepmd-kit
pip install ase
pip install quippy-ase
pip install mendeleev
pip install toml