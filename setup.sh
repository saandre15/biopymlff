#!/bin/bash

virtualenv venv && source venv/bin/activate
pip install deepmd-kit
pip install ase
pip install quippy-ase
