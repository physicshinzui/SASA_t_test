#!/bin/bash
set -Ceu

python sasa.py -f *.pdb
head *dat
