#!/bin/bash
#PBS -N npxxy_dist_pickl
#PBS -e /home/willf/utils/log.err
#PBS -o /home/willf/utils/log.out
#PBS -l nodes=1
#PBS -l walltime=72:00:00
#PBS -m ea
#PBS -M willf@stanford.edu
#PBS -l mem=20000mb

python npxxy_dist_pickl.py
