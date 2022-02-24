#! /usr/bin/env bash
export GAUSS_SCRDIR=/scratch1/08249/sa748948/lsqc/24846154/
mkdir -p /scratch1/08249/sa748948/lsqc/24846154/
make -s -j10 -f 4znn_01_gebf_pm6.make
rm -rf /scratch1/08249/sa748948/lsqc/24846154
