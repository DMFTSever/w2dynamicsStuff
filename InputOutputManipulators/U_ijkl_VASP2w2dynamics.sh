#!/bin/bash

# VASP cRPA convention:
# W_abcd = < a,c | \hat{W} | b,d >
#
# w2dynamics convention:
# \hat{W} = U_ijkl c^\dagger_i c^\dagger_j c_l c_k
#
# => W_abcd =        < a,c | U_ijkl c^\dagger_i c^\dagger_j c_l c_k | b,d >
#    W_abcd = U_ijkl < a,c |        c^\dagger_i c^\dagger_j c_l c_k | b,d >
#    W_abcd = U_ijkl \delta_ai \delta_cj \delta_lb \delta_kd
#
# => W_abcd = U_acdb This line is implemented (W in U out)
# => U_ijkl = W_iljk

infile=${1}
outfile=${2}

awk 'NR==1{print} NR!=1{print $1, $3, $4, $2, $5, $6}' ${infile} > ${outfile}
