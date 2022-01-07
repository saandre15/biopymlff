# make -s -j[N] -f 4znn_01_gebf_dft.make
fname = 4znn_01_gebf_dft
keyname = ${fname}.keys
files = 4znn_01_gebf_dft_1.gebf \
	4znn_01_gebf_dft_2.gebf \
	4znn_01_gebf_dft_3.gebf \
	4znn_01_gebf_dft_4.gebf \
	4znn_01_gebf_dft_5.gebf \
	4znn_01_gebf_dft_6.gebf \
	4znn_01_gebf_dft_7.gebf \
	4znn_01_gebf_dft_8.gebf \
	4znn_01_gebf_dft_9.gebf \
	4znn_01_gebf_dft_10.gebf \
	4znn_01_gebf_dft_11.gebf \
	4znn_01_gebf_dft_12.gebf \
	4znn_01_gebf_dft_13.gebf \
	4znn_01_gebf_dft_14.gebf \
	4znn_01_gebf_dft_15.gebf

all: $(files)
$(files):%.gebf:%.sxyz
	gauinput.py ${keyname} $(@:%.gebf=%.sxyz) g16 second
	g16 $(@:%.gebf=%.gjf) $(@:%.gebf=%.log)
	gebf.py $(@:%.gebf=%.chk) >& ${fname}.pr
