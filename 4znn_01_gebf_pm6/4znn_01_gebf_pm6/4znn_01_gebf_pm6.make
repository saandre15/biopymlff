# make -s -j10 -f 4znn_01_gebf_pm6.make
fname = 4znn_01_gebf_pm6
keyname = ${fname}.keys
files = 
all: $(files)
$(files):%.gebf:%.sxyz
	gauinput.py ${keyname} $(@:%.gebf=%.sxyz) g16 second
	g16 $(@:%.gebf=%.gjf) $(@:%.gebf=%.log)
	gebf.py $(@:%.gebf=%.chk) >& ${fname}.pr
