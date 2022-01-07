# make -s -j10 -f NCCONCCOCCCNCCOCCCNCCOCCNCCNNCCONCCOCCCNCCOCOCNCCOCOCNCCOCCCNCCOCO.make
fname = NCCONCCOCCCNCCOCCCNCCOCCNCCNNCCONCCOCCCNCCOCOCNCCOCOCNCCOCCCNCCOCO
keyname = ${fname}.keys
files = 
all: $(files)
$(files):%.gebf:%.sxyz
	gauinput.py ${keyname} $(@:%.gebf=%.sxyz) g16 second
	g16 $(@:%.gebf=%.gjf) $(@:%.gebf=%.log)
	gebf.py $(@:%.gebf=%.chk) >& ${fname}.pr
