idev
CPU_COUNT=$(grep -c ^processor /proc/cpuinfo)
export OMP_NUM_THREADS=$CPU_COUNT 

# ENV_TOML="env.toml"
# # Load our programs
# load_gaussian="stoml " + $ENV_TOML + " gaussian.load_gaussian"
# eval $load_gaussian
# load_lsqc="stoml " + $ENV_TOML + " lsqc.load_lsqc"
# eval $load_lsqc
# load_amber="stoml " + $ENV_TOML + " amber.load_amber"
# eval $load_amber 
# load_mopac="stoml " + $ENV_TOML + " mopac.load_mopac"

load_gaussian

pipenv shell
