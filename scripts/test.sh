source load_env.sh 

# Figure out why Langevin is not working
python -m unittest ../biopymlff/calculators/gebf_pm6_test.py

# NOTE: Run this when all modules are completely written
# python -m unittest discover -s "../biopymlff/" -p "*_test.py"