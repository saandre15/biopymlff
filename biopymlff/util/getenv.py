import toml
import os

def getenv():
    return toml.load(os.getcwd() + "/env.toml")