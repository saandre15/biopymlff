from setuptools import setup, find_packages

setup(
    name="biopymlff",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "tensorflow",
        "deepmd-kit",
        "ase",
        "quippy-ase",
        "mendeleev",
        "toml",
        "networkx",
        "tk",
        "multipledispatch",
        "pysmiles",
        "dscribe",
        "pytraj"
    ]
)