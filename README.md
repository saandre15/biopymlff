# biopymlff(WIP)

## About biopymlff
A library for developing near constant scaling quantum mechanical quality Machine Learning Force Fields for large biochemical systems.

For more information checkout the [documentation](localhost).

## Dependencies
**External**
* lsqc
* amber20
* gaussian
* mopac
* psi4
* openbabel
* stoml

**Internal**
* See setup.py 

## Installation
```bash
git clone https://github.com/saandre15/biopymlff
cd biopymlff
./setup.py
```

## ML Model Snapshot Downloads

New biomolecule structures are downloaded from the PDB databank. Then the structures are fragmented and fed into the biopymlff program on Frontera at UTAustin and Titan at UTDallas. Everyday the models are snapshotted and publically updated on Github for public use.

### Proteins
| Name | Version | Status | Computational Dataset Size | Download | Experimental Dataset Size | Download |
|------| ------- | ---------------- | -------- | --------- | ------- | --------- |
| GEBF_GAP | v0.1 | WIP | 0 | Download | 0 | Download |
| GEBF_DP | v0.1 | WIP | 0 | Download | 0 | Download |
| GEBF_HDNN | v0.1 | WIP | 0 | Download | 0 | Download |


### Lipids
| Name | Version | Status | Computational Dataset Size | Download | Experimental Dataset Size | Download |
|------| ------- | ---------------- | -------- | --------- | ------- | --------- |
| GEBF_GAP | v0.1 | WIP | 0 | Download | 0 | Download |
| GEBF_DP | v0.1 | WIP | 0 | Download | 0 | Download |
| GEBF_HDNN | v0.1 | WIP | 0 | Download | 0 | Download |

### DNA/RNA
| Name | Version | Status | Computational Dataset Size | Download | Experimental Dataset Size | Download |
|------| ------- | ---------------- | -------- | --------- | ------- | --------- |
| GEBF_GAP | v0.1 | WIP | 0 | Download | 0 | Download |
| GEBF_DP | v0.1 | WIP | 0 | Download | 0 | Download |
| GEBF_HDNN | v0.1 | WIP | 0 | Download | 0 | Download |

### Water Model
| Name | Version | Status | Computational Dataset Size | Download | Experimental Dataset Size | Download |
|------| ------- | ---------------- | -------- | --------- | ------- | --------- |
| GEBF_GAP | v0.1 | WIP | 0 | Download | 0 | Download |
| GEBF_DP | v0.1 | WIP | 0 | Download | 0 | Download |
| GEBF_HDNN | v0.1 | WIP | 0 | Download | 0 | Download |


## Contact
* Torabifard Group(UT Dallas)

## Citation
See CITATION.cff

## References
* [1] Cheng, Zheng, et al. “Building Quantum Mechanics Quality Force Fields of Proteins with the Generalized Energy-Based Fragmentation Approach and Machine Learning.” Physical Chemistry Chemical Physics, 2021, https://doi.org/10.1039/d1cp03934b. 
* [2] Ko, Tsz Wai, et al. “A Fourth-Generation High-Dimensional Neural Network Potential with Accurate Electrostatics Including Non-Local Charge Transfer.” Nature Communications, vol. 12, no. 1, 2021, https://doi.org/10.1038/s41467-020-20427-2. 
* [3] Aric A. Hagberg, Daniel A. Schult and Pieter J. Swart, “Exploring network structure, dynamics, and function using NetworkX”, in Proceedings of the 7th Python in Science Conference (SciPy2008), Gäel Varoquaux, Travis Vaught, and Jarrod Millman (Eds), (Pasadena, CA USA), pp. 11–15, Aug 2008
* [4] Bartók, Albert P., and Gábor Csányi. “Gaussian Approximation Potentials: A Brief Tutorial Introduction.” International Journal of Quantum Chemistry, vol. 115, no. 16, 2015, pp. 1051–1057., https://doi.org/10.1002/qua.24927. 
* [5] Han, Jiequn, Linfeng Zhang, and Roberto Car. "Deep potential: A general representation of a many-body potential energy surface." arXiv preprint arXiv:1707.01478 (2017).
* [6] Kocer, Emir, Tsz Wai Ko, and Jörg Behler. "Neural Network Potentials: A Concise Overview of Methods." arXiv preprint arXiv:2107.03727 (2021).
* [7] Ko, Tsz Wai, et al. "A fourth-generation high-dimensional neural network potential with accurate electrostatics including non-local charge transfer." Nature communications 12.1 (2021): 1-11.
* [8] Salakhutdinov, Russ. “Lecture 2.” STA 4273H: Statistical Machine Learning. 20 Nov. 2021. 
* [9] Rugar, Daniel, and Paul Hansma. "Atomic force microscopy." Physics today 43.10 (1990): 23-30.
* [10] Mahoney, Michael W., and Petros Drineas. "CUR matrix decompositions for improved data analysis." Proceedings of the National Academy of Sciences 106.3 (2009): 697-702.