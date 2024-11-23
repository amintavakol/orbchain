[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)  
The code for the OrbChain model used in [AI for Interpretable Chemistry: Predicting Radical
Mechanistic Pathways via Contrastive Learning](https://proceedings.neurips.cc/paper_files/paper/2023/file/0ca70969597da7166128f7755c64ffd5-Paper-Conference.pdf) paper from NeurIPS 2023.

_Note: this implementation is based on [OEChem Toolkit](https://docs.eyesopen.com/toolkits/python/oechemtk/index.html). Make sure to obtain the license before using it._

## Overview  
**OrbChain** is a python library to model mechanistic chemial reactions.
It was introduce in [This paper](https://proceedings.neurips.cc/paper_files/paper/2023/file/0ca70969597da7166128f7755c64ffd5-Paper-Conference.pdf) to model mechanistic reactions for **chemical reaction prediction**.

### Key Features  
- It models a mechanistic reaction (i.e., a reaction with a single transition state) as the interaction between two Molecular Orbitals (MOs). 
- OrbChain(1): Given a mechanistic reaction as _reactants>>products arrows_ , it finds the interacting MOs.
- OrbChain(2): Given a set of reactant molecules, it will generate possible reactions in the form of _reactants>>products arrows_.

---

## Installation  
Follow these steps to set up the project:  

   ```bash
   git clone https://github.com/amintavakol/orbchain.git
   cd repository
   pip install -r requirements.txt
   ```

---

## Example Usage

OrbChain(1):
    `python find_reactive_pair.py <input_reaction>`


OrbChain(2):
    `python expand.py <input_reactants> <reaction_type>`

The `<reaction_type>` is either polar or radical.


---
## Citation

If you use this repository or find it helpful, please consider citing our work:
```bibtext
@article{tavakoli2024ai,
    title={AI for interpretable chemistry: predicting radical mechanistic pathways via contrastive learning},
    author={Tavakoli, Mohammadamin and Baldi, Pierre and Carlton, Ann Marie and Chiu, Yin Ting and Shmakov, Alexander and Van Vranken, David},
    journal={Advances in Neural Information Processing Systems},
    volume={36},
    year={2024}
    }
```
