# vampire
[English] | [中文](README-cn.md)
## CIF2UCF 

This project is related to [vampire](https://github.com/richard-evans/vampire) which can use `cif` file to generate `ucf` file.

### How to install

After `git clone` this repo, using command `pip install -r requirements.txt` to install dependencie packages.

***

### How to use
+ config in `cif2ucf.py` then run `./cif2ucf.py`
```python
    cif_path = "" # cif file path
    ucf_path = "" #generate ucf file path
    atom_type=2 # how many atom types in crystal
    mat = [0, 1, 1, 0] # Mn2Au demo
    hc = [0, 0, 0, 0]
    lc = [0, 0, 0, 0]
    dimension = 3 
    isotropic="isotropic" # whether isotropic or not
    super_matrix = [3, 3, 3] #make supercell matrix  I only test [3,3,3]
    exchange = [J_1 , J_2, J_3] #set exchange enerage
    n = 3 #set neighbor num if n=1 will only write the most neighbor atom 
    magmom = {"atom1": 3.0, "atom2": 0.0} #set atom magmon
```

### Tips

> + This project is still in very early stage.
> + I only test `Mn2Au.cif` to generate `Mn2Au.ucf` to calculate curie point, still need more tests.
> + Welcome to contribute



