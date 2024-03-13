# vampire
[English] | [中文](README-cn.md)
## UnitCellFile_2 

### How to install

After `git clone` this repo, using command `pip install -r requirements.txt` to install dependencie packages.

***

### How to use
+ From this website [materialPoject](https://materialsproject.org) to download `cif` file then add the `cif` file path in `cif_path`.
```python
cif_path = ' ' #cif file path
struct = get_struct(path=cif_path).get_from_local()#get cif file path from local
#or
struct = get_struct(name_id="",api="").get_from_MPR()#get cif file path with api

```

+ need to add parameters manually

```python
mat = [[1, 0, 0], 
		[0, 0, 0], 
		[0, 0, 0],
		[1, 0, 0], 
		[1, 0, 0],
		[1, 0, 0]]
#mat, corresponds to the atomic numbers in the mat file, format as follows

energy = [-10.88E-21,-14.62E-21, 4.18E-21]
#exchangeEnergy，corresponds to J_1,J_2,J_3

atom_type_num = 2
#atom_type，corresponds to the number of atom types in the package

interaction_type = 'isotropic' 
#isotropic and anisotropic

Dimension = 3 
#The dimension parameter can be modified, this program can generate 2D or 3D ucf files

save_path = '/Users/xbunax/Downloads/*.ucf'
#File save path

ucf = gen_ucf(struct, # file structure
              magmom, # magnetic moment of atoms
              suppercell, # supercell
              energy, # interaction
              atom_name, # lattice name
              save_path, # save path
              cif_path, # cif file path
              dimension, # dimension
              atom_type_num, # number of atom types in the lattice
              mat,
              interaction_type # type of interaction
              )

```

+ generate ucf file

```python
ucf = gen_ucf(struct, magmom, suppercell, energy, atom_name, save_path, cif_path, dimension, atom_type_num, mat, interaction_type)
# Get ucf object
ucf.write_ucf_custom(n) # n=1 writes only the nearest neighbor, n=2 writes the nearest and next-nearest neighbors, and so on
```
***

