# vampire
[English](README.md) | [中文]
## CIF2UCF 

### How to install

`git clone`本仓库后，安装python配置好环境后，在terminal输入`pip install -r requirements.txt` 安装运行所需的依赖库后即可。

***

### How to use
+ 在`cif2ucf.py`中配置好参数后直接在终端运行`./cif2ucf.py`

```python 
    cif_path = "" # cif文件路径 
    ucf_path = "" #ucf文件路径
    atom_type=2 # 晶胞中有多少原子种类
    mat = [0, 1, 1, 0] 
    hc = [0, 0, 0, 0]
    lc = [0, 0, 0, 0]
    dimension = 3 
    isotropic="isotropic" # 各向同性 or 各向异性
    super_matrix = [3, 3, 3] #扩包举证
    exchange = [-10.88E-21, -14.62E-21, 4.18E-21] #设置相互作用强度，分别为最近邻、次近邻、次次近邻
    n = 3 #设置到几级近邻例如n=1是只写入最近邻
    magmom = {"atom1": 3.0, "atom2": 0.0} #设置原子磁矩
    crystal = get_structure(cif_path).from_file()
    origin_crystal = init_structure(
        crystal, magmom, super_matrix).init_magmom()
    super_crystal = make_supercell(origin_crystal, np.diag(super_matrix))
    writeucf = writeucf(super_crystal, origin_crystal, ucf_path, magmom)
    writeucf.write_Unit_cell_size() 
    writeucf.write_Unit_cell_Vector(dimension)
    writeucf.write_Atoms_num(atom_type,mat, lc, hc)
    writeucf.write_interactions(exchange, n, isotropic)
```
### Tips
> + 本项目仍在非常早期阶段，还需要不断完善
> + 我只测试了用`Mn2Au.cif`生成`Mn2Au.ucf`计算curie point
> + 欢迎提issue，贡献代码
