# vampire
[English](README.md) | [中文]
## UnitCellFile_2

### How to install

`git clone`本仓库后，安装python配置好环境后，在terminal输入`pip install -r requirements.txt` 安装运行所需的依赖库后即可。

***

### How to use
+ 从[materialPoject](https://materialsproject.org) 下载所需要的晶格文件后添加到`cif_path`路径中
```python
cif_path = '/Users/xbunax/Downloads/Mn2Au-2.cif'
struct = get_struct(path=cif_path).get_from_local()#从本地获取晶格结构
#或者
struct = get_struct(name_id="",api="").get_from_MPR()#利用api获取晶格结构

```

+ 需要手动添加参数

```python
mat = [[1, 0, 0], 
		[0, 0, 0], 
		[0, 0, 0],
		[1, 0, 0], 
		[1, 0, 0],
		[1, 0, 0]]
#mat，这一项是对应mat文件中的原子序数,格式如下

energy = [-10.88E-21,-14.62E-21, 4.18E-21]
#exchangeEnergy，这里对应J_1,J_2,J_3

atom_type_num = 2
#atom_type，这里对应原包中原子种类数

interaction_type = 'isotropic' 
#各向异性和各项同性

Dimension = 3 
#可以修改维度参数，本程序可以生成2维或者3维的ucf文件

save_path = '/Users/xbunax/Downloads/*.ucf'
#文件保存路径

ucf = gen_ucf(struct, #文件结构
			  magmom, #原子磁矩
			  suppercell,#扩包
			  energy,#相互作用
			  atom_name,#晶格名称
			  save_path,#保存路径
			  cif_path,#cif文件路径
			  dimension,#维度
			  atom_type_num,#晶格中原子种类个数
			  mat,
			  interaction_type#相互作用种类
			  )

```

+ 生成UCF

```python
ucf = gen_ucf(struct, magmom, suppercell, energy, atom_name, save_path, cif_path, dimension, atom_type_num,mat,interaction_type)
#获取ucf对象
ucf.write_ucf_custom(n)#n=1只写入最近邻，n=2写入最近邻和次最近邻，n=3同理
```
***


