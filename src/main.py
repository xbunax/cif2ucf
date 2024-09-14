# from ase.build import make_supercell
# from ase.visualize import view
# from ase.io import write
import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from writeucf import writeucf
from structure import get_structure, init_structure


def main(cif_path: str, ucf_path: str, atom_type: int, mat: list, hc: list, lc: list,
         dimension: int, magmom: dict, super_matrix: list, isotropic: str, n: int, exchange: list):
    crystal = get_structure(cif_path).from_file()
    origin_crystal = init_structure(
        crystal, magmom, super_matrix).init_magmom()
    # write('origin_crystal.cif')
    # super_crystal = make_supercell(origin_crystal, np.diag(super_matrix))
    write = writeucf(origin_crystal, ucf_path, magmom)
    write.write_Unit_cell_size()
    write.write_Unit_cell_Vector(dimension)
    write.write_Atoms_num(atom_type, mat, lc, hc)
    write.write_interactions(exchange, n, isotropic)


if __name__ == "__main__":
    cif_path = "../cif/Tm3Fe5O12.cif"
    ucf_path = "../demo/Tm3Fe5O12.ucf"
    atom_type = 1
    mat = [0 for i in range(40)]
    hc = [0 for i in range(40)]
    lc = [0 for i in range(40)]
    dimension = 3
    magmom = {"Fe": 3.0, "Tm": 0.0, "O": 0.0}
    super_matrix = [3, 3, 3]
    isotropic = "isotropic"
    n = 1
    exchange = [3.6e-22, 1.462e-20, -4.18e-21]
    main(cif_path,
         ucf_path,
         atom_type,
         mat, hc, lc,
         dimension,
         magmom,
         super_matrix,
         isotropic,
         n,
         exchange)
