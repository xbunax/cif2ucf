from ase.build import make_supercell
import numpy as np
from writeucf import writeucf
from structure import get_structure, init_structure


def main(cif_path: str, ucf_path: str, atom_type: int, mat: list, hc: list, lc: list,
         dimension: int, magmom: dict, super_matrix: list, isotropic: str, n: int, exchange: list):
    crystal = get_structure(cif_path).from_file()
    origin_crystal = init_structure(
        crystal, magmom, super_matrix).init_magmom()
    super_crystal = make_supercell(origin_crystal, np.diag(super_matrix))
    write = writeucf(super_crystal, origin_crystal, ucf_path, magmom)
    write.write_Unit_cell_size()
    write.write_Unit_cell_Vector(dimension)
    write.write_Atoms_num(atom_type, mat, lc, hc)
    write.write_interactions(exchange, n, isotropic)


if __name__ == "__main__":
    cif_path = "/Users/xbunax/Downloads/Mn3Pt.cif"
    ucf_path = "Mn3Pt.ucf"
    atom_type = 3
    mat = [0, 1, 1, 0]
    hc = [0, 0, 0, 0]
    lc = [0, 0, 0, 0]
    dimension = 3
    magmom = {"Mn": 3.0, "Pt": 0.0}
    super_matrix = [3, 3, 3]
    isotropic = "isotropic"
    n = 2
    exchange = [-3.60E-21, 2.20E-22]
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
