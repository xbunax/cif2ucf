from ase.build import make_supercell
import numpy as np
from writeucf import writeucf
from structure import get_structure, init_structure


if __name__ == "__main__":
    cif_path = ""
    ucf_path = ""
    atom_type = 2
    mat = [0, 1, 1, 0]
    hc = [0, 0, 0, 0]
    lc = [0, 0, 0, 0]
    dimension = 3
    magmom = {"Mn": 3.0, "Au": 0.0}
    super_matrix = [3, 3, 3]
    isotropic = "isotropic"
    n = 3
    exchange = [-10.88E-21, -14.62E-21, 4.18E-21]
    crystal = get_structure(cif_path).from_file()
    origin_crystal = init_structure(
        crystal, magmom, super_matrix).init_magmom()
    super_crystal = make_supercell(origin_crystal, np.diag(super_matrix))
    writeucf = writeucf(super_crystal, origin_crystal, ucf_path, magmom)
    writeucf.write_Unit_cell_size()
    writeucf.write_Unit_cell_Vector(dimension)
    writeucf.write_Atoms_num(atom_type, mat, lc, hc)
    writeucf.write_interactions(exchange, n, isotropic)
