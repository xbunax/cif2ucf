from ase.build import make_supercell
import numpy as np
from writeucf import writeucf
from structure import get_structure, init_structure


if __name__ == "__main__":
    cif_path = "/Users/xbunax/Downloads/mp-30409 Mn2Au.cif"
    ucf_path = "Mn2Au.ucf"
    atom_type=2
    mat = [0, 1, 1, 0]
    hc = [0, 0, 0, 0]
    lc = [0, 0, 0, 0]
    dimension = 3
    magmom = {"Mn": 3.0, "Au": 0.0}
    crystal = get_structure(cif_path).from_file()
    super_matrix = [3, 3, 3]
    exchange = [-10.88E-21, -14.62E-21, 4.18E-21]
    origin_crystal = init_structure(
        crystal, magmom, super_matrix).init_magmom()
    super_crystal = make_supercell(origin_crystal, np.diag(super_matrix))
    writeucf = writeucf(super_crystal, origin_crystal, ucf_path, magmom)
    writeucf.write_Unit_cell_size()
    writeucf.write_Unit_cell_Vector(dimension)
    writeucf.write_Atoms_num(atom_type,mat, lc, hc)
    writeucf.write_interactions(exchange, 3, isotropic="isotropic")


# crystal = read('/Users/xbunax/Downloads/Mn3Pt.cif')
# magmom = {"Mn": 3.0, "Pt": 0.0}
# crystal = init_crystal(crystal, magmom)
# super_matrix = [3, 3, 3]
# supercell, distance_matrix = get_supercell_and_position(crystal, super_matrix)
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# center_index = get_center_index(crystal, [3, 3, 3])
# for atom in supercell:
#     x, y, z = atom.position
#     if atom.index in center_index:
#         ax.scatter(x, y, z, color='r', s=100)
#         ax.text(x, y, z, atom.index)
#     else:
#         ax.scatter(x, y, z, color='b')
#         ax.text(x, y, z, atom.index)
# plt.show()
#
# # get_index_and_position(52,[3,3,3],atoms)
# for i in center_index:
#     # print(get_index_and_position(i,[3,3,3],atoms))
#     # print(distance_matrix[i])
#     distance_dict = index_values_to_dict(distance_matrix[i])
#     sorted_dict = dict(sorted(distance_dict.items(), key=lambda item: item[0]))
#     print(i, sorted_dict)
