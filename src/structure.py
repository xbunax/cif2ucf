from ase.io import read
import numpy as np
from ase.build import make_supercell


class get_structure:
    def __init__(self, cif_path):
        self.path = cif_path

    def from_file(self):
        return read(self.path)


class init_structure:
    def __init__(self, structure, magmom, super_matrix):
        self.structure = structure
        self.magmom = magmom
        self.super_matrix = super_matrix
        self.structure_without_magmom_atom = self.init_magmom()

    def init_magmom(self):
        magnetic_moments = []
        for atom in self.structure:
            for i in self.magmom:
                if atom.symbol == i:
                    magnetic_moments.append(self.magmom[i])
        self.structure.set_initial_magnetic_moments(magnetic_moments)
        indexes_to_remove = [i for i, magmom in enumerate(
            self.structure.get_initial_magnetic_moments()) if magmom == 0.0]
        for index in sorted(indexes_to_remove, reverse=True):
            del self.structure[index]
        structure_without_magmom = self.structure
        return structure_without_magmom

    def get_supercell_and_distance(self):
        supercell_matrix = np.diag(self.super_matrix)
        supercell_matrix = np.diag(self.super_matrix)
        supercell = make_supercell(
            self.structure_without_magmom_atom, supercell_matrix)
        positions = supercell.get_positions()
        distance_matrix = np.sqrt(
            ((positions[:, np.newaxis, :] - positions[np.newaxis, :, :]) ** 2).sum(axis=2))
        distance_matrix = np.round(distance_matrix, 3)
        return supercell, distance_matrix


# class supercell:
#     def __init__(self, structure, ):
#         self.structure = structure
#         self.super_matrix = super_matrix
#
#     def get_supercell_and_distance(self):
#         supercell_matrix = np.diag(self.super_matrix)
#         supercell_matrix = np.diag(self.super_matrix)
#         supercell = make_supercell(self.structure, supercell_matrix)
#         positions = supercell.get_positions()
#         distance_matrix = np.sqrt(
#             ((positions[:, np.newaxis, :] - positions[np.newaxis, :, :]) ** 2).sum(axis=2))
#         distance_matrix = np.round(distance_matrix, 3)
#         return supercell, distance_matrix
#
#     def get_index_and_position(self, atom_index):
#         num_of_crystal = len(self.structure)
#         sum_crystal = self.super_matrix[0] * \
#             self.super_matrix[1]*self.super_matrix[2]
#         if atom_index > sum_crystal*num_of_crystal-1:
#             print("atom index out of index")
#             return 0
#         center_position = np.array([1, 1, 1])
#         idx = self.super_matrix[0]*self.super_matrix[1]*num_of_crystal
#         idy = self.super_matrix[1]*num_of_crystal
#         idz = num_of_crystal
#         x = atom_index // idx
#         y = int((atom_index - x*idx) / idy)
#         z = (atom_index - x*idx - y*idy) / idz
#         position = np.array([x, y, z], dtype=int) - center_position
#         if (atom_index+1) % len(self.structure) == 0:
#             return num_of_crystal-1, position
#         return (atom_index+1) % num_of_crystal-1, position
