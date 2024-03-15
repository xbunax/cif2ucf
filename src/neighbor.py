from structure import get_structure, init_structure
import numpy as np


class neighbor:
    def __init__(self, super_crystal, origin_crystal, magmom, supercell_matrix=np.array([3, 3, 3])):
        self.super_crystal = super_crystal
        self.origin_crystal = origin_crystal
        self.supercell_matrix = supercell_matrix
        self.magmom = magmom

    def get_relative_index_and_position(self, atom_index):
        len_of_crystal = len(self.origin_crystal)
        sum_crystal = self.supercell_matrix[0] * \
            self.supercell_matrix[1]*self.supercell_matrix[2]
        if atom_index > sum_crystal*len_of_crystal-1:
            print("atom index out of index")
            return 0
        center_position = np.array([1, 1, 1])
        idx = self.supercell_matrix[0]*self.supercell_matrix[1]*len_of_crystal
        idy = self.supercell_matrix[1]*len_of_crystal
        idz = len_of_crystal
        x = atom_index // idx
        y = int((atom_index - x*idx) / idy)
        z = (atom_index - x*idx - y*idy) / idz
        position = np.array([x, y, z], dtype=int) - center_position
        if (atom_index+1) % len(self.origin_crystal) == 0:
            return len_of_crystal-1, position
        return (atom_index+1) % len_of_crystal-1, position

    def get_center_crystal_index(self):
        origin_crystal_num_in_super_crystal = np.prod(self.supercell_matrix)
        len_of_origin_crystal = len(self.origin_crystal)
        center_index_start = int(
            origin_crystal_num_in_super_crystal/2)*len_of_origin_crystal
        center_index = [i for i in range(
            center_index_start, center_index_start+len_of_origin_crystal)]
        return center_index

    def index_values_to_dict(self, lst):
        result_dict = {}
        for index, value in enumerate(lst):
            if value in result_dict:
                result_dict[value].append(index)
            else:
                result_dict[value] = [index]
        return result_dict

    def get_neighbor_index(self, n):
        supercell, distance_matrix = init_structure(
            self.origin_crystal, self.magmom, self.supercell_matrix).get_supercell_and_distance()
        center_index = self.get_center_crystal_index()
        len_of_origin_crystal = len(self.origin_crystal)
        neighbor_index_list = np.empty(
            (len_of_origin_crystal, n), dtype=object)
        for i in center_index:
            sorted_distance_matrix = distance_matrix[i].copy()
            sorted_distance_matrix = np.unique(sorted_distance_matrix)
            sorted_distance_matrix.sort()
            neighbor_distance_list = sorted_distance_matrix[1:n+1]
            distance_dict = self.index_values_to_dict(distance_matrix[i])
            sorted_distance_dict = dict(
                sorted(distance_dict.items(), key=lambda item: item[0]))
            print(sorted_distance_dict)
            for j in range(len(neighbor_distance_list)):
                neighbor_index_list[i-center_index[0]
                                    ][j] = np.array(sorted_distance_dict[neighbor_distance_list[j]])
        return neighbor_index_list

    def get_neighbor_relative_index_and_position(self, n):
        neighbor_list = self.get_neighbor_index(n)
        row, col = neighbor_list.shape
        n = 0
        neighbor_index_position = []
        for i in range(col):
            neighbor_list_col = neighbor_list[:, i]
            for j in range(len(neighbor_list_col)):
                for k in neighbor_list_col[j]:
                    n += 1
                    index, position = self.get_relative_index_and_position(k)
                    neighbor_index_position.append([n, j, index, *position, i])
        # print(neighbor_index_position)
        return neighbor_index_position
