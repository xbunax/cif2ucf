from neighbor import neighbor
from structure import init_structure
import numpy as np
import ase


class writeucf:
    def __init__(self, origin_crystal: ase.atoms.Atoms, ucf_path: str, magmom: dict):
        init = init_structure(
            origin_crystal, magmom, super_matrix=[3, 3, 3])
        self.super_crystal, self.distance = init.get_supercell_and_distance()
        self.origin_crystal = origin_crystal
        self.ucf_path = ucf_path
        self.magmom = magmom

    def write_Unit_cell_size(self):
        unit_cell_size = self.origin_crystal.cell
        with open(self.ucf_path, 'a') as f:
            f.writelines("# Unit cell size:\n")
            for i in range(len(unit_cell_size)):
                f.write(str(unit_cell_size[i][i]) + "\t")
        return "unit_cell_size_write_Success"

    def write_Unit_cell_Vector(self, dimension=3):
        Unit_cell_Vector = np.eye(dimension)
        with open(self.ucf_path, "a") as f:
            f.writelines("\n# Unit cell vectors:\n")
            for i in range(len(Unit_cell_Vector)):
                for j in range(len(Unit_cell_Vector[i])):
                    f.write(str(Unit_cell_Vector[i][j]) + "\t")
                f.write("\n")
        return "write_Unit_cell_Vector_success"

    def write_Atoms_num(self, atom_type: int, mat: list, lc: list, hc: list):
        Atoms_num = len(self.origin_crystal)
        Atoms_position = self.origin_crystal.get_scaled_positions()
        with open(self.ucf_path, "a") as f:
            f.write("# Atoms num, id cx cy cz mat lc hc \n")
            f.write(str(Atoms_num)+"\t" + str(atom_type) + "\n")
            for i in range(len(Atoms_position)):
                f.write(f"{i}\t")
                for j in Atoms_position[i]:
                    f.write(f"{j}\t")
                f.write(f"{mat[i]}\t{lc[i]}\t{hc[i]}\n")
        return "write_Atoms_num_Success"

    def write_interactions(self, exchange: list, n: int, isotropic: str):
        exchange_list = neighbor(
            self.super_crystal, self.origin_crystal, self.magmom).get_neighbor_relative_index_and_position(n)
        print(exchange_list)
        with open(self.ucf_path, "a") as f:
            f.writelines(
                "#Interactions n exctype, id i j dx dy   dz        Jij\n")
            f.writelines(f"{len(exchange_list)}\t{isotropic}\n")
            for i in range(len(exchange_list)):
                for j in range(len(exchange_list[i])):
                    if j == len(exchange_list[i])-1:
                        f.write(f"{exchange[exchange_list[i][j]]}\t")
                    else:
                        f.write(f"{exchange_list[i][j]}\t")
                f.write("\n")
        return "write_interactions_success"
