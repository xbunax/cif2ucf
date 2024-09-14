from ase.io import read
import ase
import numpy as np
from ase.build import make_supercell


class get_structure:
    def __init__(self, cif_path: str):
        self.path = cif_path

    def from_file(self) -> ase.atoms.Atoms:
        return read(self.path)


class init_structure:
    def __init__(self, crystal: ase.atoms.Atoms, magmom: dict, super_matrix: list):
        self.crystal = crystal
        self.magmom = magmom
        self.super_matrix = super_matrix
        self.structure_without_magmom_atom = self.init_magmom()

    def normalize_cell(self, cell: ase.Atoms) -> ase.Atoms:
        """
        对给定的原子结构进行坐标归一化，并返回新的 Atoms 对象。

        参数:
        - cell (Atoms): ASE 中的原子结构。

        返回:
        - Atoms: 归一化坐标后的新的原子结构。
        """
        # 提取原子坐标
        positions = cell.get_positions()

        # 找到坐标的最大值和最小值
        min_coords = np.min(positions, axis=0)
        max_coords = np.max(positions, axis=0)

        # 避免除以零
        diff = max_coords - min_coords
        diff[diff == 0] = 1

        # 归一化坐标到 [0, 1] 范围
        normalized_positions = (positions - min_coords) / diff
        normalized_positions = np.around(normalized_positions, 3)

        # 创建一个新的 Atoms 对象，保留原有的元素信息
        normalized_cell = cell.copy()
        normalized_cell.set_positions(normalized_positions)

        return normalized_cell

    def suppercell_with_lattice_constant(self):
        lattice_constant = self.structure_without_magmom_atom.cell.lengths()
        print(lattice_constant)
        position = self.structure_without_magmom_atom.get_scaled_positions()
        position_with_lattice_parameter = position*lattice_constant
        self.structure_without_magmom_atom.set_positions(
            position_with_lattice_parameter)
        return self.structure_without_magmom_atom

    def init_magmom(self) -> ase.atoms.Atoms:
        magnetic_moments = []
        for atom in self.crystal:
            for i in self.magmom:
                if atom.symbol == i:
                    magnetic_moments.append(self.magmom[i])
        self.crystal.set_initial_magnetic_moments(magnetic_moments)
        indexes_to_remove = [i for i, magmom in enumerate(
            self.crystal.get_initial_magnetic_moments()) if magmom == 0.0]
        for index in sorted(indexes_to_remove, reverse=True):
            del self.crystal[index]
        structure_without_magmom = self.crystal
        # structure_without_magmom_normal = self.normalize_cell(
        #     structure_without_magmom)
        # structure_without_magmom_change = self.change_position_cell(
        #     structure_without_magmom_normal, [2**0.5/2, 2**0.5/2, 1])
        return structure_without_magmom

    def get_supercell_and_distance(self) -> [ase.atoms.Atoms, np.array]:
        supercell_matrix = np.diag(self.super_matrix)
        self.suppercell_with_lattice_constant()
        supercell = make_supercell(
            self.structure_without_magmom_atom, supercell_matrix)
        positions = supercell.get_positions()
        distance_matrix = np.sqrt(
            ((positions[:, np.newaxis, :] - positions[np.newaxis, :, :]) ** 2).sum(axis=2))
        distance_matrix = np.round(distance_matrix, 3)
        print(distance_matrix)
        return supercell, distance_matrix

    def change_position_cell(self, cell: ase.Atoms, change_tensor: np.array) -> ase.Atoms:
        position = np.array(cell.get_positions())
        position[:, 0] = position[:, 0]*change_tensor[0]
        position[:, 1] = position[:, 1]*change_tensor[1]
        position[:, 0], position[:, 2] = position[:, 2], position[:, 0]
        position = np.around(position, 4)
        change_cell = cell.copy()
        change_cell.set_positions(position)
        return change_cell
