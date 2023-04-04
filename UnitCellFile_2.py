import pymatgen as mg
import numpy as np
from mp_api.client import MPRester
from collections import Counter
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class get_struct:

    def __init__(self, name_id=None, path=None, api=None):
        self.name_id = name_id
        self.path = path
        self.api = api

    def get_from_MPR(self):
        mpr = MPRester(self.api)
        struct = mpr.get_structure_by_material_id(self.name_id)
        return struct

    def get_from_local(self):
        struct = mg.core.Structure.from_file(self.path)
        return struct


class gen_ucf:

    def __init__(self, structure, magmom, supercell, energy, atom_name, save_path, cif_path,  dimension, atom_type_num,mat,interaction_type):
        self.structure = structure
        self.magmom = magmom
        self.supercell = supercell
        self.energy = energy
        self.atom_name = atom_name
        self.ucf_path = save_path
        self.cif_path = cif_path
        self.dimension=dimension
        self.atom_type_num=atom_type_num
        self.mat=mat
        self.interaction_type=interaction_type

    def get_struct_prim(self):
        struct_prim = get_struct(path=self.cif_path)
        self.set_magmon(struct_prim.get_from_local())
        return struct_prim.get_from_local()

    def get_magmon(self):
        return self.magmom

    def get_structpath(self):
        return self.cif_path

    def get_suppercell(self):
        return self.supercell

    def make_supercell(self, struct):
        struct.make_supercell(self.supercell)
        return struct

    def create_dimension(self):
        if self.dimension == 3:
            diagonal = np.array([1.0, 1.0, 1.0])
            vectors = np.diag(diagonal)
        else:
            diagonal = np.array([1.0, 1.0, 0.0])
            vectors = np.diag(diagonal)
        return vectors.tolist()

    def set_magmon(self, struct):
        for i in range(len(struct)):
            for j in self.get_magmon():
                if struct[i].specie.name == j:
                    struct[i].properties.update({"magmom": self.magmom[j]})
        for i in range(len(struct)):
            if struct[i].properties["magmom"] == 0:

                struct.remove_species([struct[i].specie.name])
                break

        return struct

    def get_magmom_num(self, struct):
        n = 0
        for i in range(len(struct)):
            if struct[i].properties["magmom"] != 0:
                n += 1
        return n

    def get_coord(self, struct):
        self.set_magmon(struct)
        coord = np.zeros([len(struct), 3])
        for i in range(len(struct)):
            coord[i] = struct[i].frac_coords
        return np.around(coord, 2)

    def get_distence(self,struct,center_index):
        self.set_magmon(struct)
        self.make_supercell(struct)
        distance_matrix=struct.distance_matrix
        distance_dict=[[] for i in range(len(center_index))]
        distance=[[] for i in range(len(center_index))]
        for i in range(len(center_index)):
            index = []
            atom_distance = []
            for j in range(len(distance_matrix[center_index[i]])):
                index.append(j)
                atom_distance.append(np.around(distance_matrix[center_index[i]][j],2))
            distance[i]=list(sorted(set(np.around(distance_matrix[i],2))))
            distance_dict[i]=dict(zip(index,atom_distance))
        print(Counter(distance_dict[0].values()))
        return distance_dict,distance

    def trans(self, arr, axix):
        n = len(arr) / self.supercell[axix]
        x = []
        for i in range(int(n), 2 * int(n)):
            x.append(list(arr)[i])
        return np.array(x)

    def get_center(self):
        vectorxyz = np.zeros([3, 1])
        struct = get_struct(path=self.cif_path).get_from_local()
        struct = self.make_supercell(struct)
        self.set_magmon(struct)
        coord_suppercell = np.around(struct.frac_coords,2)
        coord_prim = self.get_coord(self.get_struct_prim())
        xlib1 = sorted(set(coord_suppercell[:, 0]))
        ylib1 = sorted(set(coord_suppercell[:, 1]))
        zlib1 = sorted(set(coord_suppercell[:, 2]))
        xlib2 = np.array(sorted(set((1 / self.supercell[0]) * coord_prim[:, 0])))
        ylib2 = np.array(sorted(set((1 / self.supercell[1]) * coord_prim[:, 1])))
        zlib2 = np.array(sorted(set((1 / self.supercell[2]) * coord_prim[:, 2])))
        vectorxyz[0] = sum(self.trans(xlib1, 0) - np.around(xlib2, 2)) / len(self.trans(xlib1, 0) - np.around(xlib2, 2))
        vectorxyz[1] = sum(self.trans(ylib1, 1) - np.around(ylib2, 2)) / len(self.trans(ylib1, 1) - np.around(ylib2, 2))
        vectorxyz[2] = sum(self.trans(zlib1, 2) - np.around(zlib2, 2)) / len(self.trans(zlib1, 2) - np.around(zlib2, 2))
        vectorxyz = np.around(vectorxyz, 2)
        center_coord = np.zeros([len(coord_prim), 3])
        # print(np.shape(self.trans(xlib1, 0)))
        # center_coord[:,0],center_coord[:,1],center_coord[:,2]=self.trans(xlib1, 0).reshape(4,1),\
        #                                                       self.trans(ylib1, 1).reshape(4,1),\
        #                                                       self.trans(zlib1, 2).reshape(4,1)
        center_coord[:, 0], center_coord[:, 1], center_coord[:, 2] = (1 / self.supercell[0]) * coord_prim[:, 0] + \
                                                                     vectorxyz[0], \
                                                                     (1 / self.supercell[1]) * coord_prim[:, 1] + \
                                                                     vectorxyz[1], \
                                                                     (1 / self.supercell[2]) * coord_prim[:, 2] + \
                                                                     vectorxyz[2]
        center_coord = np.around(center_coord, 2)
        vectorxyz = np.around(vectorxyz, 2)
        center_index=[]
        for i in center_coord:
            for j in range(len(struct)):
                coord=coord_suppercell[j]
                if np.allclose(i,coord,rtol=0.01,atol=0.1):
                    center_index.append(j)
        return center_coord, vectorxyz, center_index

    def get_index(self, list):
        move_tensor = np.array(
            [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [1, 1, 0], [1, -1, 0], [-1, -1, 0], [-1, 1, 0],
             [0, 0, 0],
             [1, 0, 1], [0, 1, 1], [-1, 0, 1], [0, -1, 1], [1, 1, 1], [1, -1, 1], [-1, -1, 1], [-1, 1, 1],
             [0, 0, 1],
             [1, 0, -1], [0, 1, -1], [-1, 0, -1], [0, -1, -1], [1, 1, -1], [1, -1, -1], [-1, -1, -1], [-1, 1, -1],
             [0, 0, -1]])
        center_coord, vectorxyz, center_index = self.get_center()
        for i in move_tensor:
            for j in range(len(center_coord)):
                k=np.around(i.reshape(1, 3) * vectorxyz.reshape(1, 3) + list, 2)
                if np.allclose(k,center_coord[j],rtol=0.01,atol=0.1):
                    return j, np.negative(i)

    def get_key(self, dict, value):
        return [k for k, v in dict.items() if v == value]

    def distance(self):
        center_coord, vector, center_index = self.get_center()
        struct1 = get_struct(path=self.cif_path).get_from_local()
        self.set_magmon(struct1)
        distance_dict,distance=self.get_distence(struct1,center_index)
        neighbournnii = [[] for i in range(len(center_coord))]
        neighbournnnii = [[] for i in range(len(center_coord))]
        neighbournnnnii = [[] for i in range(len(center_coord))]
        for i in range(len(center_index)):
            neighbournnii[i]=self.get_key(distance_dict[i],distance[i][1])
            neighbournnnii[i]=self.get_key(distance_dict[i],distance[i][2])
            neighbournnnnii[i]=self.get_key(distance_dict[i],distance[i][3])
        return neighbournnii, neighbournnnii, neighbournnnnii

    def get_neighbour_coord(self):
        struct2 = get_struct(path=self.cif_path).get_from_local()
        coords = self.get_coord(self.make_supercell(struct2))
        neighbournnii, neighbournnnii, neighbournnnnii = self.distance()
        neighbournn_coord = [[] for i in range(len(self.get_coord(self.get_struct_prim())))]
        neighbournnn_coord = [[] for i in range(len(self.get_coord(self.get_struct_prim())))]
        neighbournnnn_coord = [[] for i in range(len(self.get_coord(self.get_struct_prim())))]
        for i in range(len(neighbournnii)):
            for j in neighbournnii[i]:
                neighbournn_coord[i].append(coords[j].tolist())
            for j in neighbournnnii[i]:
                neighbournnn_coord[i].append(coords[j].tolist())
            for j in neighbournnnnii[i]:
                neighbournnnn_coord[i].append(coords[j].tolist())
        neighbournn_coord = np.array(neighbournn_coord)
        neighbournnn_coord = np.array(neighbournnn_coord)
        neighbournnnn_coord = np.array(neighbournnnn_coord)
        return neighbournn_coord, neighbournnn_coord, neighbournnnn_coord

    def get_neighbour(self):
        neighbournn_coord, neighbournnn_coord, neighbournnnn_coord = self.get_neighbour_coord()
        neighbournn = [[] for i in range(len(self.get_coord(self.get_struct_prim())))]
        neighbournnn = [[] for i in range(len(self.get_coord(self.get_struct_prim())))]
        neighbournnnn = [[] for i in range(len(self.get_coord(self.get_struct_prim())))]
        for i in range(len(neighbournn_coord)):
            for j in neighbournn_coord[i]:
                target_atom, neighbournn_direction = self.get_index(j)
                neighbournn_direction="\t".join(str(h) for h in neighbournn_direction)
                k = "{}\t{}\t{}\t{}".format(i, target_atom, neighbournn_direction, self.energy[0])
                neighbournn[i].append(k)
        for k in range(len(neighbournnn_coord)):
            for l in neighbournnn_coord[k]:
                target_atom, neighbournnn_direction = self.get_index(l)
                neighbournnn_direction="\t".join(str(i) for i in neighbournnn_direction)
                p = "{}\t{}\t{}\t{}".format(k, target_atom, neighbournnn_direction, self.energy[1])
                neighbournnn[k].append(p)
        for m in range(len(neighbournnnn_coord)):
            for n in neighbournnnn_coord[m]:
                target_atom, neighbournnnn_direction = self.get_index(n)
                neighbournnnn_direction="\t".join(str(s) for s in neighbournnnn_direction)
                k = "{}\t{}\t{}\t{}".format(m, target_atom, neighbournnnn_direction, self.energy[2])
                neighbournnnn[m].append(k)
        return neighbournn, neighbournnn, neighbournnnn

    def write_atom_type(self):
        atom_num=[[] for i in range(len(self.get_coord(self.get_struct_prim())))]
        for i in range(len(self.mat)):
            str_mat="\t".join(str(i) for i in self.mat[i])
            struct_prim=self.get_struct_prim()
            str_frac_coord="\t".join(str(i) for i in np.around(struct_prim.frac_coords[i],2))
            atom_num[i]="{}\t{}\t{}".format(i,str_frac_coord,str_mat)
        with open(self.ucf_path, "a+") as f:
            f.writelines("# Atoms num, id cx cy cz mat lc hc"+'\n')
            f.write(str(len(self.get_coord(self.get_struct_prim())))+'\t')
            f.write(str(self.atom_type_num)+'\n')
            for i in atom_num:
                f.write(str(i))
                f.write('\n')
        return True

    def write_unit_cell_size(self):
        struct3=self.get_struct_prim()
        with open(self.ucf_path, 'a+') as f:
            f.writelines('# Unit cell size:' + '\n')
            for i in struct3.lattice.abc:
                f.write(str(i) + '\t')
            f.write('\n')
        print('unit_cell_size write success')
        return True

    def write_Unit_cell_vectors(self):
        with open(self.ucf_path, 'a+') as f:
            f.writelines('# Unit cell vectors:' + '\n')
            for i in range(len(self.create_dimension())):
                for j in self.create_dimension()[i]:
                    f.write(str(j) + '\t')
                f.write('\n')
        print('Unit_cell_vectors success')
        return True

    def write_interaction(self):
        self.creat_atom_type()
        neighbournn, neighbournnn, neighbournnnn = self.get_neighbour()
        neighbour = neighbournn + neighbournnn + neighbournnnn
        sum_interaction=0
        for i in neighbour:
            sum_interaction+=len(i)
        n=0
        with open(self.save_path, "a+") as f:
            f.writelines('#Interactions n exctype, id i j dx dy   dz        Jij' + '\n')
            f.write(str(sum_interaction)+'\t')
            f.write(self.interaction_type+'\n')
            for i in range(len(neighbour)):
                for j in neighbour[i]:
                    n+=1
                    f.write(str(n)+'\t')
                    f.write(j)
                    f.write("\n")
        print('interaction write success')
        print('ucf creat success')

    def write_ucf(self):
        self.write_unit_cell_size()
        self.write_Unit_cell_vectors()
        self.write_atom_type()
        self.write_interaction()

    def write_ucf_custom(self,n):
        self.write_unit_cell_size()
        self.write_Unit_cell_vectors()
        self.write_atom_type()
        neighbournn, neighbournnn, neighbournnnn = self.get_neighbour()
        if n==1:
            neighbour = neighbournn
        elif n==2:
            neighbour = neighbournn+neighbournnn
        elif n==3:
            neighbour = neighbournn + neighbournnn+neighbournnnn
        sum_interaction=0

        for i in neighbour:
            sum_interaction+=len(i)
        n=0
        with open(self.ucf_path, "a+") as f:
            f.writelines('#Interactions n exctype, id i j dx dy   dz        Jij' + '\n')
            f.write(str(sum_interaction)+'\t')
            f.write(self.interaction_type+'\n')
            for i in range(len(neighbour)):
                for j in neighbour[i]:
                    n+=1
                    f.write(str(n)+'\t')
                    f.write(j)
                    f.write("\n")
        print('interaction write success')
        print('ucf creat success')

    def plot_figure(self):
        struct=get_struct(path=self.cif_path).get_from_local()
        fig=plt.figure()
        ax = Axes3D(fig)
        self.make_supercell(struct)
        self.set_magmon(struct)
        center_coord,vectorxyz,center_index = self.get_center()
        neighbournn, neighbournnn, neighbournnnn = self.get_neighbour_coord()
        ax.scatter(struct.frac_coords[:,0],struct.frac_coords[:,1],struct.frac_coords[:,2],c='g')
        ax.scatter(center_coord[:,0],center_coord[:,1],center_coord[:,2],c='r',marker='+',s=50)
        ax.scatter(neighbournn[0][:,0],neighbournn[0][:,1],neighbournn[0][:,2],c='b',marker='o',s=30)
        ax.scatter(neighbournn[1][:,0],neighbournn[1][:,1],neighbournn[1][:,2],c='b',marker='*',s=30)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.show()


def main():
    cif_path = "/Users/xbunax/Downloads/CrI3.cif"
    structure = get_struct(path=cif_path)
    struct = structure.get_from_local()
    energy = [-2.72e-22, -1.6e-24,  1.92e-23]
    suppercell = [3, 3, 3]
    magmom = {"Cr": 3.1, "I": 0}
    atom_name = 'CrI3'
    atom_type_num=2
    interaction_type='isotropic'
    save_path = "/Users/xbunax/Downloads/CrI3.ucf"
    dimension=3
    mat=[[0,0,0],
         [1,0,0]]
    ucf = gen_ucf(struct, magmom, suppercell, energy, atom_name, save_path, cif_path, dimension, atom_type_num,mat,interaction_type)
    neighbournn,neighbournnn,neighbournnnn=ucf.get_neighbour_coord()
    print('neighbournn:',neighbournn,'\n',
          'neighbournnn:',neighbournnn,'\n',
          'neighbournnnn',neighbournnnn)
    # ucf.write_ucf_nn()
    # ucf.plot_figure()
    ucf.write_ucf_custom(3)
    # print(ucf.set_magmon(struct))


if __name__ == '__main__':
    main()
