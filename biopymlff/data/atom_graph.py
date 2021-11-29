import os
import random
import math

from multipledispatch import dispatch

from ase.atoms import Atoms
from ase.atom import Atom

from ase.geometry.analysis import Analysis

import networkx as nx

from typing import Callable

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

# matplotlib.use('TkAgg')

import numpy as np

from biopymlff.data.atom_graph_node import AtomGraphNode


class AtomGraph():
    
    counter = 1

    electrons = {
        'H': 1,
        'C': 6,
        'N': 7,
        'O': 8,
        'F': 9,
    }

    # @dispatch(Atoms)
    def __init__(self, atoms: Atoms):
        self.atoms = atoms
        self.graph = self.to_graph(atoms)

    # @dispatch(str)
    # def __init__(self, file: str):
    #     atom: Atoms = None
    #     bond_list = None
    #     if ".mol" in file or ".mol2" in file:
    #         try:
    #             with open(file) as file:
    #                 lines = file.readlines()
    #                 atom_mode = False
    #                 bonding_mode = False
                    
    #                 atom_list = []
    #                 bond_list = []

    #                 for line in lines:
    #                     if "@<TRIPOS>ATOM" in line: atom_mode = True; bonding_mode = False ; continue
    #                     if "@<TRIPOS>BOND" in line: bonding_mode =  True; atom_mode = False ; continue
    #                     if atom_mode:
    #                         vals = line.split()
    #                         symbol = vals[5].split(".")[0]
    #                         print(symbol)
    #                         x = float(vals[2])
    #                         y = float(vals[3])
    #                         z = float(vals[4])
    #                         charge = float(vals[8])
    #                         atom  = Atom(symbol=symbol, position=(x, y, z), charge=charge)
    #                         atom_list.append(atom)
    #                     if bonding_mode:
    #                         vals = line.split()
    #                         idx = int(vals[1]) - 1
    #                         idy = int(vals[2]) - 1
    #                         bond_type = vals[3]
    #                         if bond_type == "am": bond_type = 4
    #                         if bond_type == "ar": bond_type = 5
    #                         if bond_type == "du": bond_type = 6
    #                         if bond_type == "un": bond_type = 7
    #                         if bond_type == "nc": bond_type = 8
    #                         bond_type = int(bond_type)
    #                         bond_list.append((idx, idy, bond_type))
                    
    #                 atom = Atoms(atom_list)

    #         except IOError: print("Failed to read " + file)
    #     else: raise IOError("Make sure the extension type is parsable ")

    #     print(atom)
    #     print(bond_list)

    #     self.__init__(atom, bond_list)

    def reset(self):

        graph_list = [node for node in self.graph.nodes]
        
        for node in graph_list:
            node.setVisited(False)
        

    def show(self):
        node_xyz = []
        edge_xyz = []
        it = self.graph.nodes
        for val in it:
            # print("ATOM")
            # print(val.getAtom())
            node_xyz.append([val.getAtom().x, val.getAtom().y, val.getAtom().z])
            neighbors = self.graph.neighbors(val)
            for adjacent in neighbors:
                # print("ADJACENT")
                # print(adjacent.getAtom())
                edge_xyz.append([[val.getAtom().x, val.getAtom().y, val.getAtom().z], [adjacent.getAtom().x, adjacent.getAtom().y, adjacent.getAtom().z]])

        fig: Figure = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        np_node_xyz=np.array(node_xyz)
        np_edge_xyz=np.array(edge_xyz)

        # print(node_xyz)
        # print(edge_xyz)
        
        ax.scatter(*np_node_xyz.T, s=100, ec="w")
        
        for edge in np_edge_xyz:
            ax.plot(*edge.T, color="tab:gray")
        
        def _format_axes(ax):
            """Visualization options for the 3D axes."""
            # Turn gridlines off
            ax.grid(False)
            # Suppress tick labels
            for dim in (ax.xaxis, ax.yaxis, ax.zaxis):
                dim.set_ticks([])
            # Set axes labels
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")


        _format_axes(ax)
        fig.tight_layout()
        # plt.savefig(os.getcwd() + "/test.png")
        plt.show()

    def fragment_by_bond_as_atoms_list(self, symbolA: str, symbolB: str) -> list:
        self.reset()
        fragments = self.fragments_by_bond_as_indexes(symbolA, symbolB)
        graph_list = [node for node in self.graph.nodes]
        atoms_list = []
        print(fragments)
        for fragment in fragments:
            atom_list = []
            for index in fragment:
                print(index)
                atom: Atom = graph_list[index].getAtom()
                atom_list.append(atom)
            atoms: Atoms = Atoms(atom_list)
            atoms_list.append(atoms)
            print(atoms.get_positions())
        return atoms_list
                    
    # Returns a list of atoms
    def fragments_by_bond_as_indexes(self, symbolA: str, symbolB: str) -> list:
        self.reset()
        all_atoms = [atom for atom in self.atoms]
        graph_list = [node for node in self.graph.nodes]
        not_explored_atoms = [index for index in range(0, len(all_atoms))]
        node_fragments = []
        fragments = []
        system_size = len(not_explored_atoms)

        def traversal_fn(cur: AtomGraphNode, next: AtomGraphNode):
            curSymbol = cur.getAtom().symbol
            nextSymbol = next.getAtom().symbol
            node_fragments.append(cur)
            return not ((curSymbol == symbolA and nextSymbol == symbolB) or (curSymbol == symbolB and nextSymbol == symbolA))

        # While there is still atoms that are not explored
        while len(not_explored_atoms) > 0:
            # Select a random atom within the not explored
            index = math.floor(random.random() * len(not_explored_atoms))
            # Traverse and add the atom to the current fragment
            self.traverse(graph_list[not_explored_atoms[index]], traversal_fn)
            node_fragments = list(set(node_fragments))
            not_explored = []
            # Resets the current fragment
            fragment = []
            print("Current Fragment Length " + str(len(node_fragments)))
            print(not_explored_atoms)
            for index in not_explored_atoms:
                has_appended = False
                
                for node in node_fragments:
                    # If the node in the graph list matches the node in the current fragment then add it to the explored atom
                    if graph_list[index].getAtom() == node.getAtom():
                        print("atom append to fragment")
                        
                        fragment.append(index)
                        has_appended = True
                        break
                # Otherwise, if the atom is not in the fragment then add it to the not explored atom
                if not has_appended: not_explored.append(index)
            # Set the current not explored atoms to the new not explored atoms
            not_explored_atoms = not_explored

            # Append the fragment into the entire fragments array
            fragments.append(fragment)
            
            # Reset the node fragments
            node_fragments = []
            
        return fragments

    def to_graph(self, atoms: Atoms) -> nx.Graph:
        """bonds = list of (idx, idy, type: AtomGraphEdgeType)"""
        analysis = Analysis(atoms)
        graph_list = [AtomGraphNode(atom) for atom in atoms]
        bond_list = []
        unique_atoms = list(set(atoms.get_chemical_symbols()))
        for atomI in unique_atoms:
            for atomJ in unique_atoms:
                bonds = analysis.get_bonds(atomI, atomJ)
                for bond in bonds[0]:
                    bond_list.append(bond)
        bond_list = list(set(bond_list))
        G = nx.Graph()
        
        if len(atoms.get_positions()) == 1:
            G.add_node(graph_list[0])
        
        for bond in bond_list:
            start = bond[0]
            end = bond[1]
            bond_type = -1
            try:
                bond_type = bond[2]
            except IndexError: bond_type = 7
            a = graph_list[start]
            b = graph_list[end]
            if not G.__contains__(a):
                self.counter+=1
                G.add_node(a)
            if not G.__contains__(b):
                self.counter+=1
                G.add_node(b)
            G.add_edge(a, b, weight=1, bond_type=bond_type)
        return G

    # fn: given our current and next atom should we continue?
    def traverse(self, cur: AtomGraphNode, fn: Callable[[AtomGraphNode, AtomGraphNode], bool]):
        # print("Cur " + cur.getAtom().symbol)
        # print("Cur Coord " + str(cur.getAtom().position))
        atoms = self.graph.neighbors(cur)
        cur.setVisited(True)
        for _next in atoms:
            # print("Next " + _next.getAtom().symbol)
            # print("Next Coord" + str(_next.getAtom().position))
            # We are adding the current node multiple times
            should_continue = fn(cur, _next) and (not _next.isVisited())
            if should_continue: 
                self.traverse(_next, fn)

    def get_spin_multiplicity(self):
        # TODO
        # unpaired_electron_count = self.get_unpaired_electron_count()
        # return unpaired_electron_count + 1
        return 2 if self.get_electron_count() % 2 == 1 else 1

    def get_charges(self):
        # TODO
        return 0

    # def get_unpaired_electron_count(self):
    #     self.reset()
    #     total_unpaired_electron = []
    #     graph_list = [node for node in self.graph.nodes]

    #     def traversal_fn(a: AtomGraphNode, b: AtomGraphNode):
    #         unpaired_electron = self.organic_unpaired_electrons[a.getAtom().symbol]
    #         print("ELECTRON " + str(unpaired_electron))
    #         total_unpaired_electron.append(unpaired_electron)
    #         return True

    #     self.traverse(graph_list[0], traversal_fn)

    #     return sum(total_unpaired_electron)

    def get_electron_count(self):
        self.reset()
        total_electron = []
        graph_list = [node for node in self.graph.nodes]

        def traversal_fn(a: AtomGraphNode, b: AtomGraphNode):
            unpaired_electron = self.electrons[a.getAtom().symbol]
            print("ELECTRON " + str(unpaired_electron))
            total_electron.append(unpaired_electron)
            return True

        self.traverse(graph_list[0], traversal_fn)

        return sum(total_electron)
    
    def size(self):
        return self.graph.number_of_nodes()

    def get_graph(self): return self.graph


