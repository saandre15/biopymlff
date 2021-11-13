from ase.atoms import Atoms
from ase.atom import Atom

import networkx as nx

from typing import Callable

from biopymlff.util.AtomGraphNode import AtomGraphNode


class AtomGraph():

    def __init__(self, atoms: Atoms):
        self.graph = self.to_graph(atoms)
        self.graph.neighbors(n)

    # Returns a list of atoms
    def fragments_by_bond(self, symbolA: str, symbolB: str) -> list:
        all_atoms = to_list(mol)
        explored_atoms = []
        not_explored_atoms = range(1, len(all_atoms))
        fragments = []
        cur_fragment = []
        system_size = len(not_explored_atoms)

        def traversal_fn(cur: AtomGraphNode, next: AtomGraphNode):
            curSymbol = cur.getAtom().get_chemical_symbol()
            nextSymbol = next.getAtom().get_chemical_symbol()
            cur_fragment.append(curSymbol)
            return not ((curSymbol == symbolA and nextSymbol == symbolB) or (curSymbol == symbolB and nextSymbol == symbolA))

        while len(not_explored_atoms) > 0:
            index = 0
            self.traverse(AtomGraphNode(not_explored_atoms[index]), traversal_fn)
            fragments.append(cur_fragment)
            cur_fragment = []

    def to_graph(atoms: Atoms) -> nx.Graph:
        analysis = Analysis(mol)
        bond_list = []
        unique_atoms = list(set(atoms.get_chemical_symbols()))
        for atomI in unique_atoms:
            for atomJ in unique_atoms:
                bonds = analysis.get_bonds(atomI, atomJ)
                for bond in bonds[0]:
                    bond_list.append(bond)
        bond_list = list(set(bond_list))
        G = nx.Graph()
        for bond in bond_list:
            start = bond.split("-")[0]
            end = bond.split("-")[1]
            a = AtomGraphNode(atoms[start])
            b = AtomGraphNode(atoms[end])
            G.add_node(a)
            G.add_node(b)
            G.add_edge(a, b, weight=1)

    # fn: given our current and next atom should we continue?
    def traverse(self, cur: AtomGraphNode, fn: Callable[[AtomGraphNode, AtomGraphNode], bool]):
        cur.visited = True
        atoms = self.graph.neighbors(cur)
        for _next in atoms:
            should_continue = fn(cur, _next) or (not _next.visited())
            if should_continue:
                self.traverse(_next, fn)

    def get_spin_multiplcity(self):
        return 1 if self.get_electron_count() % 2 == 1 else 2

    def get_electron_count():
        pass # 

