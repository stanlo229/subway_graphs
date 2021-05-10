from typing import List

from .molecule_node import MoleculeNode
from .reaction_node import ReactionNode


class Graph:
    def __init__(
        self, molecule_nodes: List[MoleculeNode], reaction_nodes: List[ReactionNode]
    ):
        self.molecule_nodes = molecule_nodes
        self.reaction_nodes = reaction_nodes

    @staticmethod
    def from_json_files(molecule_json_path, reaction_json_path):
        pass
