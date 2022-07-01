from collections import OrderedDict, defaultdict
from .node import Node
import pptree


class BaseTree(object):
    """
    A tree class, for extracting features from genealogical trees in the future.
    Arguments:
        args:
    # >>> tree = BaseTree();
    # >>> tree = popgen.utils.treeutils.from_newick('((1,2),(3,4));') # or
    """
    def __init__(self):
        self.root = None
        self._nodes = OrderedDict()

    def __contains__(self, node):
        # if a node is in this tree
        if isinstance(node, Node):
            node = node.identifier
        if node in self._nodes:
            return True
        else:
            return False

    def __getitem__(self, identifier) -> Node:
        # get a node from node list
        return self._nodes[identifier]

    def __len__(self):
        # obtain number of nodes in a tree
        return len(self._nodes)

    def create_node(self, identifier=None, name=None, parent=None) -> Node:
        node = Node(identifier=identifier, name=name)
        self.add_node(node, parent)
        return node

    def add_node(self, node, parent=None):
        if node.identifier in self._nodes:
            raise Exception('cannot add the node that has already been in the tree.')
        if parent is None:
            if self.root:
                raise Exception('root has already existed and parent cannot be none.')
            self.root = node
            self._nodes[node.identifier] = node
            # /* set root and make its level as 0 */
            node.set_parent(None)
            return
        pid = parent.identifier if isinstance(parent, Node) else parent
        if not pid in self._nodes:
            raise Exception('parent not found in this tree.')
        self._nodes[node.identifier] = node
        # /* link node with its parent */
        node.set_parent(self[pid])
        self[pid].add_child(node)
        return

    def get_all_nodes(self):
        return self._nodes

    def get_leaves(self):
        # leaves = []
        # for nid in self.get_all_nodes():
        #     if self[nid].is_leaf():
        #         leaves.append(nid)
        # return leaves
        leaves = [node.identifier for node in self.root.get_leaves()]
        return leaves

    def get_splits(self):
        splits = set()
        for nid in self.get_all_nodes():
            if not self._nodes[nid].is_leaf() and not self._nodes[nid].is_root():
                splits.add(frozenset([node.identifier for node in self._nodes[nid].get_leaves()]))
        return splits

    def to_dict(self):
        # return a dict for the whole tree.
        pass

    def output(self, output_format='newick', branch_lengths=False):
        def _newick(node, branch_lengths):
            if node.is_leaf():
                return node.name
            fstr = '(' + ','.join(['{newick}:{branch}'
                                  .format(newick=_newick(child, branch_lengths), branch=child.branch)
                                   if branch_lengths else _newick(child, branch_lengths)
                                   for child in node.get_children()]) + ')'
            return fstr

        def newick():
            return _newick(self.root, branch_lengths) + ';'

        funcs = {'newick': newick}
        return funcs[output_format]()

    def print(self):
        pptree.print_tree(self.root, "_children", horizontal=False)

