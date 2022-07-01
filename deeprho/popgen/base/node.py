from collections import OrderedDict
import os
import binascii
import time

class Node(object):
    """
    A node class, was implemented for constructing genealogical trees.
    Arguments:
        args:
    >>> node = Node()
    >>> node = Node(1) # or
    >>> node = Node(identifier=1, name='ken', branch=10) # or
    """
    def __init__(self, identifier=None, name=None, branch=0):
        self.identifier = identifier
        if name is None:
            self.name = str(self.identifier)
        else:
            self.name = name
        self.parent = {}
        self.branch = branch
        self.children = OrderedDict()
        self._children = []  # using for printing tree

    @property
    def identifier(self):
        return self._identifier

    @identifier.setter
    def identifier(self, nid):
        if nid is None or nid is '':
            nid = str(binascii.hexlify(os.urandom(5)).decode('utf-8'))
        self._identifier = nid

    @property
    def branch(self):
        return self._branch

    @branch.setter
    def branch(self, value):
        if self.parent is None:
            raise Exception('node has no parent, branch cannot be assigned.')
        if isinstance(value, float) or isinstance(value, int):
            self._branch = value
        else:
            raise Exception('branch must be a number.')

    def is_root(self):
        return self.get_level() == 0

    def is_leaf(self):
        return len(self.children) == 0

    def get_ancestors(self):
        parents = []
        node = self.parent
        while node:
            parents.append(node)
            node = node.parent
        return parents

    def get_descendants_dict(self):
        # simple bfs for searching descendants
        descendants = {}
        queue = [self]
        while queue:
            node = queue.pop(0)
            for child in node.get_children():
                descendants.update({child.identifier: child})
                queue.append(child)
        return descendants

    def get_descendants(self):
        # simple bfs for searching descendants
        descendants = []
        queue = [self]
        while queue:
            node = queue.pop(0)
            for child in node.get_children():
                descendants.append(child)
                queue.append(child)
        return descendants

    def add_child(self, node):
        """
        Adding child to a specific node. Potentially running slow due to inner traversal for ensuring no conflicts
        happened. (should be optimized in the future)
        """
        if node.identifier in [d.identifier for d in self.get_descendants()]:
            raise Exception('node %s has already been added.' % node.identifier)
        if node.identifier in [a.identifier for a in self.get_ancestors()]:
            raise Exception('parent %s cannot be added as a child.' % node.identifier)
        self.children[node.identifier] = node
        self._children.append(node)

    def get_children(self):
        return self._children

    def get_leaves(self):
        leaves = []
        for node in self.get_descendants():
            if node.is_leaf():
                leaves.append(node)
        return leaves

    def get_level(self):
        level = 0
        node = self.parent
        while node:
            node = node.parent
            level += 1
        return level

    def set_parent(self, parent):
        if parent is None:
            self.parent = None
        else:
            self.parent = parent

    def set_branch(self, value):
        self.branch = value


if __name__ == '__main__':
    st = time.time()
    pre = Node()
    for i in range(10000):
        node = Node()
        pre.add_child(node)
        node.set_parent(pre)
        pre = node

    print(time.time()- st)
