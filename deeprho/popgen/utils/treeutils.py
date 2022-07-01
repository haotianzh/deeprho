from ..base import BaseTree, Node
import warnings
import time

"""
author: Haotian
created_at: 12/4/2021
description: utils for trees
"""


class TraversalGenerator(object):
    """
    A generator class used for tree traversal
    Arguments:
        order: traversal order
    >>> generator = TraversalGenerator(order='post')
    >>> for node in generator(tree):
            # to do something
            pass
    """

    def __init__(self, order='post'):
        self.order = order
        self.iterator = None

    def __call__(self, tree):
        # calling function
        self.tree = tree
        self.iterator = iter(self._method())
        return self.iterator

    def _pre(self):
        # pre-order traverse
        warnings.warn("no implementation currently")

    def _in(self):
        # mid-order traverse
        warnings.warn("no implementation currently")

    def _post(self):
        node = self.tree.root
        traverse_nodes = []
        while len(traverse_nodes) != len(self.tree):
            while node.get_children():
                flag = True
                for child in node.get_children():
                    if child.identifier not in traverse_nodes:
                        node = child
                        flag = False
                        break
                if flag:
                    break
            traverse_nodes.append(node.identifier)
            yield node
            node = node.parent

    @property
    def order(self):
        return self._order

    @order.setter
    def order(self, method):
        valid_methods = {'pre': self._pre, 'in': self._in, 'post': self._post}
        if method in valid_methods:
            self._order = method
            self._method = valid_methods[method]
        else:
            raise Exception("order should be in ['pre', 'mid', 'post']")


def from_node(node: Node) -> BaseTree:
    """ Build a tree by directly setting a root """
    tree = BaseTree()
    tree.root = node
    tree._nodes = node.get_descendants_dict()
    return tree


def from_newick(newick: str) -> BaseTree:
    """ Build a tree according to a newick-format string """
    def _isvalid(s):
        checking_stack = []
        for ch in s:
            checking_stack.append(ch) if ch == '(' else None
            if ch == ')':
                if checking_stack:
                    checking_stack.pop()
                else:
                    return False
        return True if not checking_stack and ch == ';' else False

    def _next(i):
        stop_words = [',', ')', ';']
        if newick[i] in stop_words:
            return None, 0, i
        br = 0
        j = i + 1
        while newick[j] not in stop_words:
            j += 1
        if ':' in newick[i:j]:
            nid, br = newick[i:j].split(':')
        else:
            nid = newick[i:j]
        return nid.strip(), float(br), j

    newick = newick.strip()
    assert isinstance(newick, str), Exception('newick should be a string.')
    assert _isvalid(newick), Exception('invalid newick string.')
    nodes = []
    level, i, key = 0, 0, ''
    while not newick[i] == ';':
        if newick[i] in [',', ' ']:
            i += 1
            continue
        if newick[i] is '(':
            level += 1
            i += 1
            continue
        if newick[i] is ')':
            if not nodes:
                raise Exception('newick: bad parsing.')
            identifier, branch, end = _next(i + 1)
            node = Node(identifier=identifier, branch=branch)
            while nodes:
                if nodes[-1][1] == level:
                    child = nodes.pop()
                    node.add_child(child[0])
                    child[0].set_parent(node)
                else:
                    break
            level -= 1
            nodes.append((node, level))
            i = end
            continue
        identifier, branch, end = _next(i)
        node = Node(identifier=identifier, branch=branch)
        nodes.append((node, level))
        i = end

    root = nodes.pop()
    tree = from_node(root[0])
    return tree
