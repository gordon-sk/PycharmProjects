"""
Created on Apr 11, 2012

@author: rothgh, heminggs
"""


class BST:
    """ Binary search tree implementation
        each node in the tree has a data object and left & right children (also trees) """

    def __init__(self, data):
        """ Constructor
            create a root node with the given data. """
        self.data = data
        self.left = None
        self.right = None

    def print_tree(self):
        """ Print data in the tree in increasing order using in-order traversal. """
        if self.left:   # same as: if self.left is not None:
            self.left.print_tree()
            print("left")
        print(self.data)
        if self.right:
            self.right.print_tree()
            print("right")

    def format_tree(self, indent=0):
        """ Print data in the tree in a sideways tree format. """
        if self.right:
            self.right.format_tree(indent+3)
        print(' ' * indent, self.data)
        if self.left:   # same as: if self.left is not None:
            self.left.format_tree(indent+3)

    def search(self, target):
        """ Search the tree for target and return True if found, else return False.
            Recursive version of search. """
        if target == self.data:
            return True
        elif target < self.data:
            if self.left:
                return self.left.search(target)
            else:
                return False
        else:  # target > self.data
            if self.right:
                return self.right.search(target)
            else:
                return False

    def insert(self, data):
        """ Insert data into the binary tree by adding a new node. """
        if self.data is None:     # if root was deleted, or someone built an empty root
            self.data = data
        elif data == self.data:
            raise ValueError("Value already exists")
        elif data < self.data:
            if self.left is None:
                self.left = BST(data)
            else:
                self.left.insert(data)
        else:
            if self.right is None:
                self.right = BST(data)
            else:
                self.right.insert(data)

    def __eq__(self, other_tree):
        """ Compare two trees for equality """
        if other_tree is None:
            return False
        if self.data != other_tree.data:
            return False
        result = True
        if self.left:
            result = (self.left == other_tree.left)
        elif other_tree.left:
            return False
        if result and self.right:
            result = (self.right == other_tree.right)
        elif other_tree.right:
            return False
        return result

    def print_larger(self, data):
        if self.data > data:
            if self.right:
                self.right.print_larger(data)
            print(self.data)
            if self.left:
                self.left.print_larger(data)
        else:
            if self.right:
                self.right.print_larger(data)



if __name__ == '__main__':
    #create a simple tree
    tree = BST(10)
    tree.insert(8)
    tree.insert(7)
    tree.insert(3)
    tree.insert(5)
    tree.insert(6)
    tree.insert(1)
    tree.insert(12)
    tree.insert(11)
    tree.insert(15)

    tree.print_larger(5)
    print()

    #print it
    tree.print_tree()
    print()
    print()
    tree.format_tree()

    print(tree.search(3))
    print(tree.search(1))
    print(tree.search(12))
    print(tree.search(10))
    print(tree.search(15))
    print(tree.search(7.5))
    print(tree.search(13))

    tree2 = BST(10)
    tree2.insert(8)
    tree2.insert(7)
    tree2.insert(3)
    tree2.insert(5)
    tree2.insert(6)
    tree2.insert(1)
    tree2.insert(12)
    tree2.insert(11)
    tree2.insert(15)

    print()
    print(tree == tree2)

    tree2.insert(100)

    print(tree==tree2)