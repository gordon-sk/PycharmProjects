"""
Created on Apr 19, 2012

@author: rothgh, heminggs
"""

import sys
from lec17_queues import Queue

##############################################################################
# Vertex class
##############################################################################


class Vertex:
    """Representation of a graph vertex"""
    
    def __init__(self, num):
        """constructor expects an identifier for this vertex"""
        self.id = num
        self.adj = {}             # list of adjacent vertices with edge weights
        self.color = 'white'      # some algorithms want to mark a vertex
        self.dist = sys.maxsize   # some algorithms determine a distance for each vertex
        self.pred = None          # some algorithms determine a predecessor

    def add_neighbor(self, nbr, weight=0):
        """add a neighbor (with edge weight) to this vertex"""
        self.adj[nbr] = weight

    def __str__(self):
        nbrs = []
        for nbr in self.get_adjacent():
            nbrs.append(nbr.get_id())
        if self.pred:   
            pred = self.pred.get_id()
        else:
            pred = None
        return str(self.id) + ":color " + self.color + \
            ":dist " + str(self.dist) + ":pred [" + str(pred) + "]" + ':nbrs ' + str(nbrs)

    # accessors and mutators
    def get_id(self):
        return self.id
    
    def get_adjacent(self):
        return self.adj.keys()
    
    def get_weight(self, nbr):
        return self.adj[nbr]
    
    def get_color(self):
        return self.color
    
    def set_color(self, color):
        self.color = color
        
    def get_distance(self):
        return self.dist
    
    def set_distance(self, dist):
        self.dist = dist
        
    def get_pred(self):
        return self.pred
    
    def set_pred(self, pred):
        self.pred = pred
        
    def reset(self):
        self.color = 'white'     
        self.dist = sys.maxsize
        self.pred = None




##############################################################################
# Graph class
##############################################################################

class Graph:
    """A graph representation using adjacency lists"""
    
    def __init__(self):
        """Constructor for graph class"""
        self.vertices = {}         # use a dictionary to map vertex names to Vertex object
        self.num_vertices = 0
        
    def add_vertex(self, key):
        """Add a vertex with give key to the graph"""
        if key in self:
            raise ValueError("Vertex already exists")
        self.num_vertices += 1
        new_vertex = Vertex(key)
        self.vertices[key] = new_vertex
        return new_vertex
    
    def get_vertex(self, key):
        """Get a vertex with a given key from the graph"""
        if key in self.vertices:
            return self.vertices[key]
        else:
            return None
        
    def __contains__(self, n):
        return n in self.vertices
    
    def has_key(self, key):
        """Determine if the graph has a vertex with a certain key"""
        return key in self.vertices
    
    def add_edge(self, f, t, c=0):
        """Add an edge to the graph from the f vertex to the t vertex with cost c"""
        if f not in self.vertices:
            self.add_vertex(f)
        if t not in self.vertices:
            self.add_vertex(t)
        self.vertices[f].add_neighbor(self.vertices[t], c)
        
    def get_vertices(self):
        """Get a list of all the vertices of the graph"""
        return self.vertices.keys()
    
    # def __iter__(self):
    #     """Define an iterator for the graph"""
    #     return self.vertices.itervalues()

    def BFTraversal(self, start):
        node = self.set_vertex(self, start):
        node.color = 'black'
        print(node)
        q = Queue()
        [q.enqueue(n) for n in node.get_adj()]
        while not q.is_empty():
            node = q.dequeue()
            if node.get_color() == 'white':
                print(node)
                node.set_color('black')
                [q.enqueue() for n in node.get_adjacent()]

    def DFTraversal(self, start):
        s = Stack()
        node = self.set_vertex(start)
        node.set_color('black')
        print(node)
        [s.push(n) for n in node.get_adjacent()]
        while not s.is_empty():
            node = s.pop()
            if node.get_color == 'white':
                print(node)
                node.set_color('black')
                [s.push(n) for n in node.get_adjacent()]


##############################################################################
# Depth-first search
##############################################################################

def dfs(g, vert_key, dist=0):
    pass

    
##############################################################################
# Breadth-first search
##############################################################################

def bfs(g, vert_key):
    pass


##############################################################################
# Djikstra's Path Finding Algorithm
##############################################################################


def dijkstra(g, start_key, second_key):
    return None


if __name__ == '__main__':
    
    gr = Graph()   # create a graph
    
    # add some edges
    gr.add_edge(1, 2, 3)
    gr.add_edge(1, 3, 5)
    gr.add_edge(1, 4, 2)
    gr.add_edge(2, 3, 6)
    gr.add_edge(3, 5, 2)
    gr.add_edge(4, 5, 1)
    gr.add_edge(5, 6, 1)
    gr.add_edge(6, 'a', 10)