## Implementing the datastructures for nodes, edges... "topology" descrbibed in 2016 multilevel hp-refinement paper cited by Notaros and Harmom
# 4 classes for edge node face elements and topology parent class which defines most aspects of the RBS procedure. 
# deactivate all topological components of the overlay mesh whose adjacecy list contains elements of different levels.
import numpy as np


class topology:
    #store base mesh, and mesh levels in general
    def __init__(self, id,active, level=0, parent=None):
        self.id = id
        self.children = []
        self.level = level
        self.parent = parent
        self.active = True
    @property
    def is_leaf(self):
        #if len(self.children) == 0:
        return len(self.children) == 0



class node(topology):
    # In 1D a node is a face -> important for H(div)
    def __init__(self, id, coordinate, level=0, parent=None):
        super().__init__(id, level, parent)
        self.coordinate = coordinate
  

        
        

class edge(topology):
    def __init__(self, id, nodes, level, parent):
        super().__init__(id, level, parent)
        self.nodes = nodes

#class face(topology):


class element(topology):
    def __init__(self, id,nodes,p_order, level, parent):
        super().__init__(id, level, parent)
        self.p_order = p_order
        self.nodes = nodes

    def deactivate(self):
        if self.active == False:
            p_order = 0 
        return p_order

    def boundary_cond(self):

        return []

class mesh:
    def __init__(self,L, N, p, layer=0):
        self.L = L
        self.N = N
        self.p = p
        self.layer = layer
        self.node_set = []
        self.elems = []

        

    def gen_init_mesh(self):
        nodes = np.linspace(0,self.L,self.N)
        n_index = np.linspace(0,len(nodes))
        #lems = []
        #ode_set = []
        for i in range(len(nodes)):
            self.node_set.append(node(i, nodes[i],level=self.layer))
        for j in range(1,len(nodes)):
            self.elems.append(element(j-1, [self.node_set[j-1], self.node_set[j]], p_order = self.p, level=self.layer, parent=None))
        return self.node_set, self.elems

    

class refiner:
    def __init__(self, mesh, marked_elem_h, marked_elem_p):
        self.mesh = mesh
        self.marked_elem_h = marked_elem_h
        self.marked_elem_p = marked_elem_p

    def refine_h(self):
        elements, nodes = self.mesh
        level_new = self.mesh.layer +1
        for i in self.marked_elem_h:
            parent = self.mesh.elems[i]
            newpos = parent.nodes.coordinate[1]  + parent.nodes.coordinate[0] / 2
            new_id = self.mesh.nodeset[-1].id
            child1_id = self.mesh.elems[-1].id + 1
            child2_id = self.mesh.elems[-1].id + 2
            new_node= node(new_id, newpos,level = level_new)
            self.mesh.nodeset.append(new_node)
            child1 = element(child1_id, [parent.nodes.coordinate[0] , new_node], p_order = parent.p_order, parent=parent, level=level_new, active=True)
            child2 = element(child2_id,[new_node,parent.nodes.coordinate[1]],  p_order = parent.p_order,level=level_new,parent=parent, active=True)
            self.mesh.elems.append(child1)
            self.mesh.elems.append(child2)
            parent.children.extend([self.mesh.elems[-2],self.mesh.elems[-1]])
            parent.active = False
            return self.mesh

    def refine_p(self):
        
        for i in self.marked_elem_p:
            self.mesh.elems[i].p_order +=1
            return self.mesh




            
mesh1 = mesh(10,20, 2)
nodes, elems = mesh1.gen_init_mesh()




# next need way to assing dofs on each element according the p order of the element,  refining p is easy now - simply p+=1!


def h_refine(marked, mesh):
    elements = mesh.elems










# manually creating a small mesh to larn how the classes work together for 1D
n0 = node(0, 0)
n1 = node(1,1)

n2 = node(2,2)
n3 = node(3,3)
n4= node(4,4)

e1 = [n0,n1]
e2 = [n1,n2]
e3 = [n2, n3]
e4= [n3, n4]

# creating each element
E1 = element(0, e1, p_order=2,level=0,parent=None)
E2 = element(1, e2, p_order=2, level=0,parent=None)
E3 = element(2, e3, p_order = 2, level=0, parent=None)
E4 = element(3, e4, p_order=2, level=0, parent=None)

# manual h-refinement:
n5 = node(5,1.5,level=1)
e5 = [n1,n5]
e6 = [n5,n2]
E5  = element(4,e5, p_order=2, level=1,parent =E2)
E6 = element(5, e6,p_order=2, level=1,parent=E2)

E2.children = [E5,E6]
E2.active= False

