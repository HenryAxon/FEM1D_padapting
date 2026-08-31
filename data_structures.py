## Implementing the datastructures for nodes, edges, faces, solids... "topology" as descrbibed in 2016 multilevel hp-refinement paper from Nils Zander et al. 
# deactivate all topological components of the overlay mesh whose adjacecy list contains elements of different levels.
# Currently can create mesh and refine in h and p, NOT ADAPTIVE. 
# 
# NEED TO IMPLEMENT: 
# 1) degrees of freedom deactivation on elements
# 2) basis function implementation on elements
# 3) the function space continuity requirements of H(div),normal/Raviart THomas space and H(curl), tangential/Nedelec space, which should be relatively simple, but the goal 
# is to stitch the two together without an if else statement,
# presumably somehow based on information that is included from the initial setup of the problem? 
# 4) reread the Notaros implementation and compare to Ross code structure to learn the anisotropic version and identify where the 2D version will be implemented.



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
        self.adjacent_elem = []
  

    
class edge(topology):
    def __init__(self, id, nodes, level, parent):
        super().__init__(id, level, parent)
        self.nodes = nodes

#class face(topology):


class element(topology):
    def __init__(self, id,nodes,p_order, level, parent):
        super().__init__(id=id, active=True, level=level, parent=parent)
        self.p_order = p_order
        self.nodes = nodes
        self.adjacency_list = []

    def deactivate(self):
        if self.active == False:
            p_order = 0 
        return p_order





class mesh:
    def __init__(self,L, N, p, layer=0):
        self.L = L
        self.N = N
        self.p = p
        self.layer = layer
        self.node_set = []
        self.elems = []
        #self.active_elems = []
        

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
    @property
    def active_elements(self):
        return [e for e in self.elems if e.active]

    @property
    def active_nodes(self):
        return [n for n in self.node_set if n.active]

    @property
    def jacobian_finder(self):
        jacobians = []
        for i in range(len(self.active_elements)):
            nodes = self.active_elements[i].nodes
            coord_1 = nodes[0].coordinate
            coord_2 = nodes[1].coordinate
            jacobians.append((coord_2 - coord_1) / 2)
        return jacobians

    @property
    def gaussian_integral_points(self):
        '''
        property to find the gaussian integral points on each element to then perform integration to form the 
        K, M matrices for solving problems on the active elements
        '''
        gauss_info = []
        
        for i in range(len(self.active_elements)):
            [x,w] = np.polynomial.legendre.leggauss(20)
            gauss_info.append([x,w])

        return gauss_info




    def adjacency_info(self):
        for element_item in self.elems:
            element_item.adjacency_list = []

        for index, element_item in enumerate(self.elems):
            if index > 0:
                element_item.adjacency_list.append(self.elems[index - 1])
            if index < len(self.elems) - 1:
                element_item.adjacency_list.append(self.elems[index + 1])

        for node_item in self.node_set:
            node_item.adjacent_elem = [
                element_item.id
                for element_item in self.elems
                if node_item in element_item.nodes
        ]

    def deactivate(self):
        '''
        Meant to deactivate the child elements overlapping nodes to allow the parent elements nodes to 
        have precedence in degree of freedom.
        '''
        for node_item in self.node_set:
            adjacent = node_item.adjacent_elem
            if adjacent.level < node_item.level and node_item in adjacent.node_set:
                node_item.active = False
                node_item.p_order = 0
    
    def plot_mesh(self, ax=None, annotate=True, show_nodes=False):
        if ax is None:
            fig, ax = plt.subplots()
        for element in self.elems:
            element.plot(ax=ax, annotate=annotate, show_nodes=show_nodes)
        return ax

            

# all we are really doing in 1D is technically just 1 edge refinement so to speak, and there is no directionality - which of course limits the use of this code, given that the 
# actual goals are to produce effective anisotropic hp refinement. 

    

class refiner:
    def __init__(self, mesh, marked_elem_h, marked_elem_p):
        self.mesh = mesh
        self.marked_elem_h = marked_elem_h
        self.marked_elem_p = marked_elem_p

    def refine_h(self):
        #elements, nodes = self.mesh
        level_new = self.mesh.layer +1
        for i in self.marked_elem_h:
            parent = self.mesh.elems[i]
            left = parent.nodes[0]
            right = parent.nodes[1]

            print("LEFT:", left.id, left.coordinate)
            print("RIGHT:", right.id, right.coordinate)

            newpos = (left.coordinate + right.coordinate) / 2

            print("NEW:", newpos)
            new_id = len(self.mesh.node_set)
            child1_id = len(self.mesh.elems)
            child2_id = len(self.mesh.elems) + 1
            new_node= node(new_id, newpos, level = level_new)
            self.mesh.node_set.append(new_node)
            child1 = element(child1_id, [parent.nodes[0] , new_node], p_order = parent.p_order, parent=parent, level=level_new)
            child2 = element(child2_id, [new_node,parent.nodes[1]],  p_order = parent.p_order,level=level_new,parent=parent)
            self.mesh.elems.append(child1)
            self.mesh.elems.append(child2)
            parent.children.extend([self.mesh.elems[-2],self.mesh.elems[-1]])
            parent.active = False
            parent.p_order = 0
        self.mesh.adjacency_info()
        return self.mesh

    def refine_p(self):
        
        for i in self.marked_elem_p:
            self.mesh.elems[i].p_order +=1
        return self.mesh



def basis_eval(j, p, s):
    """Hierarchical shape functions on the reference element [-1, 1]."""

    if j == p+1:
        f_base = 0.5 * (1 + s)
        df_base = 0.5
    elif j == 1:
        f_base = 0.5 * (1 - s) * (1 + s) ** (j-1)
        df_base = -0.5
    else:
        f_base = 0.5 * (1 - s) * (1 + s) ** (j-1)
        df_base = 0.5 * (-(1 + s) ** (j-1) + (j-1) * (1 - s) * (1 + s) ** (j - 2))
        
        #print(f"Basis function {j}: f_base = {f_base}, df_base = {df_base}")
    return f_base, df_base





def K_integrals(p, i, j, alpha, beta, jacob):
    #'Integrals via gaussian quadrature for stiffness matrix', using 20 points purely to match Cam Key results
    x, w = np.polynomial.legendre.leggauss(20)
    integral = 0.0
    for q in range(len(x)):
        #f, df = basis_eval( n, x[q])
        # mass-like term uses the jacobian, stiffness-like term divides by it
        # (physical derivative = reference derivative / jacobian)
        integral += w[q] * ( basis_eval(i, p, x[q])[1] * basis_eval(j, p, x[q])[1])
    integral = -integral*alpha/jacob
    return integral
 
 
def M_integrals(p, i, j, beta, jacob):
    #'perform integrals via gaussian quadrature for mass matrix'
    x, w = np.polynomial.legendre.leggauss(20)
    integral = 0.0
    for q in range(len(x)):
        #f, df = basis_eval(i, n, x[q])
        #print(f)
        integral += w[q] *  basis_eval(i, p, x[q])[0] * basis_eval(j, p, x[q])[0]
    integral = integral*beta*jacob
    return integral
 
 
def K_fill_local(p_order, alpha, beta, jacob):
    #'Fill each element stiffness matrix'
    K_loc = np.zeros((p_order+1, p_order+1))
    for i in range(1, p_order+2):
        for j in range(1, p_order+2):
            K_loc[i-1, j-1] = K_integrals(p_order, i, j, alpha, beta, jacob)
    return K_loc
 
 
def M_fill_local(p_order, beta, jacob):
    #'Filling each block for each element of M_ij'
    M_loc = np.zeros((p_order+1, p_order+1))
    for i in range(1, p_order+2):
        for j in range(1, p_order+2):
            M_loc[i-1, j-1] = M_integrals(p_order, i, j, beta, jacob)
    return M_loc
 
 
def K_fill_global(M, meshy, alpha=1, beta=1):
    #'Stiffness matrix builder '
    nodes, elements, jac, p_orders, dof = meshy
    K_globe = np.zeros((len(nodes), len(nodes)))
    for e in range(M):
        elem = elements[e]
        current_jacob = jac[e]
        current_p_order = p_orders[e]
        K_loc = K_fill_local(current_p_order, alpha, beta, current_jacob)
        for i in range(current_p_order+1):
            for j in range(current_p_order+1):
                I = elem[i]
                J = elem[j]
                K_globe[I, J] += K_loc[i, j]
    return K_globe
 
 
def M_fill_global(M, meshy, beta=1):
    #"Mass matrix builder"
    nodes, elements, jac, p_orders,dof = meshy
    M_globe = np.zeros((len(nodes), len(nodes)))
    for e in range(M):
        elem = elements[e]
        current_jacob = jac[e]
        current_p_order = p_orders[e]
        M_loc = M_fill_local(current_p_order, beta, current_jacob)
        for i in range(current_p_order+1):
            for j in range(current_p_order+1):
                I = elem[i]
                J = elem[j]
                M_globe[I, J] += M_loc[i, j]
    return M_globe
 
 
def weights_unknown(N):
    return np.zeros(N)
 
 
def G_assemble(g, meshy):
    #"excitation vector ==0 unless other g specified"
    nodes, elements, jac, p_orders, dof = meshy
    G = np.zeros(dof)
    if g != 0:
        G[:] = g
    return G
 
 







mesh1 = mesh(4,5, 2)
mesh1.gen_init_mesh()
marked_elem_p = [0,3]
marked_elem_h = [1]

refine = refiner(mesh1, marked_elem_h, marked_elem_p)

refine.refine_h()
refine.refine_p()

activee = mesh1.active_elements
activen = mesh1.active_nodes




for e in refine.mesh.active_elements:
    print(
        "Element:", e.id,
        "nodes:", [n.id for n in e.nodes],
        "coordinates:", [n.coordinate for n in e.nodes],
        "p:", e.p_order,
        "level:", e.level,
        "active:", e.active,
        "adjacency list", e.adjacency_list
    )


for e in activen:
    print(
        "Element", e.id,
        "nodes:", [n.id for n in e.nodes],
        "p:", e.p_order,
        "level:", e.level,
        "active:", e.active
    )

