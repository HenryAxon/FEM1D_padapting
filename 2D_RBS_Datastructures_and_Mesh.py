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
import matplotlib.pyplot as plt


class topology:
    '''
    Most generally the storage of topology information, primarily meant to be inherited by each type of topology
    '''
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
    '''
    Nodes are the 1D interface where continuity is enforced it mostly just stores this info to be used as the endpoints for edges
    '''
    # In 1D a node is a face -> important for H(div)
    def __init__(self, id, coordinatex,coordinatey, level=0, parent=None):
        super().__init__(id, active=True, level=level, parent=parent)
        self.coordinatex = coordinatex
        self.coordinatey = coordinatey
        self.adjacent_elem = []
  

    
class edge(topology):
    '''
    The edges are the interfaces between elements in 2D.
    '''
    def __init__(self, id, nodes, level=0, parent=None):
        super().__init__(id, active=True, level=level, parent=parent)
        self.nodes = nodes
        self.adjacency = []

class face(topology):
    '''
    In 3D, we need faces as these are what continuity is enforced over.
    '''
    def __init__(self, id, nodes,edges, level, parent):
        super().__init__(id, active=True, level=level, parent=parent)
        self.edges = edges
        self.nodes = nodes
        self.adjacent = []



#class face(topology):


class element(topology):
    '''
    Defines the elements in a general sense, computing edges, nodes, faces (in 3D).
    it inherits from topology and assembles an element and needs to handle full element deacivation.
    
    '''
    def __init__(self, id,nodes,edges,faces,p_order_u, p_order_v, level, parent):
        super().__init__(id=id, active=True, level=level, parent=parent)
        self.p_order_u = p_order_u
        self.p_order_v = p_order_v
        self.nodes = nodes
        self.edges = edges
        self.faces = faces
        self.adjacency_list = []

    def deactivate(self):
        if self.active == False:
            p_order_u = 0
            p_order_v = 0
        return p_order_u, p_order_v

    def plot(self, ax, annotate=True, show_nodes=False):
        # nodes are stored as [ (i,j), (i,j+1), (i+1,j), (i+1,j+1) ]
        # reorder into a closed loop: bottom-left, bottom-right, top-right, top-left
        bl, br, tl, tr = self.nodes[0], self.nodes[1], self.nodes[2], self.nodes[3]
        loop = [bl, br, tr, tl, bl]
        xs = [n.coordinatex for n in loop]
        ys = [n.coordinatey for n in loop]
 
        ax.plot(xs, ys, 'k-', linewidth=0.8)
 
        if show_nodes:
            node_xs = [n.coordinatex for n in self.nodes]
            node_ys = [n.coordinatey for n in self.nodes]
            ax.plot(node_xs, node_ys, 'o', color='tab:blue', markersize=3)
 
        if annotate:
            cx = sum(n.coordinatex for n in self.nodes) / len(self.nodes)
            cy = sum(n.coordinatey for n in self.nodes) / len(self.nodes)
            ax.annotate(str(self.id), (cx, cy), ha='center', va='center', fontsize=8)    





class mesh:
    '''Uses elememts to stitch together a mesh'''
    def __init__(self,L, N, p_u, p_v, level=0):
        self.L = L
        self.N = N
        self.p_u = p_u
        self.p_v = p_v
        self.level = level
        self.node_set = []
        self.edges_u = []
        self.edges_v = []
        self.elems = []
        #self.active_elems = []
        
        self._next_node_id = 0
        self._next_edge_id = 0
        self._next_elem_id = 0    

        self.nodes_by_id = {}
        self.edges_by_id = {}
        self.elems_by_id = {}        

            
    def next_node_id(self):
        i = self._next_node_id
        self._next_node_id += 1
        return i

    def next_edge_id(self):
        i = self._next_edge_id
        self._next_edge_id += 1
        return i

    def next_elem_id(self):
        i = self._next_elem_id
        self._next_elem_id += 1
        return i

    
    def gen_init_mesh(self):
        nodesx = np.linspace(0,self.L,self.N)
        nodesy = np.linspace(0, self.L, self.N)
        [nodesX, nodesY] = np.meshgrid(nodesx, nodesy)
        node_index = lambda i, j: i * self.N + j
        edge_u_index = lambda i, j: i * (self.N - 1) + j
        edge_v_index = lambda i, j: i * self.N + j


        node_at = {}
        edge_u_at = {}
        edge_v_at = {}
        #n_index = np.linspace(0,len(nodes))
        #lems = []
        #ode_set = []
        for i in range(self.N):
            for j in range(self.N):
                self.node_set.append(node(self.next_node_id(), nodesX[i,j],nodesY[i,j],level=self.level))
                n = self.node_set[-1]
                self.nodes_by_id[n.id] = n
                node_at[i, j] = n
        for i in range(self.N):
            for j in range(self.N-1):
                e=edge(self.next_edge_id(), [self.node_set[node_index(i,j)], self.node_set[node_index(i,j+1)]], level=self.level, parent=None)
                self.edges_u.append(e)
                self.edges_by_id[e.id] = e
                edge_u_at[i, j] =e
        for i in range(self.N-1):
            for j in range(self.N):
                e= edge(self.next_edge_id(), [self.node_set[node_index(i,j)], self.node_set[node_index(i+1,j)]], level=self.level, parent=None)
                self.edges_v.append(e)
                self.edges_by_id[e.id] = e
                edge_v_at[i, j] = e

        for m in range(self.N-1):
            for n in range(self.N-1):
                e=element(self.next_elem_id(), nodes=[self.node_set[node_index(m,n)],self.node_set[node_index(m,n+1)], self.node_set[node_index(m+1,n)],self.node_set[node_index(m+1,n+1)]], edges=[self.edges_u[edge_u_index(m,n)],self.edges_u[edge_u_index(m+1,n)],self.edges_v[edge_v_index(m,n)],self.edges_v[edge_v_index(m,n+1)]],faces = [], p_order_u=self.p_u, p_order_v=self.p_v, level=self.level,parent=None)
                self.elems.append(e)
                self.elems_by_id[e.id] = e
        return self.node_set, self.edges_u, self.edges_v, self.elems
    @property
    def active_elements(self):
        return [e for e in self.elems if e.active]

    @property
    def active_nodes(self):
        return [n for n in self.node_set if n.active]

    @property
    def active_edges_u(self):
        return [e for e in self.edges_u if e.active]

    @property
    def active_edges_v(self):
        return [e for e in self.edges_v if e.active]

    # @property
    # def jacobian_finder(self):
    #     jacobians = []
    #     for i in range(len(self.active_elements)):
    #         nodes = self.active_elements[i].nodes
    #         coord_1 = nodes[0].coordinate
    #         coord_2 = nodes[1].coordinate
    #         jacobians.append((coord_2 - coord_1) / 2)
    #     return jacobians

    # @property
    # def gaussian_integral_points(self):
    #     '''
    #     property to find the gaussian integral points on each element to then perform integration to form the 
    #     K, M matrices for solving problems on the active elements
    #     '''
    #     gauss_info = []
        
    #     for i in range(len(self.active_elements)):
    #         [x,w] = np.polynomial.legendre.leggauss(20)
    #         gauss_info.append([x,w])

    #     return gauss_info




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
        for edge in self.edges_u:
            edge.adjacent_elem = [
                element_item.id
                for element_item in self.elems
                if edge in element_item.edges
            ]
        for edge in self.edges_v:
            edge.adjacent_elem = [
                element_item.id
                for element_item in self.elems
                if edge in element_item.edges_v
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
 #               node_item.p_order = 0
        for edge in self.edges:
            adjacent = edge.adjacent_elem
            if adjacent.level < edge.level and edge in adjacent.edges:
                edge.active = False
#                edge.p_order = 0
    
    def plot_mesh(self, ax=None, annotate=True, show_nodes=False):
        if ax is None:
            fig, ax = plt.subplots()
        for elem in self.active_elements:
            elem.plot(ax=ax, annotate=annotate, show_nodes=show_nodes)
        ax.set_aspect('equal')
        return ax

            

# all we are really doing in 1D is technically just 1 edge refinement so to speak, and there is no directionality - which of course limits the use of this code, given that the 
# actual goals are to produce effective anisotropic hp refinement. 
 

    





class refiner:
    '''
    So this takes the mesh and then can refine in h and p and ideally save the correct parents and children to each entity and element.
    '''
    def __init__(self, mesh, marked_elem_h_t,marked_elem_h_u,marked_elem_h_v, marked_elem_p):
        self.mesh = mesh
        self.marked_elem_h_t = marked_elem_h_t
        # self.marked_elem_h_u = marked_elem_h_u
        # self.marked_elem_h_v = marked_elem_h_v
        self.marked_elem_p = marked_elem_p

    def child_nodes(self,parent):
        left_bot = parent.nodes[0]
        left_top = parent.nodes[1]
        right_top = parent.nodes[3]
        right_bot = parent.nodes[2]

        newpos_l = (left_bot.coordinatex, left_bot.coordinatey + left_top.coordinatey /2)
        newpos_r = (right_bot.coordinatex, left_bot.coordinatey + left_top.coordinatey /2)
        newpos_t = (left_bot.coordinatex + left_top.coordinatex /2, left_top.coordinatey)
        newpos_b = (left_bot.coordinatex + left_top.coordinatex /2, left_bot.coordinatey)

        new_nodem = node(
            self.mesh.next_node_id(),
            (left_bot.coordinatex + right_top.coordinatex) / 2,
            (left_bot.coordinatey + right_top.coordinatey) / 2,
            level=self.mesh.level + 1,
        )
        self.mesh.node_set.append(new_nodem)
        self.mesh.nodes_by_id[new_nodem.id] = new_nodem


        return new_nodem

    def split_edge(self, edge_in):
        # Split the edge into two child edges
        left_node = edge_in.nodes[0]
        right_node = edge_in.nodes[1]

        new_pos = ((left_node.coordinatex + right_node.coordinatex) / 2, (left_node.coordinatey + right_node.coordinatey) / 2)
        new_node = node(self.mesh.next_node_id(), new_pos[0], new_pos[1], level=edge_in.level + 1)
        self.mesh.node_set.append(new_node)
        self.mesh.nodes_by_id[new_node.id] = new_node
            

        child_edge1 = edge(self.mesh.next_edge_id(), [left_node, new_node], level=edge_in.level + 1)
        child_edge2 = edge(self.mesh.next_edge_id(), [new_node, right_node], level=edge_in.level + 1)
        if edge_in in self.mesh.edges_u:
            self.mesh.edges_u.append(child_edge1)
            self.mesh.edges_u.append(child_edge2)
            self.mesh.edges_by_id[child_edge1.id] = child_edge1
            self.mesh.edges_by_id[child_edge2.id] = child_edge2
        elif edge_in in self.mesh.edges_v:
            self.mesh.edges_v.append(child_edge1)
            self.mesh.edges_v.append(child_edge2)
            self.mesh.edges_by_id[child_edge1.id] = child_edge1
            self.mesh.edges_by_id[child_edge2.id] = child_edge2

        return new_node, child_edge1, child_edge2






    def refine_h(self):
        #elements, nodes = self.mesh
        level_new = self.mesh.level +1
        for i in self.marked_elem_h_t:
            parent = self.mesh.elems[i]



            center_node = self.child_nodes(parent)

            bottom_midpoint, bottom_left_edge, bottom_right_edge = self.split_edge(parent.edges[0])
            top_midpoint, top_left_edge, top_right_edge = self.split_edge(parent.edges[1])
            left_midpoint, left_bottom_edge, left_top_edge = self.split_edge(parent.edges[2])
            right_midpoint, right_bottom_edge, right_top_edge = self.split_edge(parent.edges[3])


        # manually create the new edges that are not children of any edge but are on the newly created midpoints from splitting edges
            child_edge1 = edge(self.mesh.next_edge_id(), [left_midpoint, center_node], level=level_new)
            child_edge2 = edge(self.mesh.next_edge_id(), [center_node, right_midpoint], level=level_new)
            child_edge3 = edge(self.mesh.next_edge_id(), [bottom_midpoint, center_node], level=level_new)
            child_edge4 = edge(self.mesh.next_edge_id(), [center_node, top_midpoint], level=level_new)
            
            self.mesh.edges_u.append(child_edge1)
            self.mesh.edges_u.append(child_edge2)
            self.mesh.edges_v.append(child_edge3)
            self.mesh.edges_v.append(child_edge4)
            self.mesh.edges_by_id[child_edge1.id] = child_edge1
            self.mesh.edges_by_id[child_edge2.id] = child_edge2
            self.mesh.edges_by_id[child_edge3.id] = child_edge3
            self.mesh.edges_by_id[child_edge4.id] = child_edge4

            # create the 4 new child elements 

            child_elem1 = element(self.mesh.next_elem_id(), nodes=[parent.nodes[0], bottom_midpoint, left_midpoint, center_node], edges=[bottom_left_edge, child_edge1, left_bottom_edge, child_edge3], faces=[], p_order=parent.p_order, level=level_new, parent=parent)
            child_elem2 = element(self.mesh.next_elem_id(), nodes=[bottom_midpoint, parent.nodes[1], center_node, right_midpoint], edges=[bottom_right_edge, child_edge2, child_edge3, right_bottom_edge], faces=[], p_order=parent.p_order, level=level_new, parent=parent)
            child_elem3 = element(self.mesh.next_elem_id(), nodes=[left_midpoint, center_node, parent.nodes[2], top_midpoint], edges=[child_edge1, top_left_edge, left_top_edge, child_edge4], faces=[], p_order=parent.p_order, level=level_new, parent=parent)
            child_elem4 = element(self.mesh.next_elem_id(), nodes=[center_node, right_midpoint, top_midpoint, parent.nodes[3]], edges=[child_edge2, top_right_edge, child_edge4, right_top_edge], faces=[], p_order=parent.p_order, level=level_new, parent=parent)

            # save the new child elemetns to the mesh
            self.mesh.elems.append(child_elem1)
            self.mesh.elems.append(child_elem2)
            self.mesh.elems.append(child_elem3)
            self.mesh.elems.append(child_elem4)

            # add to the dictionary of elements by id
            self.mesh.elems_by_id[child_elem1.id] = child_elem1
            self.mesh.elems_by_id[child_elem2.id] = child_elem2
            self.mesh.elems_by_id[child_elem3.id] = child_elem3
            self.mesh.elems_by_id[child_elem4.id] = child_elem4

            parent.children.extend([self.mesh.elems[-4],self.mesh.elems[-3],self.mesh.elems[-2],self.mesh.elems[-1]])
            parent.active = False
            parent.p_order_u = 0
            parent.p_order_v = 0

    def refine_p(self):
        
        for i in self.marked_elem_p:
            self.mesh.elems[i].p_order_u +=1
            self.mesh.elems[i].p_order_v +=1
        return self.mesh


mesh1 = mesh(L=10, N=6, p_u=2, p_v=2)
mesh1.gen_init_mesh()
refined_1 = refiner(mesh1, marked_elem_h_t=[0,1], marked_elem_h_u=[], marked_elem_h_v=[], marked_elem_p=[2])
refined_1.refine_h()
refined_1.refine_p()


ax1 = refined_1.mesh.plot_mesh(annotate=True, show_nodes=True)
plt.savefig("mesh_preview.png", dpi=150)
print("saved mesh_preview.png")
# in 1d these are essentially the same, but writing the skeleton can still be useful for the 2D implementation




## now need to implement the degrees of freedom adn basis functions for Raviart Thomas spaces to be able to solve a H(div) type problem eventually. 
# implement raviart thomas basis elements by rotation of the 2D Nedelec basis functions of a given order.
class error_indication:
    '''
    Will eventually peform the adjoint refinement marking and all that.
    '''
    def __init__(self, mesh, error_threshold):
        self.mesh = mesh
        self.error_threshold = error_threshold
        self.marked_p = []
        self.marked_h_u = []
        self.marked_h_v = []
        self.marked_h_t = []




class vector_basis_function:
    '''
    Computes basis (need to add derivatives) 
    '''
    def __init__(self, mesh, order_u, order_v):
        self.mesh = mesh
        self.order_u = order_u
        self.order_v = order_v
        self.basis_functions = {}

    def basis_functions(self,u, v, order_u, order_v):
        # this actually should be basis function agnostic. From notaros review paper H(div) basis constrution on quadrilatiral [-1,1] for both dimensions
        # this only created the 1D component in either u or v direction, evaluate the same basis functions for the v direction then you must take the product of 
        # each direction with a constant parameter in the otehr direction to get f_ij for u and v and then multiply them in the summation with the unkowns to get the 
        # final full basis in terms of the solution.  
        # 
        if order_v == 0: 
            v_basis = (1 - v) * u**order_u
        elif order_v == 1:
            v_basis = (v + 1) * u **order_u
        elif order_v >= 2 and order_v % 2 == 0:
            v_basis = (v**order_v - 1) * u **order_u
        else:
            v_basis = (v**order_v - v) * u **order_u

        # now the cases where we have the u direction basis functions

        if order_u == 0:
            u_basis = (1 - u) * v**order_v
        elif order_u == 1:
            u_basis = (u + 1) * v **order_v
        elif order_u >= 2 and order_u % 2 == 0:
            u_basis = (u**order_u - 1) * v **order_v 
        else:
            u_basis = (u**order_u - u) * v **order_v

        return u_basis, v_basis


    def interior_dofs(self, element):
        for i in range(element.p_order_u + 1):
            for j in range(element.p_order_v + 1):
                dof_id = f"{element.id}_{i}_{j}"
                self.dof_map[dof_id] = (element, i, j)

    def edge_dofs(self, element):
            # assign p orders per edge according to if the edge is horizontal or vertical,then use this to assign teh dofs accordingly
        for edge in element.edges[0:1]:
            # this is the u oriented edges
            dof_id = f"{element.id}_edge_{edge.id}"
            dof_count = element.p_order_u + 1
        for edge in element.edges[2:3]:
            # this is the v oriented edges
            dof_id = f"{element.id}_edge_{edge.id}"
            dof_count = element.p_order_v + 1
  

class global_dof_handling:
    def __init__(self,mesh):
        self.mesh = mesh
        self.global_dof_map = {}

    def assign_global_dofs(self):
        # traverse each elemetn and determine the shared local dofs on each edge between adjacent elements and assing global dof to that edge interface.
        for element in self.mesh.active_elements:
            for edge in element.edges:
                main_edge_dofs = edge.edge_dofs()
                # find the shared edge
                edge.adjacent_elem = [e for e in self.mesh.active_elements if edge in e.edges and e != element]
                if edge in edge.adjacent_elem.edges:
                    self.global_dof_map[edge.id] = min(main_edge_dofs, edge.adjacent_elem.edge_dofs())




# class integration:
#     def __init__(self, mesh, order_u, order_v):
#         self.mesh = mesh
#         self.order_u = order_u
#         self.order_v = order_v


#     def gaussian_quadrature(self):

class stiffness_matrix:
    def __init__(self, mesh, degrees_of_freedom):
        self.mesh = mesh
        self.degrees_of_freedom = degrees_of_freedom


    def assemble_stiffness_matrix(self):


    def local_stiffness_matrix(self, element):



class excitation_vector:
    def __init__(self, mesh)
        self.mesh = mesh

    def assemble_excitation_vector(self):


class boundary_conditions:
    def __init__(self, mesh, boundary_type, boundary_value):
        self.mesh = mesh
        self.boundary_type = boundary_type
        self.boundary_value = boundary_value



        
class FEM_solver:
    def __init__(self, stiffness, excitation, boundary, solution_type):
        self.stiffness = stiffness
        self.excitation = excitation
        self.boundary = boundary
        self.solution_type = solution_type

    def eigen_solution(self):



    def scattering_solution(self):







