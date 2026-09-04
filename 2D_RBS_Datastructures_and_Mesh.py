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
    def __init__(self, id, coordinatex,coordinatey, level=0, parent=None):
        super().__init__(id, active=True, level=level, parent=parent)
        self.coordinatex = coordinatex
        self.coordinatey = coordinatey
        self.adjacent_elem = []
  

    
class edge(topology):
    def __init__(self, id, nodes, level=0, parent=None):
        super().__init__(id, active=True, level=level, parent=parent)
        self.nodes = nodes
        self.adjacency = []

class face(topology):
    def __init__(self, id, nodes,edges, level, parent):
        super().__init__(id, active=True, level=level, parent=parent)
        self.edges = edges
        self.nodes = nodes
        self.adjacent = []



#class face(topology):


class element(topology):
    def __init__(self, id,nodes,edges,faces,p_order, level, parent):
        super().__init__(id=id, active=True, level=level, parent=parent)
        self.p_order = p_order
        self.nodes = nodes
        self.edges = edges
        self.faces = faces
        self.adjacency_list = []

    def deactivate(self):
        if self.active == False:
            p_order = 0 
        return p_order

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
    def __init__(self,L, N, p, level=0):
        self.L = L
        self.N = N
        self.p = p
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
                e=element(self.next_elem_id(), nodes=[self.node_set[node_index(m,n)],self.node_set[node_index(m,n+1)], self.node_set[node_index(m+1,n)],self.node_set[node_index(m+1,n+1)]], edges=[self.edges_u[edge_u_index(m,n)],self.edges_u[edge_u_index(m+1,n)],self.edges_v[edge_v_index(m,n)],self.edges_v[edge_v_index(m,n+1)]],faces = [], p_order=self.p,level=self.level,parent=None)
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
                node_item.p_order = 0
        for edge in self.edges:
            adjacent = edge.adjacent_elem
            if adjacent.level < edge.level and edge in adjacent.edges:
                edge.active = False
                edge.p_order = 0
    
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
            parent.p_order = 0

    def refine_p(self):
        
        for i in self.marked_elem_p:
            self.mesh.elems[i].p_order +=1
        return self.mesh


mesh1 = mesh(L=10, N=6, p=2)
mesh1.gen_init_mesh()
refined_1 = refiner(mesh1, marked_elem_h_t=[0,1], marked_elem_h_u=[], marked_elem_h_v=[], marked_elem_p=[2])
refined_1.refine_h()
refined_1.refine_p()


ax1 = refined_1.mesh.plot_mesh(annotate=True, show_nodes=True)
plt.savefig("mesh_preview.png", dpi=150)
print("saved mesh_preview.png")
# in 1d these are essentially the same, but writing the skeleton can still be useful for the 2D implementation


class func_space:
    def __init__(self, element):
        self.element = element

    def get_local_dofs(self, element):
        dofs = []

        # vertex DOFs
        for node in element.nodes:
            dofs.append(node)

        # interior DOFs
        for i in range(element.p_order - 1):
            dofs.append((element, i))

        return dofs


    # def get_global_dof(self):
    #     neighbors = self.element.adjacent_list
    #     nodes = self.element.nodes
    #     # look into the different dofs and find the overlapping ones. in 1D this is literally just finding what nodes are what global dof
    #     globalList = []
    #     for i in range(len(neighbors)):
    #         nodes_neigh = neighbors[i].nodes
    #         for j in range(len(nodes)):
    #             globalList.append(nodes_neigh[i] if nodes_neigh[i] == nodes[i])

class DOFManager:

    def __init__(self):
        self.global_dofs = {}
        self.next_id = 0

    def get_global_dof(self, local_dof):

        if local_dof not in self.global_dofs:
            self.global_dofs[local_dof] = self.next_id
            self.next_id += 1

        return self.global_dofs[local_dof]

class H1_cont(func_space):
    
    def basis(self, element, s):

        p = element.p_order

        return [
            basis_eval(j, p, s)[0]
            for j in range(1, p + 2)
        ]

    def local_dofs(self, element):

        # vertex DOFs
        dofs = [
            ("node", element.nodes[0].id),
            ("node", element.nodes[1].id)
        ]

        # element-interior DOFs
        for j in range(1, element.p_order):
            dofs.append(("element", element.id, j))

        return dofs

class raviart_thomas_cont(func_space):
    #associate tangential component on geometric component
    def enforce_curl(self,element,s):
        p = element.p_order
        basis = []




        return tans


class nedelec_cont(func_space):
    # "associate" normal component on geometric component
    def enforce_div(self, element, s):
        p = element.p_order
        basis = []
        return divs 
    


# analytic solution is sin(x), so we should see that here.
#print(f'solutions: {solutionstore}')



