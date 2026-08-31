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



# technically these will only be used when doing 2D version of the algorithm. the 1d is all the same. 
    @property
    def tangent(self):

        x0 = np.asarray(self.nodes[0].coordinate)
        x1 = np.asarray(self.nodes[1].coordinate)

        t = x1 - x0

        return t / np.linalg.norm(t)

    @property
    def normal(self):
        # as in paper - rotate the Nedelec element (t) 90 degrees to get the normal element!!!
        t = self.tangent
        return np.array([t[1], -t[0]])


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
            child2 = element(child2_id, [new_node, parent.nodes[1]],  p_order = parent.p_order,level=level_new,parent=parent)
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


    def get_global_dof(self):
        neighbors = self.element.adjacent_list
        nodes = self.element.nodes
        # look into the different dofs and find the overlapping ones. in 1D this is literally just finding what nodes are what global dof
        globalList = []
        for i in range(len(neighbors)):
            nodes_neigh = neighbors[i].nodes
            for j in range(len(nodes)):
                globalList.append(nodes_neigh[i] if nodes_neigh[i] == nodes[i])

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
    


    


class dof_manager:
    def __init__(self, mesh):
        self.mesh = mesh

    def compute_total_dofs(self):




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


    ndof = dof_manager.num_dofs
    K_globe = np.zeros((ndof, ndof))
    #K_globe = np.zeros((len(nodes), len(nodes)))
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


 
def boundary(G, K, L, btype):
    #""'Neumann BC as specified in the problem statement. Dirichlet not yet defined, will need to do so'
    #A bit confused on implementing and derivign different BC. Think I am rusty
    G = G.copy()
    K = K.copy()
    if btype == 'dirichlet':
        G[0] = 1
        G[-1] = -np.cos(L)
    elif btype == 'neumann':
        G[0] = 1
        G[-1] = -np.cos(L)
    return G, K
 
 
def waveguide(K, Mmat):
    # This is not necessarily fully incorporated into the overally p-refining scheme
    # eigenproblem K x = lambda Mmat x 
    vals, vecs = scipy.linalg.eig(K, Mmat)
    return vals, vecs
 
 
def scattering(K, G,L):
    soln = scipy.linalg.solve(K, G)
    #error = error_scatter(soln,L)
    return soln
 
 
def HW(K, G,L):
    #soln = scipy.linalg.solve(K, G)
    #error = error_scatter(soln,L)
    soln = np.linalg.inv(K).dot(G)
    
    return soln




 
def local_errors(approx, exact, L, M,n):
    """Compute local errors for each element."""
    nodes, elements, jac, p_orders = mesh(L, M, n)
    local_errors = np.zeros(M)
    for e in range(M):
        elem = elements[e]
        approx_local = approx[elem]
        exact_local = exact[elem]
        local_errors[e] = np.sqrt(np.sum((approx_local - exact_local) ** 2))
    return local_errors


def p_refine(mesh_info, local_errors, threshold):
    sort_indices = np.argsort(local_errors,descending=True)
    sorted_errors = np.sort(local_errors,descending=True)
    top_to_refine = sort_indices[:len(sorted_errors)//3]
    refine = local_errors[top_to_refine]
    new_mesh_info = mesh_info
    mesh_elem = mesh_info[1]
    for i in range(len(refine)):
        if refine[i] > threshold:
            new_mesh_info[3][mesh_elem.index(top_to_refine[i])] += 1  
        else:
            break
    return new_mesh_info


def dual_weight_residual_global(solution, higher_order_soln):
    '''Computing DWR for the whole domain, using the solution and enriched solution to error estimate.
    In this case the operator L is self adjoint, so dual is same as primal, but the excitation is 
    then the divergence of the primal solution, giving the "charge" density to the primals E field.'''

    
    soln = solution
    enriched_soln = higher_order_soln
    interpolated_soln = np.interp(np.linspace(0, 1, len(enriched_soln)), np.linspace(0, 1, len(soln)), soln)
    error_estimate = np.gradient(enriched_soln) - np.gradient(interpolated_soln)
    inner_product = np.dot(enriched_soln, interpolated_soln)
    DWR = np.abs(inner_product) * error_estimate

    return DWR



def basis_scale_plotting(solution_in, L, M, mesh):
    '''From Cam Key code (adapted) scales the solution to the physical output!'''
    nodes, elems, jac, p_order, dof = mesh
    s = np.linspace(-1, 1, 100)
    plot_samps = np.zeros(M * 100)
    s_len = len(s)
    for e in range(M):
        elem = elems[e]
        for j_local in range(1, p_order[e] + 2):
            coeff = solution_in[elem[j_local - 1]]
            for k in range(s_len):
                plot_samps[e * s_len + k] += coeff * basis_eval(j_local, p_order[e], s[k])[0]
    return plot_samps




def p_refine_DWR_global(type, L, M, p_orders, g=0):
    # Use the DWR to refine the mesh with p refinement - no longer only globally
    p_orders = np.asarray(p_orders, dtype=int)
    low_order, low_order_plotting, mesher1 = main(p_orders, type, L, M, g)
    high_order, high_order_plotting, mesher2 = main(p_orders + 1, type, L, M, g)
    DWR = dual_weight_residual_global(low_order_plotting, high_order_plotting)
    error_indicators = DWR * high_order_plotting
    print(f'length of error: {len(error_indicators)}')
    return DWR, error_indicators, low_order_plotting, mesher1




# need to compute the low order initial basis funcitons across the global domain in the very first initial mesh. Then refinement is when we then superimpose the local p refinements 
# on top. the only dof for the low order basis is on the boundaries between elements so the local are not suhc a problem.save these in the mesh?? 






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



