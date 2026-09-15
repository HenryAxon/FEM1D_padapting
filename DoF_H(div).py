# this file will hopefully implement the DoF details necesary for a basic RT/H(div) 2D solver.




class vector_basis_function:
    '''
    Computes basis (need to add derivatives) 
    '''
    def __init__(self, mesh, order_u, order_v):
        self.mesh = mesh
        self.order_u = order_u
        self.order_v = order_v
        self.basis_functions = {}

    def basis_functions(self,u, v):
        #Notaros/Ilic higher order basis functions
        # 
        if self.order_v == 0: 
            v_basis = (1 - v) * u**self.order_u
        elif self.order_v == 1:
            v_basis = (v + 1) * u **self.order_u
        elif self.order_v >= 2 and self.order_v % 2 == 0:
            v_basis = (v**self.order_v - 1) * u **self.order_u
        else:
            v_basis = (v**self.order_v - v) * u **self.order_u

        # now the cases where we have the u direction basis functions

        if self.order_u == 0:
            u_basis = (1 - u) * v**self.order_v
        elif self.order_u == 1:
            u_basis = (u + 1) * v **self.order_v
        elif self.order_u >= 2 and self.order_u % 2 == 0:
            u_basis = (u**order_u - 1) * v **self.order_v 
        else:
            u_basis = (u**self.order_u - u) * v **self.order_v

        return u_basis, v_basis

    def derivatives(self,u,v):
         # this is really the divergence I suppose if you add the resulting components
         if order_v == 0: 
            v_basis = (- v) * u ** order_u
        elif order_v == 1:
            v_basis = (v) * u **order_u
        elif order_v >= 2 and order_v % 2 == 0:
            v_basis = ((order_v) * v**(order_v-1)) * u **self.order_u
        else:
            v_basis = ((order_v) * v**(order_v-1) - 1) * u **order_u

        # now the cases where we have the u direction basis functions

        if order_u == 0:
            u_basis = (- u) * v**order_v
        elif order_u == 1:
            u_basis = (u) * v **order_v
        elif order_u >= 2 and order_u % 2 == 0:
            u_basis = ((order_u)*u**(order_u-1)) * v **order_v 
        else:
            u_basis = ((order_u)*u**(order_u-1) - 1) * v **order_v

        return du_basis, dv_basis   

         

    def contravariant_piola(self):
         # Piola transform for a 2D quad element, contravariant maps physical to reference while preserving flux
         # for affine quads I think I can just worry about mapping the nodes on each edge. and taking jacobians for integration.
        for k in self.mesh.elements:
            edges = [i for i in k.edges]
            nodes = [j.nodes for j in edges]

            top_node_y = max(nodes[1])
            bottom_y = min(nodes[1])
            right_x = max(nodes[0])
            left_x = min(nodes[0])

            map_y = (top_node_y - bottom_y) # transforms the quad to reference square
            map_x = right_x - left_x

        


    def standard_piola(self, u, v, mesh):
         # this version is the standard transformation mapping the reference to physical u->x, v->y

         
             

  

class global_dof_handling:
    def __init__(self,mesh):
        self.mesh = mesh
        self.global_dof_map = {}

    def assign_global_dofs(self):
        # traverse each element and determine the shared local dofs on each edge between adjacent elements and assing global dof to that edge interface.
        for element in self.mesh.active_elements:
            for edge in element.edges:
                main_edge_dofs = edge.edge_dofs()
                # find the shared edge
                edge.adjacent_elem = [e for e in self.mesh.active_elements if edge in e.edges and e != element]
                if edge in edge.adjacent_elem.edges:
                    self.global_dof_map[edge.id] = min(main_edge_dofs, edge.adjacent_elem.edge_dofs())

class local:
    def __init__(self,mesh):
        self.mesh = mesh


    def local_dofs(self, element):
        # local DoF computation on an element, computes the normal trace aka just the H(div) hierarchical basis 
        # for each edge, and then the inner degrees of freedom.
        dofs = []

        # Edge DOFs
        for edge_item in element.edges:
            if edge_item in element.edges[:2]:
                order = element.p_order_u
            else:
                order = element.p_order_v

            for index in range(order + 1):
                dofs.append(("edge", edge_item.id, index, "normal"))

        # Interior DOFs
        for i in range(element.p_order_u + 1):
            for j in range(element.p_order_v + 1):
                dofs.append(("element", element.id, i, j,"none"))

        return dofs

class deactivate_dof:
    def __init__(self, mesh, dofs):
         self.mesh = mesh
         self.dofs = dofs

    def inner_dof(self):
         for i in self.mesh.elements:
            if i.is_leaf == False:
                for j in self.dofs[i.id]
                    if j[0] == "element":
                        j = 0 

    def edge_dof(self):
         # checks for each leaf that i
        for i in self.mesh.elements:
            for k in i.adjacent_elems
                if k.level == i.level:
                    pass
                elif k.level < i.level
                    for j in self.dofs[i.id]:
                        if j[0] == "edge"
                            j = 0
                else:
                     pass


                      
                
         
