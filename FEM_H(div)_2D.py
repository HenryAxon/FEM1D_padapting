# main script to run the FEM code


# 1) first need to generate initial mesh


# 2) remember to collect DoF in a suitable way and compute the basis functions on each element, 


# 3) compute the initial solution for unrefined mesh, fill FEM matrices for general impedences and excitations etc


# 4) solve for the error indicator - via global refinement in p to get adjoint error indicator? 


# 5) mark elements for h, p refinement types, and refine


# 6) collect DoFs and activate/deactivate parents nad children DoFs and geoemtric components


# 7) loop steps 4 - 6 until the absolute error across the whole system is less than a predtermined threshold, or the number of iterations is exceeded

