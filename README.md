# NDLAr-netgen
NDLAr drift field simulation studies using netgen/ngsolve

This project is an effort to replace the cumbersome GUI methods currently in use for solving drift and weighting fields within DUNE NDLAr.

The field solving procedure has three main steps, which are defined in three scripts:

# Geometry Definition
The first step is defining the 3D geometry in a boundary representation.  This is usually done by constructing a geometry from primitives or from an external CAD software.  An example of this step is shown in `generate_pixelColumn_geometry.py`.  Here, the OCC interface in netgen is used to build a multi-body geometry from primitives which represents a drift volume very close to a pixel pad.

`python generate_pixelColumn_geometry.py`

# Meshing
The next step is to transform the boundary representation into a 3D unstructured grid for FEM calculations.  In the case of a multiple-body geometry, it's important that the input step file contains a compound solid (which can be done using the `Glue` function in netgen's OCC implementation, in the previous step).  The parameters of meshing are very flexible, and this step should be modified to suit the requirements of a particular problem.

`python mesher.py pixelColumn.step pixelColumn --format "Elmer Format"`

# Solving
Finally, the mesh can be loaded and a finite element space can be constructed from it.  The specification of a PDE to solve over the domain of the space is done using a weak formulation, and implemented as a bilinear form using the `u` and `v` trial and test functions.  Currently, the solver is set up to solve the ordinary Poisson equation (with no materials).  This program should also be modified to suit the users particular boundary conditions and charge density functions.  #this program is currently broken and will cease after one iteration!  Please use Elmer for now instead [or fix it :)]#

`python solver.py pixelColumn result`
