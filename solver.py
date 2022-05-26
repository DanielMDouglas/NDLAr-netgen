# solve the Poisson equation -Delta u = f

from ngsolve import *

ngsglobals.msg_level = 0

from netgen.meshing import *

def do_solving(mesh):
    fes = H1(mesh, order=3, dirichlet=[3]) # pixel top

    # define trial- and test-functions
    u = fes.TrialFunction()
    v = fes.TestFunction()

    # the right hand side
    f = LinearForm(fes)
    # f += 32 * (y*(1-y)+x*(1-x)) * v * dx # add inhomogeneity (free charge)
    f += v*dx
    f += 0.5*v*ds(definedon='default') # neumann boundary conditions

    # the bilinear-form 
    a = BilinearForm(fes, symmetric=True)
    a += grad(u)*grad(v)*dx # this is the ordinary Poisson equation

    a.Assemble()
    f.Assemble()

    # define function for boundaries labelled dirichlet
    g = sin(x*pi)*sin(y*pi)
    # this can be an arbitrary function.  Go nuts!
    # g = 1
    # g = z
    # g = x*y*z
    # def g(x, y, z):
    #     if am_in_ring3(x,y,z):
    #         return 1
    #     else:
    #         return 0

    # the solution field 
    gfu = GridFunction(fes)
    gfu.Set(g, BND)

    c = Preconditioner(a, 'local')
    c.Update()
    solvers.BVP(bf=a, lf=f, gf=gfu, pre=c)

    return mesh, gfu

def main(args):
    mesh = comp.Mesh(ImportMesh(args.infileName))
    
    mesh, gfu = do_solving(mesh)

    # save the field to a vtk file with gfu labelled 'potential'
    vtk = VTKOutput(ma = mesh,
                    coefs = [gfu],
                    names = ['potential'],
                    filename = args.outfileName,
                    subdivision=3)
    vtk.Do()
    

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description = 'use netgen to mesh a geometry (in step format) and export')
    parser.add_argument('infileName',
                        help = 'input mesh data')
    parser.add_argument('outfileName',
                        help = 'output solution data')
    
    args = parser.parse_args()

    main(args)
