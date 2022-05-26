# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0

from ngsolve import *

ngsglobals.msg_level = 0

# from netgen.csg import *
from netgen.occ import *

pixelPitch = 4.434
driftDist = 12.
pixelWidth = 3.
pixelHeight = 0.4
pixelCornerRadius = 0.5

backplaneHeight = 1.
groundplaneHeight = 0.1

viaRadius = 0.5

LArShape = Box(Pnt(-pixelPitch/2, -pixelPitch/2, 0),
               Pnt(pixelPitch/2, pixelPitch/2, driftDist))

cornerCyls = [Cylinder(Pnt(-pixelWidth/2 + pixelCornerRadius,
                           -pixelWidth/2 + pixelCornerRadius,
                           0),
                       Z,
                       r = pixelCornerRadius,
                       h = pixelHeight),
              Cylinder(Pnt(pixelWidth/2 - pixelCornerRadius,
                           -pixelWidth/2 + pixelCornerRadius,
                           0),
                       Z,
                       r = pixelCornerRadius,
                       h = pixelHeight),
              Cylinder(Pnt(-pixelWidth/2 + pixelCornerRadius,
                           pixelWidth/2 - pixelCornerRadius,
                           0),
                       Z,
                       r = pixelCornerRadius,
                       h = pixelHeight),
              Cylinder(Pnt(pixelWidth/2 - pixelCornerRadius,
                           pixelWidth/2 - pixelCornerRadius,
                           0),
                       Z,
                       r = pixelCornerRadius,
                       h = pixelHeight)]

rectA = Box(Pnt(-pixelWidth/2 + pixelCornerRadius,
                -pixelWidth/2,
                0),
            Pnt(pixelWidth/2 - pixelCornerRadius,
                pixelWidth/2,
                pixelHeight))
rectB = Box(Pnt(-pixelWidth/2,
                -pixelWidth/2 + pixelCornerRadius,
                0),
            Pnt(pixelWidth/2,
                pixelWidth/2 - pixelCornerRadius,
                pixelHeight))

pixelShape = sum(cornerCyls)
pixelShape *= Box(Pnt(-pixelPitch/2, -pixelPitch/2, 0),
                  Pnt(pixelPitch/2, pixelPitch/2, pixelHeight))
pixelShape += rectA
pixelShape += rectB

LArShape -= pixelShape

backplaneShape = Box(Pnt(-pixelPitch/2, -pixelPitch/2, -backplaneHeight),
                     Pnt(pixelPitch/2, pixelPitch/2, 0))

via = Cylinder(Pnt(0, 0, -backplaneHeight),
               Z,
               r = viaRadius,
               h = backplaneHeight)

backplaneShape -= via

groundplaneShape = Box(Pnt(-pixelPitch/2, -pixelPitch/2, -(backplaneHeight + groundplaneHeight)),
                       Pnt(pixelPitch/2, pixelPitch/2, -backplaneHeight))

# geo = CSGeometry()

# geo.Add (driftVolume)

# geo = OCCGeometry(LArShape)
geo = OCCGeometry(Glue([pixelShape, LArShape, backplaneShape, groundplaneShape, via]))

mesh = Mesh(geo.GenerateMesh(maxh=1.0))
print (mesh.GetBoundaries())
# print (mesh.GetBodies())
# mesh.Save("cube.vol")

# H1-conforming finite element space
# fes = H1(mesh, order=3, dirichlet=[1,2,3,4])
print ('making fes')
fes = H1(mesh, order=3, dirichlet=[3]) # seems to be pixel top
# fes = H1(mesh, order=3, dirichlet=[14]) # to drift volume
print ('made fes')

# define trial- and test-functions
u = fes.TrialFunction()
v = fes.TestFunction()

# the right hand side
f = LinearForm(fes)
# f += 32 * (y*(1-y)+x*(1-x)) * v * dx
f += v*dx
f += 0.5*v*ds(definedon='default') # neumann boundary conditions

# the bilinear-form 
a = BilinearForm(fes, symmetric=True)
# a += grad(u)*grad(v)*dx
a += grad(u)*grad(v)*dx

a.Assemble()
f.Assemble()

g = sin(x*pi/pixelPitch)*sin(y*pi/pixelPitch)
# g = 1
# g = z
# g = x*y*z
# def g(x, y, z):
#     if am_in_ring3(x,y,z):
        
        
# the solution field 
gfu = GridFunction(fes)
gfu.Set(g, BND)

c = Preconditioner(a, 'local')
c.Update()
solvers.BVP(bf=a, lf=f, gf=gfu, pre=c)

# r = f.vec.CreateVector()
# r.data = f.vec - a.mat*gfu.vec

# # gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
# gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r
# # print (u.vec)


# plot the solution (netgen-gui only)
# Draw (gfu)
# Draw (-grad(gfu), mesh, "Flux")

# exact = 16*x*(1-x)*y*(1-y)
# print ("L2-error:", sqrt (Integrate ( (gfu-exact)*(gfu-exact), mesh)))

# gfu.Save('thing.vtk')

vtk = VTKOutput(ma = mesh,
                coefs = [gfu],
                names = ['potential'],
                filename = 'result',
                subdivision=3)
vtk.Do()
