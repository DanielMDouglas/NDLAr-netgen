# generate a mesh from an input step file

from ngsolve import *

from netgen.occ import *

def do_meshing(geo):
    # the meshing step is quite flexible
    # edit me!
    mesh = Mesh(geo.GenerateMesh(maxh=1.0))

    return mesh

def main(args):
    geo = OCCGeometry(args.infileName)

    mesh = do_meshing(geo)

    mesh.ngmesh.Export(args.outfileName, args.format)    

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description = 'use netgen to mesh a geometry (in step format) and export')
    parser.add_argument('infileName',
                        help = 'input step file')
    parser.add_argument('outfileName',
                        help = 'output mesh data')
    parser.add_argument('-f', '--format',
                        default = "Elmer Format",
                        help = 'output mesh format (default: Elmer Format)')

    args = parser.parse_args()

    main(args)
