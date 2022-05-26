# generate a geometry using the opencascade interface
# if called as a standalone script, save that geometry
# to a step file named 'pixelColumn.step'

from netgen.occ import *

pixelPitch = 4.434
driftDist = 12.
pixelWidth = 3.
pixelHeight = 0.4
pixelCornerRadius = 0.5

backplaneHeight = 1.
groundplaneHeight = 0.1

viaRadius = 0.5

def generate_pixelColumn_geo():
    # the LAr section for drifting
    LArShape = Box(Pnt(-pixelPitch/2, -pixelPitch/2, 0),
                   Pnt(pixelPitch/2, pixelPitch/2, driftDist))

    # the rounded corners of the pixel pad
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

    # the rectangular shapes joining the corners
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
    
    # remove the bit of LAr where the pixel protrudes
    LArShape -= pixelShape
    
    # the G10 backplane
    backplaneShape = Box(Pnt(-pixelPitch/2, -pixelPitch/2, -backplaneHeight),
                         Pnt(pixelPitch/2, pixelPitch/2, 0))
    # a via connecting the pixel pad to the ground plane
    via = Cylinder(Pnt(0, 0, -backplaneHeight),
                   Z,
                   r = viaRadius,
                   h = backplaneHeight)
    # add a hole where this via lives
    backplaneShape -= via
    
    # the ground plane behind the G10 backplane
    groundplaneShape = Box(Pnt(-pixelPitch/2, -pixelPitch/2, -(backplaneHeight + groundplaneHeight)),
                           Pnt(pixelPitch/2, pixelPitch/2, -backplaneHeight))
    
    # make a compound solid from these shapes
    geo = Glue([pixelShape, LArShape, backplaneShape, groundplaneShape, via])

    return geo

def main(args):
    geo = generate_pixelColumn_geo()
    geo.WriteStep(args.outfileName)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description = 'use netgen to mesh a geometry (in step format) and export')
    parser.add_argument('-o', '--outfileName',
                        default = 'pixelColumn.step',
                        help = 'output step file name (default: pixelColumn.step)')
    
    args = parser.parse_args()

    main(args)
