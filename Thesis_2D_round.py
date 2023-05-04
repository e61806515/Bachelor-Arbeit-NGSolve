from netgen.occ import *
from ngsolve import *

from netgen.webgui import Draw as DrawGeo
from ngsolve.webgui import Draw
import numpy as np
import math as math

def MakeGeometry(H_L, H_M, H_Fe, H_W, Tau, PZ):
    
    r_Fe = H_Fe + H_W
    r_L = r_Fe + H_L
    inner = WorkPlane().Circle(H_W).Face()
    inner.edges.name="inner"

    rotor = WorkPlane().Circle(r_Fe).Face()
    rotor.name = "rotor"
    rotor.col = (1,0,0)

    magnets = [WorkPlane(Axes((0,0,0), n=Z, h=X))]*PZ
    d_phi = 360/PZ
    phi_M = d_phi*Tau/2
    for i in range(PZ):
        #magnets[i] = WorkPlane(Axes((0,0,0), n=Z, h=X))
        magnets[i].MoveTo(r_Fe*sin(i*d_phi*pi/180), r_Fe*cos(i*d_phi*pi/180))
        magnets[i].Arc(r_Fe, -phi_M).Rotate(90)
        magnets[i].Line(H_M).Rotate(90)
        magnets[i].Arc((r_Fe + H_M), 2*phi_M).Rotate(90)
        magnets[i].Line(H_M).Rotate(90)
        magnets[i].Arc(r_Fe, -phi_M).Rotate(-d_phi)
        magnets[i] = magnets[i].Face()
        magnets[i].name = f"magnets_{i}" #{i:"Anweisungen zur Formatierung"} Formated String
        magnets[i].col = (0, 0, 1)

    outer = WorkPlane().Circle(r_L).Face()
    outer.name = "air"

    air = (outer - rotor)
    air.name = "air"
    outer.edges.name = "outer"

    outer.col = (0,1,0)

    geo = Glue([outer-rotor- magnets[0] -magnets[1], rotor- inner, magnets[0], magnets[1]])
    #for i in range(PZ):
    #    geo = Glue([geo, magnets[i]])
    return geo

mesh = Mesh(OCCGeometry(MakeGeometry(H_L=7e-3, H_M=6e-3, H_Fe=20e-3, H_W=30e-3, Tau=5/6, PZ=12), dim = 2).GenerateMesh(maxh=1))
print(mesh.GetMaterials())
mesh.Curve(3)

Draw(mesh)

