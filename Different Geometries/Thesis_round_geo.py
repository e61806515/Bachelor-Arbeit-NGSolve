from netgen.occ import *
from ngsolve import *

#from netgen.webgui import Draw as DrawGeo
#from ngsolve.webgui import Draw
from netgen.meshing import MeshingParameters
import numpy as np
import math as math
import netgen.gui

#ngsglobals.msg_level = 5



#              Geometrie
#
#
#
def MakeGeometry(H_L, H_M, delta, r_Fe, Tau, PZ, maxh):

    print("Tau = ", Tau)
    d_Fe = 5*delta
    r_i = r_Fe - d_Fe
    r_L = r_Fe + H_L

    rect = WorkPlane(Axes((-r_L,-r_L,0), n=Z, h=X)).Rectangle(2*r_L, r_L).Face()

    inner = WorkPlane().Circle(r_i).Face() #- rect
    inner.edges.name="inner"

    rotor = WorkPlane().Circle(r_Fe).Face() #-rect
    rotor.maxh = maxh*delta
    rotor.name = "rotor"
    rotor.col = (1,0,0)

    if Tau != 1:
        magnets = [WorkPlane(Axes((0,0,0), n=Z, h=X))]*PZ
        d_phi = 360/PZ
        phi_M = d_phi*Tau/2
        for i in range(PZ):
            #if (i*d_phi < 90) or (i*d_phi > 270):
                #magnets = magnets + [WorkPlane(Axes((0,0,0), n=Z, h=X))]
                magnets[i].MoveTo(r_Fe*sin(i*d_phi*pi/180), r_Fe*cos(i*d_phi*pi/180))
                magnets[i].Arc(r_Fe, -phi_M).Rotate(90)
                magnets[i].Line(H_M).Rotate(90)
                magnets[i].Arc((r_Fe + H_M), 2*phi_M).Rotate(90)
                magnets[i].Line(H_M).Rotate(90)
                magnets[i].Arc(r_Fe, -phi_M).Rotate(-d_phi)
                magnets[i] = magnets[i].Face()
                magnets[i].maxh = 0.2e-3
                #magnets[i].vertices.maxh = 0.00001 #geht nicht??
                magnets[i].name = f"magnets_{i}" #{i:"Anweisungen zur Formatierung"} Formated String
                magnets[i].col = (0, 0, 1)


    else:
        magnets = WorkPlane().Circle(r_Fe+H_M).Face() #-rect
        magnets.name = "magnets"
        magnets.maxh = 0.5*delta
        magnets.col = (0,0,1)

        magnets = [magnets- rotor]


    outer = WorkPlane().Circle(r_L).Face() #- rect
    outer.name = "air"
    outer.edges.name = "outer"

    outer.col = (0,1,0)



    geo = Glue(([outer-rotor - sum(magnets), rotor- inner] + magnets))

    return geo


# -----------------------------------------------------------------------------
# --- Parameters
# -----------------------------------------------------------------------------

PZ = 10

mu_air = 4e-7*np.pi
mu_magnet = 1.05*mu_air
mu_rotor = mu_air*5e2

sigma_magnet = 8e5
sigma_rotor =  1.86e6

order0 = 2
tau = 3/4

Br = 1
K0= 10000
f = 100
omega = 2*np.pi*f

delta = lambda omega, sigma, mu : sqrt(2/(omega*sigma*mu))

if(sigma_rotor>1e-12):
    delta_rot = delta(omega, sigma_rotor, mu_rotor)
else:
     delta_rot = 5e-3
delta_mag = delta(omega, sigma_magnet, mu_magnet)

print(delta_rot)

r_Fe = 28e-3
#print(MakeGeometry(H_L=8e-3, H_M=6e-3, H_Fe=20e-3, H_W=30e-3, Tau=1, PZ=PZ))
mp = MeshingParameters(maxh=0.4)
#mp.RestrictH(x=0, y=0, z=1, h=0.0025)
mesh = Mesh(OCCGeometry(MakeGeometry(H_L=8e-3, H_M=6e-3, delta = 0.3e-3, r_Fe = r_Fe, Tau=tau, PZ=PZ,maxh = 0.5), dim = 2).GenerateMesh(mp=mp))
mesh.Curve(3)
#Materials
#('air', 'rotor', 'magnets_0', 'magnets_1', 'magnets_2', 'magnets_3', 'magnets_4', 'magnets_5', 'magnets_6', 'magnets_7', 'magnets_8', 'magnets_9', 'magnets_10', 'magnets_11')
#Boundaries
#{'default', 'inner', 'outer'}

print(mesh.GetMaterials())
print(set(mesh.GetBoundaries()))
regions = mesh.MaterialCF({"air":1, "rotor": 2, "magnets.*": 3})
Draw(regions, mesh, "regions")
# -----------------------------------------------------------------------------
#              Coefficient Functions Sigma und Mu
#
#                   magnet:NdFeB mu_r = 1.04 - 1.08 und sigma = 0.667e6
#                   rotor: Fe mu_r = 1e4 oder 5e3 u. sigma = 1.76e6 - 1.97e6 .... Welche Werte Herr Schmid??????
#
# -----------------------------------------------------------------------------

#version_0 = CF([0, 13, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
#version_1 = CF(list(range(len(mesh.GetMaterials()))))


mu = {"air":mu_air, "rotor": mu_rotor, "magnets.*": mu_magnet }
# [mu.update({f"magnets_{i}": mu_magnet}) for i in range(PZ)]
# mu.update({f"magnets_{i}": mu_magnet for i in range(PZ)})
#version_2 = CF([val[mat] if mat in val.keys() else 0 for mat in mesh.GetMaterials()])
muCF = mesh.MaterialCF(mu, default=mu_air)



sigma = {"air":0, "rotor": sigma_rotor, "magnets.*": sigma_magnet}
sigmaCF = mesh.MaterialCF(sigma, default=0)

# [sigma.update({f"magnets_{i}":sigma_magnet}) for i in range(PZ)]
# sigma.update({f"magnets.*": sigma_magnet })

#sigma = {"air":0, "rotor": 0}
#[sigma.update({f"magnets_{i}": 0}) for i in range(PZ)]
# sigma_v = {"air":None, "rotor": None}
# [sigma_v.update({f"magnets_{i}": 8e5}) for i in range(PZ)]
# sigma_v = {"air":0, "rotor": None, "magnets.*": sigma_magnet}
# sigma_visual = mesh.MaterialCF(sigma_v, default=None)

# versions = [muCF, sigmaCF]
# [Draw(versions[i], mesh, str(i)) for i in range(len(versions))]



