"""Eddy Current Problem in 2D round geometry.
"""
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
def MakeGeometry(H_L, H_M, H_Fe, H_W, Tau, PZ):

    print("Tau = ", Tau)

    # Radien des Rotors und der ganzen Geometrie
    r_Fe = H_Fe + H_W
    r_L = r_Fe + H_L

    rotor = WorkPlane().Circle(r_Fe).Face() #-rect
    rotor.name = "rotor"
    rotor.col = (1,0,0)

    outer = WorkPlane().Circle(r_L).Face() #- rect
    outer.name = "air"
    outer.edges.name = "outer"
    inner = WorkPlane().Circle(H_W).Face() #- rect
    inner.edges.name="inner"
    outer.col = (0,1,0)

    if False:
        magnets = WorkPlane().Circle(r_Fe+H_M).Face()
        magnets = magnets - rotor
        magnets.name = "magnets"
        geo = Glue(([outer -rotor- magnets, rotor- inner] + magnets))
    else:
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
                magnets[i].maxh = 0.001
                #magnets[i].vertices.maxh = 0.00001 #geht nicht??
                magnets[i].name = f"magnets_{i}" #{i:"Anweisungen zur Formatierung"} Formated String
                magnets[i].col = (0, 0, 1)
        geo = Glue(([outer-rotor - sum(magnets), rotor- inner] + magnets))
    #geo = Glue([geo - rect])
    #for i in range(PZ):
    #    geo = Glue([geo, magnets[i]])
    geo.edges.maxh = 0.0005
    return geo

PZ = 2
#print(MakeGeometry(H_L=8e-3, H_M=6e-3, H_Fe=20e-3, H_W=30e-3, Tau=1, PZ=PZ))
mp = MeshingParameters(maxh=0.4)
#mp.RestrictH(x=0, y=0, z=1, h=0.0025)
mesh = Mesh(OCCGeometry(MakeGeometry(H_L=8e-3, H_M=6e-3, H_Fe=25e-3, H_W=3e-3, Tau=19/20, PZ=PZ), dim = 2).GenerateMesh(mp=mp))
mesh.Curve(3)
#Materials
#('air', 'rotor', 'magnets_0', 'magnets_1', 'magnets_2', 'magnets_3', 'magnets_4', 'magnets_5', 'magnets_6', 'magnets_7', 'magnets_8', 'magnets_9', 'magnets_10', 'magnets_11')
#Boundaries
#{'default', 'inner', 'outer'}

print(mesh.GetMaterials())
print(set(mesh.GetBoundaries()))

#              Coefficient Functions Sigma und Mu
#
#                   magnet:NdFeB mu_r = 1.04 - 1.08 und sigma = 0.667e6
#                   rotor: Fe mu_r = 1e4 oder 5e3 u. sigma = 1.76e6 - 1.97e6 .... Welche Werte Herr Schmid??????
#
#version_0 = CF([0, 13, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
#version_1 = CF(list(range(len(mesh.GetMaterials()))))
mu_air = 4e-7*np.pi
mu_magnet = 1.05*mu_air
mu_rotor = mu_air*5e2
mu = {"air":mu_air, "rotor": mu_rotor}
[mu.update({f"magnets_{i}": mu_magnet}) for i in range(PZ)]
#version_2 = CF([val[mat] if mat in val.keys() else 0 for mat in mesh.GetMaterials()])
muCF = mesh.MaterialCF(mu, default=mu_air)

sigma = {"air":0, "rotor": 1.86e-15}
[sigma.update({f"magnets_{i}": 8e5}) for i in range(PZ)]
#sigma = {"air":0, "rotor": 0}
#[sigma.update({f"magnets_{i}": 0}) for i in range(PZ)]
sigmaCF = mesh.MaterialCF(sigma, default=0)
sigma_v = {"air":None, "rotor": None}
[sigma_v.update({f"magnets_{i}": 8e5}) for i in range(PZ)]
sigma_visual = mesh.MaterialCF(sigma_v)


versions = [muCF, sigmaCF]
[Draw(versions[i], mesh, str(i)) for i in range(len(versions))]


Br = 1
K0= 10000
f = 50
omega = 2*np.pi*f

# Phi berechnung
def Phi(x,y):
    return atan2(y,x)

def K(x,y, nu=1):
    K = K0*exp(-1j*nu*Phi(x,y))
    if(False):
        K1 = -K0*sin(Phi(x,y))
        K2 = -K0*sin(Phi(x,y) + 2*pi/3)*exp(1j*2*pi/3)
        K3 = -K0*sin(Phi(x,y) + 4*pi/3)*exp(1j*4*pi/3)
        K = K1 + K2 + K3
    if(False):
        K4 = -K0*sin(Phi(x,y) + pi/2)
        K5 = -K0*sin(Phi(x,y) + 2*pi/3 + pi/2)*exp(1j*2*pi/3)
        K6 = -K0*sin(Phi(x,y) + 4*pi/3 + pi/2)*exp(1j*4*pi/3)
        K = K + K4 + K5 + K6
    return K
#              Finite Elemente Raum
#
#
#

V = H1(mesh, order = 2, dirichlet = "inner", complex=True)
trial = V.TrialFunction()
test = V.TestFunction()
u = GridFunction(V)

a = BilinearForm(V, symmetric = True)
a +=  1/muCF*grad(trial)*grad(test) * dx
a += 1j*omega*sigmaCF*test * trial * dx#("rotor|magnet|air") 1j*

c = Preconditioner(a, type="direct", inverse = "sparsecholesky")

f = LinearForm(V)
f += K(x,y)*test.Trace()*ds("outer") #*cos(2/D*x) PROBLEM weil hier periodische RBs

#u.Set(coef_dirichlet, BND)
#solver

with TaskManager():
        a.Assemble()
        f.Assemble()
        bvp = BVP(bf=a, lf=f, gf=u, pre=c)
        bvp.Do()

#u.vec.data = a.mat.Inverse() * f.vec #freedofs=V.FreeDofs()
#print("a lautet ", a.mat)
#print("f lautet ", f.vec)
#print("DG lautet ", u.vec.data)
B = CF((grad(u)[1], -grad(u)[0]))       #Gradient(Komponenten) sind L2-Funktionen. Grad ist nur 2-dim,
                                        #weil Geometrie nur 2-dim
Draw(u) #u vom Typ gridfunction - Information über mesh bereits implizit enthalten
Draw(B, mesh, 'B') #B vom Typ tuple, keine Information über mesh
Draw(1/muCF*B, mesh, 'H')
Draw(Norm(1/muCF*B[0]), mesh, 'Norm Hx')
Draw(Norm(1/muCF*B[1]), mesh, 'Norm Hy')
Draw(Norm(B[0]), mesh, 'Norm Bx')
Draw(Norm(B[1]), mesh, 'Norm By')
Draw(u*1j*omega*sigmaCF, mesh, 'J')
Draw(K(x,y), mesh, 'K')
#   WIRBELSTROMVERLUSTE
#
#
p = sigma_visual*omega*omega*u*Conj(u)/2
energy = Integrate(p, mesh)
#print("Int_K = ", Integrate(sin(Phi(x,y)), definedon=mesh.Boundaries[1], ))
print("P(u, u) = ", energy)
print("P/omega = ", energy/omega)


delta_rot = sqrt(2/(omega*1.86e6*mu_rotor))
delta_mag = sqrt(2/(omega*0.667e6*mu_magnet))
print("delta_r = ", delta_rot)
print("delta_m = ", delta_mag)