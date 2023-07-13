from netgen.occ import *
from ngsolve import *

#from netgen.webgui import Draw as DrawGeo
#from ngsolve.webgui import Draw
import numpy as np
import math as math

ngsglobals.msg_level = 5



#              Geometrie
#
#
#
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
    outer.edges.name = "outer"

    outer.col = (0,1,0)

    geo = Glue([outer-rotor - sum(magnets), rotor- inner] + magnets)
    #for i in range(PZ):
    #    geo = Glue([geo, magnets[i]])
    return geo

PZ = 4 
print(MakeGeometry(H_L=8e-3, H_M=6e-3, H_Fe=20e-3, H_W=30e-3, Tau=1, PZ=PZ))  
mp = MeshingParameters(maxh=0.01)
mp.
mesh = Mesh(OCCGeometry(MakeGeometry(H_L=8e-3, H_M=6e-3, H_Fe=20e-3, H_W=30e-3, Tau=1, PZ=PZ), dim = 2).GenerateMesh(maxh=0.01))
mesh.Curve(3)
#Materials
#('air', 'rotor', 'magnets_0', 'magnets_1', 'magnets_2', 'magnets_3', 'magnets_4', 'magnets_5', 'magnets_6', 'magnets_7', 'magnets_8', 'magnets_9', 'magnets_10', 'magnets_11')
#Boundaries
#{'default', 'inner', 'outer'}

print(mesh.GetMaterials())
print(set(mesh.GetBoundaries()))

#              Coefficient Functions Sigma und Mu
#
#
#
#version_0 = CF([0, 13, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
#version_1 = CF(list(range(len(mesh.GetMaterials()))))
mu_air = 4e-7*np.pi
mu_magnet = 10*mu_air
mu_rotor = mu_air*1000
mu = {"air":mu_air, "rotor": mu_rotor}
[mu.update({f"magnets_{i}": mu_magnet}) for i in range(PZ)]
#version_2 = CF([val[mat] if mat in val.keys() else 0 for mat in mesh.GetMaterials()])
muCF = mesh.MaterialCF(mu)

sigma = {"air":0, "rotor": 10e6}
[sigma.update({f"magnets_{i}": 10e6}) for i in range(PZ)]
sigmaCF = mesh.MaterialCF(sigma)
sigma_v = {"air":0, "rotor": 0}
[sigma_v.update({f"magnets_{i}": 10e3}) for i in range(PZ)]
sigma_visual = mesh.MaterialCF(sigma_v)

versions = [muCF, sigmaCF]
[Draw(versions[i], mesh, str(i)) for i in range(len(versions))]

Br = 1
K0= 10000
omega = 50

# Phi berechnung
def Phi(x,y):
    return atan2(y,x)

def K(x,y):
     K1 = -K0*sin(Phi(x,y))
     K2 = -K0*sin(Phi(x,y) + 2*pi/3)*exp(1j*2*pi/3*omega)
     K3 = -K0*sin(Phi(x,y) + 4*pi/3)*exp(1j*4*pi/3*omega)
     return K1 + K2 + K3     

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

print("P(u, u) = ", energy)
