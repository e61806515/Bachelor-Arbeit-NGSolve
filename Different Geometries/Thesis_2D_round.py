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

    inner = WorkPlane().Circle(r_i).Face()
    inner.edges.name="inner"

    rotor = WorkPlane().Circle(r_Fe).Face()
    rotor.maxh = maxh*delta
    rotor.name = "rotor"
    rotor.col = (1,0,0)

    if Tau != 1:
        magnets = [WorkPlane(Axes((0,0,0), n=Z, h=X))]*PZ
        d_phi = 360/PZ
        phi_M = d_phi*Tau/2
        for i in range(PZ):
            magnets[i].MoveTo(r_Fe*sin(i*d_phi*pi/180), r_Fe*cos(i*d_phi*pi/180))
            magnets[i].Arc(r_Fe, -phi_M).Rotate(90)
            magnets[i].Line(H_M).Rotate(90)
            magnets[i].Arc((r_Fe + H_M), 2*phi_M).Rotate(90)
            magnets[i].Line(H_M).Rotate(90)
            magnets[i].Arc(r_Fe, -phi_M).Rotate(-d_phi)
            magnets[i] = magnets[i].Face()
            magnets[i].maxh = 0.2e-3

            magnets[i].name = f"magnets_{i}"
            magnets[i].col = (0, 0, 1)


    else:
        magnets = WorkPlane().Circle(r_Fe+H_M).Face()
        magnets.name = "magnets"
        magnets.maxh = 0.5*delta
        magnets.col = (0,0,1)

        magnets = [magnets- rotor]


    outer = WorkPlane().Circle(r_L).Face()
    outer.name = "air"
    outer.edges.name = "outer"

    outer.col = (0,1,0)



    geo = Glue(([outer-rotor - sum(magnets), rotor- inner] + magnets))

    return geo


# -----------------------------------------------------------------------------
# ---------------------Parameters
# -----------------------------------------------------------------------------

PZ = 2

mu_air = 4e-7*np.pi
mu_magnet = 1.05*mu_air
mu_rotor = mu_air*5e2

sigma_magnet = 8e5
sigma_rotor =  1.86e6

order0 = 2
tau = 0.9

Br = 1
K0= 10000
f = 1000
omega = 2*np.pi*f

delta = lambda omega, sigma, mu : sqrt(2/(omega*sigma*mu))

if(sigma_rotor>1e-12):
    delta_rot = delta(omega, sigma_rotor, mu_rotor)
else:
     delta_rot = 5e-3
delta_mag = delta(omega, sigma_magnet, mu_magnet)

print(delta_rot)

r_Fe = 28e-3
mp = MeshingParameters(maxh=0.4)
mesh = Mesh(OCCGeometry(MakeGeometry(H_L=8e-3, H_M=6e-3, delta = delta_rot, r_Fe = r_Fe, Tau=tau, PZ=PZ, maxh=0.5), dim = 2).GenerateMesh(mp=mp))
mesh.Curve(3)

print(mesh.GetMaterials())
print(set(mesh.GetBoundaries()))
# -----------------------------------------------------------------------------
# ---------------------Coefficient Functions Sigma und Mu
# -----------------------------------------------------------------------------

mu = {"air":mu_air, "rotor": mu_rotor, "magnets.*": mu_magnet }
muCF = mesh.MaterialCF(mu, default=mu_air)



sigma = {"air":0, "rotor": sigma_rotor, "magnets.*": sigma_magnet}
sigmaCF = mesh.MaterialCF(sigma, default=0)


# Phi berechnung
def Phi(x,y):
    return atan2(y,x)

def K(x,y, nu=1):
     return K0*exp(-1j*nu*PZ*Phi(x,y)/2)

# -----------------------------------------------------------------------------
# ---------------------Finite Elemente Raum
# -----------------------------------------------------------------------------

V = H1(mesh, order = order0, dirichlet = "inner", complex=True)
trial = V.TrialFunction()
test = V.TestFunction()
u = GridFunction(V)

a = BilinearForm(V, symmetric = True)
a +=  1/muCF*grad(trial)*grad(test) * dx
a += 1j*omega*sigmaCF*test * trial * dx#("rotor|magnet|air") 1j*

c = Preconditioner(a, type="direct", inverse = "sparsecholesky")

f = LinearForm(V)
f += K(x,y)*test.Trace()*ds("outer") #*cos(2/D*x) PROBLEM weil hier periodische RBs


with TaskManager():

        solvers.BVP(bf=a, lf=f, gf=u, pre=c, needsassembling=True )

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
Draw(u*1j*omega*sigmaCF*mesh.MaterialCF({'magnets.*':1}, default=None), mesh, 'Jmag')
Draw(K(x,y), mesh, 'K')


# -----------------------------------------------------------------------------
# ---------------------WIRBELSTROMVERLUSTE
# -----------------------------------------------------------------------------
E = -1j * omega * u
J = sigmaCF * E

# E ... V/m
# J ... A/m^2

# W/m^3
# W/m


p = E*Conj(J)/2
losses = Integrate(p, mesh, definedon=mesh.Materials("magnets.*"))
#print("Int_K = ", Integrate(sin(Phi(x,y)), definedon=mesh.Boundaries[1], ))
print("P(u, u) = ", losses)
print("P/omega = ", losses/omega)

delta = lambda omega, sigma, mu : sqrt(2/(omega*sigma*mu))
# def delta(ome, s, m):
#      return sqrt(...)
delta_rot = delta(omega, sigma_rotor, mu_rotor)
delta_mag = delta(omega, sigma_magnet, mu_magnet)
# delta_rot = sqrt(2/(omega*1.86e6*mu_rotor))
# delta_mag = sqrt(2/(omega*0.667e6*mu_magnet))
print("delta_r = ", delta_rot)
print("delta_m = ", delta_mag)