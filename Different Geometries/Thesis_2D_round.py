from netgen.occ import *
from ngsolve import *

#from netgen.webgui import Draw as DrawGeo
#from ngsolve.webgui import Draw
from netgen.meshing import MeshingParameters, IdentificationType
import numpy as np
import math as math
import netgen.gui

#ngsglobals.msg_level = 5



#              Geometrie
#
#
#

def MakeGeometry(H_L, H_M, maxh_rotor, maxh_mag, r_Fe, tau, PZ, maxh):
    print("maxh = ", maxh)
    print("tau = ", tau)
    print("delta_mag = ", delta_mag)
    print("maxh_mag = ", maxh * delta_mag)
    d_Fe = 8*maxh_rotor
    r_i = r_Fe - d_Fe
    r_L = r_Fe + H_L
    print("r_Fe = ", r_Fe)
    print("r_L = ", r_L)
    print("r_M = ", r_Fe + H_M)
    print("r_i = ", r_i)
    rect = WorkPlane(Axes((-r_L,-r_L,0), n=Z, h=X)).Rectangle(2*r_L, r_L).Face()

    inner = WorkPlane().Circle(r_i).Face()
    inner.edges.name="inner"

    rotor = WorkPlane().Circle(r_Fe).Face() #-rect
    rotor.maxh = maxh_rotor
    rotor.name = "rotor"
    rotor.col = (1,0,0)

    if tau != 1:
        magnets = [WorkPlane(Axes((0,0,0), n=Z, h=X))]*PZ
        d_phi = 360/PZ
        phi_M = d_phi*tau/2
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
                magnets[i].maxh = maxh_mag
                #magnets[i].vertices.maxh = 0.00001 #geht nicht??
                magnets[i].name = f"magnets_{i}" #{i:"Anweisungen zur Formatierung"} Formated String
                magnets[i].col = (0, 0, 1)


    else:
        magnets = WorkPlane().Circle(r_Fe+H_M)
        magnets = magnets.Face()
        magnets.name = "magnets"
        magnets.maxh = maxh_mag
        magnets.col = (0,0,1)

        magnets = [magnets- rotor]


    outer = WorkPlane().Circle(r_L).Face()
    outer.name = "air"
    outer.maxh = maxh_mag
    outer.edges.name = "outer"

    outer.col = (0,1,0)



    geo = Glue(([outer-rotor - sum(magnets), rotor- inner] + magnets))

    return geo


# -----------------------------------------------------------------------------
# ---------------------Parameters
# -----------------------------------------------------------------------------

mu_air = 4e-7*np.pi
mu_magnet = 1.05*mu_air
mu_rotor = mu_air*5e2

sigma_magnet = 8e5
sigma_rotor =  0

order0 = 3

savetime = 1

# -----------------------------------------------------------------------------

# parser = argparse.ArgumentParser()
# parser.add_argument("--nu", type=int, default=3, help="Value of nu")

# parser.add_argument("--PZ", type=int, default=12, help="Number of pole pairs")

# parser.add_argument("--tau", type=int, default=1, help="Value of tau")
# args = parser.parse_args()

# nu = args.nu
# PZ = args.PZ
# tau = args.tau
# -----------------------------------------------------------------------------



nu = 1
PZ = 8
tau = 1

K0= 10000
f = 1e3
omega = 2*np.pi*f

delta = lambda omega, sigma, mu : sqrt(2/(nu*omega*sigma*mu))


delta_rot = delta(omega, 1.86e6, mu_rotor)
delta_mag = delta(omega, sigma_magnet, mu_magnet)
print("Delta_mag ist ", delta_mag)

print("delta_rot", delta_rot)
print(f"Frequenz = {f} und u = {nu}\n")
r_Fe = 144.5e-3
maxh = 0.5
maxh_rotor = max(maxh*delta_rot, 1e-4)
maxh_mag = 4*maxh_rotor
mp = MeshingParameters(maxh = maxh_mag)
mesh = Mesh(OCCGeometry(MakeGeometry(H_L=8e-3, H_M=6e-3, maxh_rotor=maxh_rotor, maxh_mag=maxh_mag, r_Fe=144.5e-3, tau=tau, PZ=PZ, maxh = maxh), dim = 2).GenerateMesh())
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

#Region CF
regions = mesh.MaterialCF({"air":1, "rotor": 2, "magnet": 3})
Draw(regions, mesh, "regions")

# Phi berechnung
def Phi(x,y):
    return atan2(y,x)+np.pi/2

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
f += K(x,y, nu)*test.Trace()*ds("outer") #*cos(2/D*x) PROBLEM weil hier periodische RBs


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
try:
    A = Integrate(1, mesh, definedon=mesh.Materials("magnets.*"))/PZ
    print("Fläche = ", A)
    losses = Integrate(p, mesh, definedon=mesh.Materials("magnets.*"))
    print("P(u, u) = ", losses/PZ)
    print("P/omega = ", losses/omega)


    with open("simulations.txt", "a") as file:
        file.write(f"Periodic Sim: PZ = {PZ}, Tau = {tau}, omega = {omega}, A_magnet = {A*PZ}: losses = {PZ*losses.real}, delta_magnet = {delta_mag}\n")
        print("written to file")
    print("delta_r = ", delta_rot)
    print("delta_m = ", delta_mag)
except Exception as e:
    print("An error occurred: ", e)