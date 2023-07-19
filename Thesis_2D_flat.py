from ngsolve import *
import netgen.geom2d as ng
from netgen.occ import *
from ngsolve import *
from netgen.csg import *

#from netgen.webgui import Draw as DrawGeo
#from ngsolve.webgui import Draw
from netgen.meshing import MeshingParameters, IdentificationType
import numpy as np
import math as math
import netgen.gui
ngsglobals.msg_level = 5
#
#
#                       (outer)
#             0------------>--------------1
#             |           (air)           |
#             |      4------<------5      |
#      (left) |      |             |      |(right)
#             |      |   (magnet)  |      |
#             |      |             |      |
#             2--<---3------>------6--<---7 
#             | (i.l.)             (i.r.) |
#             |                           |(right)
#     (left)  |          (rotor)          |
#             |                           |
#             |                           |
#             |                           |
#             |                           |
#             8-------------->------------9
#                       (inner)

def MakeGeometry(d_M, d_Fe, d_L, b, tau_p):
        b_M = b*tau_p
        
        #Add points to geometry
        #geo = WorkPlane((-b_M/2, 0, 0), n=Z, h=X).Rectangle(b_M, d_M).Face()
        #  
        geo = ng.SplineGeometry()

        pnts = [(0,(d_M+d_L)), (b,(d_M+d_L)), (0,0), 
                ((b-b_M)/2, 0), ((b-b_M)/2, d_M), ((b+b_M)/2, d_M), 
                ((b+b_M)/2, 0), (b,0), (0, -d_Fe), (b, -d_Fe)]
        pnums = [geo.AppendPoint(*p) for p in pnts]

        #Add lines and domains
        # start-point, end-point, (boundary-condition,) domain on left side, domain on right side:
        lines = [ (0, 1, 0, 1, "outer"), (7, 6, 3, 1, "inner_right"), (5, 4, 2, 1, "magnet"), 
                 (3, 6, 2, 3, "magnet"), (3, 2, 3, 1, "inner_left"), (8, 9, 3, 0, "inner")]
        if tau_p < 1:
                lines.append((4, 3, 2, 1, "magnet_left"))
                lines.append((6, 5, 2, 1, "magnet_right"))

                for p1, p2, left, right, bc in lines:
                        geo.Append( ["line", pnums[p1], pnums[p2]], leftdomain=left, rightdomain=right, bc = bc, maxh=0.001)
                rotor_l = geo.Append(["line", pnums[2], pnums[8]], leftdomain=3, rightdomain=0, bc="rotor_left", maxh = 0.001)
                air_l = geo.Append(["line", pnums[2], pnums[0]], leftdomain=0, rightdomain=1, bc="air_left", maxh = 0.001)
                geo.Append(["line", pnums[9], pnums[7]], leftdomain=3, rightdomain=0, bc="rotor_right", maxh=0.001) #, copy=rotor_l)
                geo.Append(["line", pnums[1], pnums[7]], leftdomain=0, rightdomain=1, bc="air_right", maxh = 0.001) #, copy=air_l)
        else:
                for p1, p2, left, right, bc in lines:
                        geo.Append( ["line", pnums[p1], pnums[p2]], leftdomain=left, rightdomain=right, bc = bc, maxh=0.001)
                rotor_l = geo.Append(["line", pnums[2], pnums[8]], leftdomain=3, rightdomain=0, bc="rotor_left")
                air_l = geo.Append(["line", pnums[2], pnums[0]], leftdomain=0, rightdomain=1, bc="air_left")
                magnet_l = geo.Append(["line", pnums[4], pnums[3]], leftdomain=2, rightdomain=0, bc="magnet_left")
                geo.Append(["line", pnums[9], pnums[7]], leftdomain=3, rightdomain=0, bc="rotor_right", copy=rotor_l)
                geo.Append(["line", pnums[1], pnums[7]], leftdomain=0, rightdomain=1, bc="air_right", copy=air_l)
                geo.Append(["line", pnums[6], pnums[5]], leftdomain=2, rightdomain=0, bc="magnet_right", copy = magnet_l)

        
        geo.SetMaterial( 1, "air")
        geo.SetMaterial( 2, "magnet")
        geo.SetMaterial( 3, "rotor")

        return geo

b = np.pi/2 * 76.38e-3
mesh = MakeGeometry(6e-3, 26.19e-3, 2e-3, b, 3/4)
Draw(mesh)

print(mesh.GetNSplines())

mesh = mesh.GenerateMesh(maxh = 0.01)
mesh = Mesh(mesh)



#Parameters and Coefficient Functions
omega = 50
K0 = 10000

mu_air = 4e-7*np.pi
mu_magnet = 1.08*mu_air
mu_rotor = mu_air*5e3
mu = {"air":mu_air, "rotor": mu_rotor, "magnet": mu_magnet}
#version_2 = CF([val[mat] if mat in val.keys() else 0 for mat in mesh.GetMaterials()])
muCF = mesh.MaterialCF(mu, default=mu_air)

sigma = {"air":0, "rotor": 1.86e6, "magnet": 0.667e6}
sigmaCF = mesh.MaterialCF(sigma, default=0)
sigma_v = {"air":None, "rotor": None, "magnet": 0.667e6}
sigma_visual = mesh.MaterialCF(sigma_v)

versions = [muCF, sigmaCF]
[Draw(versions[i], mesh, str(i)) for i in range(len(versions))]
"""sigmaCF = CF([0, 10e3, 10e6]) #[0, 10e3, 10e6]
sigma_visual = CF([0, 10e3, 0])
mu_air = 4e-7*np.pi
mu_magnet = 10*mu_air
mu_rotor = mu_air*1000
mu = CF([mu_air, mu_magnet, mu_rotor])

Kvalues = {mat: 0 for mat in mesh.GetBoundaries()}
Kvalues["outer"] = K0
K0CF = CF([Kvalues[mat] for mat in mesh.GetBoundaries()])"""
def Phi(x):
    return asin(2*x/b-1)

def K(x):
     K1 = -K0*sin(Phi(x))
     K2 = -K0*sin(Phi(x) + 2*pi/3)*exp(1j*2*pi/3*omega)
     K3 = -K0*sin(Phi(x) + 4*pi/3)*exp(1j*4*pi/3*omega)
     return K1 + K2 + K3     

#Randwert Probleme u. (Bi-)Linearform
coef_dirichlet = CF([0,0,0,0,0,0,0,0,0,0,0,0])

#V = Periodic(H1(mesh, order = 2, dirichlet = "inner", complex=True), phase=[0,0])
V = H1(mesh, order = 2, dirichlet = "inner", complex=True)
trial = V.TrialFunction()
test = V.TestFunction()
u = GridFunction(V)

a = BilinearForm(V, symmetric = True)
a +=  1/muCF*grad(trial)*grad(test) * dx
a += 1j*omega*sigmaCF*test * trial * dx#("rotor|magnet|air") 1j*

c = Preconditioner(a, type="direct", inverse = "sparsecholesky")

f = LinearForm(V)
f += K(x)*test.Trace()*ds("outer")   #*cos(2/D*x)  # set neumann .. macht das Sinn? ...
                                                    #Wir wollen am oberen Neumann-Rand 
# .Trace() reduziert Domain-funktionen auf die Grenzen. Bei "ds"-Integralen sinnvoll. 
#"dx" für VOL-Integrale, "ds" für Surface-Integrale
#u.Set(coef_dirichlet, BND) #Wird ohne bvp-solver überschrieben wenn inhomogene BC

#u.Set(bndCF, BND) # set the drc bnds
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
#Draw(CF([1,2,3]), mesh, "materials")
Draw(K(x), mesh, 'K')
print(u.vec.Norm())

p = sigma_visual*omega*omega*u*Conj(u)/2
#energy = Integrate(p, mesh, region_wise= True)
energy = Integrate(p, mesh)

print("P(u, u) = ", energy*2)
print("P/omega = ", energy*2/omega)

delta_rot = sqrt(2/(omega*1.86e6*mu_rotor))
delta_mag = sqrt(2/(omega*0.667e6*mu_magnet))
print("delta_r = ", delta_rot)
print("delta_m = ", delta_mag)

#input()



