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
#ngsglobals.msg_level = 5
#
#
#                       (outer)
#             0------------>--------------1
#             |           (air)  1        |
#             |      4------<------5      |
#      (left) |      |             |      |(right)
#             |      |   (magnet)  |      |
#             |      |       2     |      |
#     0       2--<---3------>------6--<---7
#             | (i.l.)             (i.r.) |
#             |                           |(right)
#     (left)  |          (rotor)          |
#             |             3             |
#             |                           |
#             |                           |
#             |                           |
#             8-------------->------------9
#                       (inner)

def MakeGeometry(d_M, delta_rot, delta_mag, d_L, b, tau, maxh, faktor_d_rotor, maxh_mp):
        d_Fe = faktor_d_rotor*delta_rot
        b_M = b*tau

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

        if tau < 1:
                lines = [ (0, 1, 0, 1, "outer", maxh_mp), (7, 6, 3, 1, "inner_right", 2*maxh*delta_rot), (5, 4, 2, 1, "magnet", maxh*delta_mag),
                 (3, 6, 2, 3, "magnet", 2*maxh*delta_rot), (3, 2, 3, 1, "inner_left", 2*maxh*delta_rot), (8, 9, 3, 0, "inner", 2*maxh*delta_rot)]
                for p1, p2, left, right, bc, maxh_ in lines:
                        geo.Append( ["line", pnums[p1], pnums[p2]], leftdomain=left, rightdomain=right, bc = bc, maxh=maxh_)
                rotor_l = geo.Append(["line", pnums[2], pnums[8]], leftdomain=3, rightdomain=0, bc="rotor_left", maxh=2*maxh*delta_rot)
                air_l = geo.Append(["line", pnums[2], pnums[0]], leftdomain=0, rightdomain=1, bc="air_left", maxh=maxh_mp)
                geo.Append(["line", pnums[7], pnums[9]], leftdomain=0, rightdomain=3, bc="rotor_right", copy= rotor_l, maxh=2*maxh*delta_rot)
                geo.Append(["line", pnums[7], pnums[1]], leftdomain=1, rightdomain=0, bc="air_right", copy= air_l, maxh=maxh_mp)
                geo.Append(["line", pnums[4], pnums[3]], leftdomain=2, rightdomain=1, bc="magnet_left", maxh = maxh * delta_mag)
                geo.Append(["line", pnums[6], pnums[5]], leftdomain=2, rightdomain=1, bc="magnet_right", maxh = maxh * delta_mag)

        else:
                lines = [ (0, 1, 0, 1, "outer", maxh_mp), (5, 4, 2, 1, "magnet", maxh*delta_mag),
                 (3, 6, 2, 3, "magnet", 2*maxh*delta_rot), (8, 9, 3, 0, "inner", 2*maxh*delta_rot)]
                for p1, p2, left, right, bc, maxh_ in lines:
                        geo.Append( ["line", pnums[p1], pnums[p2]], leftdomain=left, rightdomain=right, bc = bc, maxh=maxh_)
                rotor_l = geo.Append(["line", pnums[3], pnums[8]], leftdomain=3, rightdomain=0, bc="rotor_left", maxh=2*maxh*delta_rot)
                air_l = geo.Append(["line", pnums[0], pnums[4]], leftdomain=1, rightdomain=0, bc="air_left", maxh=maxh_mp)
                magnet_l = geo.Append(["line", pnums[4], pnums[3]], leftdomain=2, rightdomain=0, bc="magnet_left", maxh = maxh * delta_mag)
                geo.Append(["line", pnums[6], pnums[9]], leftdomain=0, rightdomain=3, bc="rotor_right", copy=rotor_l, maxh=2*maxh*delta_rot)
                geo.Append(["line", pnums[1], pnums[5]], leftdomain=0, rightdomain=1, bc="air_right", copy=air_l, maxh=maxh_mp)
                geo.Append(["line", pnums[5], pnums[6]], leftdomain=0, rightdomain=2, bc="magnet_right", copy = magnet_l, maxh = maxh * delta_mag)


        geo.SetMaterial( 1, "air")
        geo.SetMaterial( 2, "magnet")
        geo.SetMaterial( 3, "rotor")
        """
        geo.SetDomainMaxH(1, 0.003)
        geo.SetDomainMaxH(2, maxh*delta_mag)
        geo.SetDomainMaxH(3, maxh*delta_rot)
        """
        return geo

#b = np.pi * 38.19e-3      #Breite b entspricht halben Umfang, also
                #Pi*r mit r = r_rot + h_magnet/2, halber Umfang ist pole pitch bei Schmid
                #liefert

d_L = 2e-3
d_M = 3*d_L
b = 60*d_L

nu = 9
PZ = 16
order0=3
tau=1
f_dr = 8
A_mags=b*tau*d_M
f = 1e6
omega = 2*np.pi*f
K0 = 10000
mu_air = 4e-7*np.pi
mu_magnet = 1.05*mu_air
mu_rotor = mu_air*5e2
mu = {"air":mu_air, "rotor": mu_rotor, "magnet": mu_magnet}

sigma_rotor = 0
sigma_magnet = 8e5
sigma = {"air":0, "rotor": sigma_rotor, "magnet": sigma_magnet}

sigma_v = {"air":None, "rotor": None, "magnet": sigma_magnet}

delta = lambda omega, sigma, mu : sqrt(2/(nu*omega*sigma*mu))

delta_rot = delta(omega, 1.86e6, mu_rotor)

delta_mag = delta(omega, sigma["magnet"], mu_magnet)
maxh = sqrt(nu)*0.7
geo = MakeGeometry(d_M=d_M, delta_rot = delta_rot, delta_mag=delta_mag, d_L=2e-3, b=b, tau=tau, maxh = maxh, faktor_d_rotor=f_dr, maxh_mp=0.1)
print(geo.GetNSplines())
mesh = geo.GenerateMesh(maxh = 0.1)
mesh = Mesh(mesh)

muCF = mesh.MaterialCF(mu, default=mu_air)

sigmaCF = mesh.MaterialCF(sigma, default=0)

sigma_visual = mesh.MaterialCF(sigma_v)

versions = [muCF, sigmaCF]
[Draw(versions[i], mesh, str(i)) for i in range(len(versions))]
Draw(sigma_visual, mesh, "sigma_v")

#Region CF
regions = mesh.MaterialCF({"air":1, "rotor": 2, "magnet": 3})
Draw(regions, mesh, "regions")

def Phi(x):
    return 2*pi*x/(PZ*b)

def K(x, nu=1):
      K = K0*exp(-1j*np.pi*x*nu/(b))
      return K
"""
def K(x, nu=1):
      return K0*exp(-1j*nu*PZ*Phi(x)/2)
"""
#Randwert Probleme u. (Bi-)Linearform

#V = Periodic(H1(mesh, order = 2, dirichlet = "inner", complex=True), phase=[0,0])
if tau <1:
        phase = [-1,-1]
else:
        phase = [-1,-1,-1]
V = Periodic(H1(mesh, order = order0, dirichlet = "inner", complex=True), phase = phase)
trial = V.TrialFunction()
test = V.TestFunction()
u = GridFunction(V)

a = BilinearForm(V, symmetric = True)
a +=  1/muCF*grad(trial)*grad(test) * dx
a += 1j*omega*sigmaCF*test * trial * dx#("rotor|magnet|air") 1j*

#c = Preconditioner(a, type="direct", inverse = "sparsecholesky")
c = Preconditioner(a, type="direct", inverse = "sparsecholesky")
f = LinearForm(V)
f += K(x, nu=nu)*test*ds("outer")   #*cos(2/D*x)  # set neumann .. macht das Sinn? ...
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

p = sigma_visual*omega*omega*u*Conj(u)/2
#losses = Integrate(p, mesh, region_wise= True)
# A = Integrate(1, mesh, definedon=mesh.Materials("magnet"))
# losses = Integrate(p, mesh, definedon=mesh.Materials("magnet"))
# print("Fläche = ", A)
# print("P(u, u) = ", PZ*losses)
# print("P/omega = ", PZ*  losses/omega)

# print("delta_r = ", delta_rot)
# print("delta_m = ", delta_mag)

# with open("simulations.txt", "a") as file:
#     file.write(f"Flat Sim: PZ = {PZ}, Tau = {tau}, f = {f}, A_magnet = {A}: losses = {losses}, delta_magnet = {delta_mag} \n")
#     print("written to file")
try:
    A = Integrate(1, mesh, definedon=mesh.Materials("magnet"))
    losses = Integrate(p, mesh, definedon=mesh.Materials("magnet"))
    print(f"Frequenz = {f}")
    print("Fläche = ", A*PZ)
    print("P(u, u) = ", PZ*losses)
    print("P/omega = ", PZ*  losses/omega)
    print("losses/Area = ", losses/A)
    print("delta_r = ", delta_rot)
    print("delta_m = ", delta_mag)
except Exception as e:
    print("An error occurred: ", e)

try:
    with open("simulations.txt", "a") as file:
        file.write(f"Flat Sim A_mag: PZ = {PZ}, Tau = {tau}, omega = {omega}, A_magnet = {A*PZ}: losses = {losses.real}, delta_magnet = {delta_mag} \n")
        print("written to file")
except Exception as e:
    print("An error occurred while writing to file: ", e)
#input()



