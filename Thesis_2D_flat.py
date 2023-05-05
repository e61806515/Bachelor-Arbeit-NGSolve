from ngsolve import *
import netgen.geom2d as ng
import numpy as np
import math as math
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

def MakeGeometry(d_M, d_Fe, d_L, D, tau_p):
        b = np.pi/2 * D
        b_M = b*tau_p
        
        #Add points to geometry
        geo = ng.SplineGeometry()
        pnts = [(0,(d_M+d_L)), (b,(d_M+d_L)), (0,0), 
                ((b-b_M)/2, 0), ((b-b_M)/2, d_M), ((b+b_M)/2, d_M), 
                ((b+b_M)/2, 0), (b,0), (0, -d_Fe), (b, -d_Fe)]
        pnums = [geo.AppendPoint(*p) for p in pnts]

        #Add lines and domains
        # start-point, end-point, (boundary-condition,) domain on left side, domain on right side:
        lines = [ (0, 1, 0, 1, "outer"), (1, 7, 0, 1, "air_right"), (7, 6, 3, 1, "inner_right"), 
                 (6, 5, 2, 1, "magnet"), (5, 4, 2, 1, "magnet"), (4, 3, 2, 1, "magnet"), 
                 (3, 6, 2, 3, "magnet"), (3, 2, 3, 1, "inner_left"), (2, 0, 0, 1, "air_left"), 
                 (2, 8, 3, 0, "rotor_left"), (8, 9, 3, 0, "inner"), (9, 7, 3, 0, "rotor_right")]

        for p1, p2, left, right, bc in lines:
                geo.Append( ["line", pnums[p1], pnums[p2]], leftdomain=left, rightdomain=right, bc = bc, maxh=0.001)
        geo.SetMaterial( 1, "air")
        geo.SetMaterial( 2, "magnet")
        geo.SetMaterial( 3, "rotor")

        #geo.SetDomainMaxH(1, 1000)
        return geo

D = 100e-3
mesh = MakeGeometry(6e-3, 20e-3, 1e-3, D, 5/6).GenerateMesh(maxh = 0.01)
mesh = Mesh(mesh)


#Parameters and Coefficient Functions
omega = 50
#sigma = {"air" : 0, "magnet" : 10e6, "rotor" : 10e3}
#sigmaCF = CF([sigma[mat] for mat in mesh.GetMaterials()])

sigmaCF = CF([0, 0, 0]) #[0, 10e3, 10e6]
sigma_visual = CF([0, 10e3, 0])
mu_air = 4e-7*np.pi
mu_magnet = 10*mu_air
mu_rotor = mu_air*1000
mu = CF([mu_air, mu_magnet, mu_rotor])
K0= 10
Kvalues = {mat: 0 for mat in mesh.GetBoundaries()}
Kvalues["outer"] = K0
K0CF = CF([Kvalues[mat] for mat in mesh.GetBoundaries()])


print(mesh.GetMaterials())
#print(mesh.GetBoundaries())

#Randwert Probleme u. (Bi-)Linearform
coef_dirichlet = CF([0,0,0,0,0,0,0,0,0,0,0,0])
#coef_dirichlet = mesh.BoundaryCF(default = 0)
#bnd_dir = {"air_left" : 0, "air_right" : 0, "inner" : 0}
#bndCF = CF([bnd_dir[bnd] for bnd in mesh.GetBoundaries()])

V = H1(mesh, order = 2, dirichlet = "inner", complex=True)

#print(V)

trial = V.TrialFunction()
test = V.TestFunction()
u = GridFunction(V)

a = BilinearForm(V, symmetric = True)
a +=  1/mu*grad(trial)*grad(test) * dx
a += 1j*omega*sigmaCF*test * trial * dx#("rotor|magnet|air") 1j*

c = Preconditioner(a, type="direct", inverse = "sparsecholesky")

f = LinearForm(V)
f += -K0*cos(2/D*x)*test.Trace()*ds("outer")   #*cos(2/D*x)  # set neumann .. macht das Sinn? ...
                                                        #Wir wollen am oberen Neumann-Rand 
# .Trace() reduziert Domain-funktionen auf die Grenzen. Bei "ds"-Integralen sinnvoll. 
#"dx" für VOL-Integrale, "ds" für Surface-Integrale
u.Set(coef_dirichlet, BND) #Wird ohne bvp-solver überschrieben wenn inhomogene BC

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
Draw(1/mu*B, mesh, 'H')
Draw(Norm(1/mu*B[0]), mesh, 'Norm Hx')
Draw(Norm(1/mu*B[1]), mesh, 'Norm Hy')
Draw(Norm(B[0]), mesh, 'Norm Bx')
Draw(Norm(B[1]), mesh, 'Norm By')
Draw(-1j*omega*sigma_visual*u, mesh, 'Jz') #Wirbelströme J = -jomega*sigma*A
Draw(CF([1,2,3]), mesh, "materials")
print(u.vec.Norm())

p = sigma_visual*omega*omega*u*Conj(u)/2
energy = Integrate(p, mesh)
#Integral über 1 dOmega muss Fläche liefern

print("P(u, u) = ", energy)

#input()



