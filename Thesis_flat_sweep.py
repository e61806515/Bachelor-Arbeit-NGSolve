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
import matplotlib.pyplot as plot

#ngsglobals.msg_level = 5
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

                for p1, p2, left, right, bc in lines:
                        geo.Append( ["line", pnums[p1], pnums[p2]], leftdomain=left, rightdomain=right, bc = bc, maxh = 0.0005)
                rotor_l = geo.Append(["line", pnums[2], pnums[8]], leftdomain=3, rightdomain=0, bc="rotor_left", maxh = 0.0005)
                air_l = geo.Append(["line", pnums[2], pnums[0]], leftdomain=0, rightdomain=1, bc="air_left", maxh = 0.0005)
                geo.Append(["line", pnums[7], pnums[9]], leftdomain=0, rightdomain=3, bc="rotor_right", copy= rotor_l, maxh = 0.0005)
                geo.Append(["line", pnums[7], pnums[1]], leftdomain=1, rightdomain=0, bc="air_right", copy= air_l, maxh = 0.0005)
                geo.Append(["line", pnums[4], pnums[3]], leftdomain=2, rightdomain=1, bc="magnet_left", maxh = 0.0005)
                geo.Append(["line", pnums[6], pnums[5]], leftdomain=2, rightdomain=1, bc="magnet_right", maxh = 0.0005)

        else:
                for p1, p2, left, right, bc in lines:
                        geo.Append( ["line", pnums[p1], pnums[p2]], leftdomain=left, rightdomain=right, bc = bc, maxh = 0.0005)
                rotor_l = geo.Append(["line", pnums[2], pnums[8]], leftdomain=3, rightdomain=0, bc="rotor_left", maxh = 0.0005)
                air_l = geo.Append(["line", pnums[2], pnums[0]], leftdomain=0, rightdomain=1, bc="air_left", maxh = 0.0005)
                magnet_l = geo.Append(["line", pnums[4], pnums[3]], leftdomain=2, rightdomain=0, bc="magnet_left", maxh = 0.0005)
                geo.Append(["line", pnums[9], pnums[7]], leftdomain=3, rightdomain=0, bc="rotor_right", copy=rotor_l, maxh = 0.0005)
                geo.Append(["line", pnums[1], pnums[7]], leftdomain=0, rightdomain=1, bc="air_right", copy=air_l, maxh = 0.0005)
                geo.Append(["line", pnums[6], pnums[5]], leftdomain=2, rightdomain=0, bc="magnet_right", copy = magnet_l, maxh = 0.0005)

        
        geo.SetMaterial( 1, "air")
        geo.SetMaterial( 2, "magnet")
        geo.SetMaterial( 3, "rotor")

        geo.SetDomainMaxH(2, 0.001)
        
        return geo

#b = np.pi * 38.19e-3      #Breite b entspricht halben Umfang, also 
                #Pi*r mit r = r_rot + h_magnet/2, halber Umfang ist pole pitch bei Schmid
                #liefert 

sweep = True

PZ = 4
A = 0.00010798
tau = 3/4
r_M = 42.19e-3
r = r_M - 1.5e-3
#b = A/((6e-3)*tau)
b = r*2*pi/PZ
geo = MakeGeometry(d_M=6e-3, d_Fe=26.19e-3, d_L=2e-3, b=b, tau_p=tau)
Draw(geo)
print(geo.GetNSplines())
mesh = geo.GenerateMesh(maxh = 0.001)
mesh = Mesh(mesh)



#Parameters and Coefficient Functions
K0 = 10000

mu_air = 4e-7*np.pi
mu_magnet = 1.08*mu_air
mu_rotor = mu_air*5e3
mu = {"air":mu_air, "rotor": mu_rotor, "magnet": mu_magnet}
muCF = mesh.MaterialCF(mu, default=mu_air)

sigma = {"air":0, "rotor": 1.86e6, "magnet": 0.667e6}
sigmaCF = mesh.MaterialCF(sigma, default=0)
sigma_v = {"air":None, "rotor": None, "magnet": 0.667e6}
sigma_visual = mesh.MaterialCF(sigma_v)

versions = [muCF, sigmaCF]

def Phi(x):
    return x*2*pi/(PZ*b)

def K(x):
    K=0 
    
    if(True):
        K1 = -K0*cos(Phi(x))
        K2 = -K0*cos(Phi(x) + 2*pi/3)*exp(1j*2*pi/3)
        K3 = -K0*cos(Phi(x) + 4*pi/3)*exp(1j*4*pi/3)
        K = K + K1 + K2 + K3
    if(True):
        K4 = -K0*cos(Phi(x) + pi/2)
        K5 = -K0*cos(Phi(x) + 2*pi/3 + pi/2)*exp(1j*2*pi/3)
        K6 = -K0*cos(Phi(x) + 4*pi/3 + pi/2)*exp(1j*4*pi/3)
        K = K + K4 + K5 + K6 
    return K/2

#Randwert Probleme u. (Bi-)Linearform
coef_dirichlet = CF([0,0,0,0,0,0,0,0,0,0,0,0])
V = Periodic(H1(mesh, order = 2, dirichlet = "inner", complex=True), phase = [-1,-1])
trial = V.TrialFunction()
test = V.TestFunction()

x_val = np.logspace(0, 6, 100)
p_values=[]
for omega in x_val: 
        if(sweep!=True):
         omega = 50
        u = GridFunction(V)

        a = BilinearForm(V, symmetric = True)
        a +=  1/muCF*grad(trial)*grad(test) * dx
        a += 1j*omega*sigmaCF*test * trial * dx#("rotor|magnet|air") 1j*

        c = Preconditioner(a, type="direct", inverse = "sparsecholesky")

        f = LinearForm(V)
        f += K(x)*test.Trace()*ds("outer")   

        with TaskManager():
                a.Assemble()
                f.Assemble()
                bvp = BVP(bf=a, lf=f, gf=u, pre=c)
                bvp.Do()

        p = sigma_visual*omega*omega*u*Conj(u)/2
        energy = Integrate(p, mesh)*PZ
        if(sweep!=True):
                print("Verluste = ", energy)
                Draw(u*1j*omega*sigmaCF, mesh, 'J')
                break
        p_values.append(energy.real)

if(sweep==True):        
        np.savetxt('p_flat.csv', p_values, delimiter=',')
        plot.plot(x_val, p_values)
        plot.title('Eddy Current 1st O. Losses Flat 180°, Tau=3/4, PZ=2')
        plot.xscale('log')
        plot.yscale('log')
        plot.xlabel('Frequency (Hz)')
        plot.ylabel('Power Losses (W/m)')
        plot.grid(True, which='both')
        plot.savefig('Losses_Flat_180_Deg_3Phase.pdf', format='pdf')
        plot.show()





