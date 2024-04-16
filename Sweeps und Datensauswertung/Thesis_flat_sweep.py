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
import matplotlib.pyplot as plot

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

def MakeGeometry(d_M, delta_rot, delta_mag, d_L, tau_p, tau, maxh):
        d_Fe = 5*delta_rot
        b = tau_p
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
                lines = [ (0, 1, 0, 1, "outer"), (7, 6, 3, 1, "inner_right"), (5, 4, 2, 1, "magnet"),
                 (3, 6, 2, 3, "magnet"), (3, 2, 3, 1, "inner_left"), (8, 9, 3, 0, "inner")]
                for p1, p2, left, right, bc in lines:
                        geo.Append( ["line", pnums[p1], pnums[p2]], leftdomain=left, rightdomain=right, bc = bc)
                rotor_l = geo.Append(["line", pnums[2], pnums[8]], leftdomain=3, rightdomain=0, bc="rotor_left")
                air_l = geo.Append(["line", pnums[2], pnums[0]], leftdomain=0, rightdomain=1, bc="air_left")
                geo.Append(["line", pnums[7], pnums[9]], leftdomain=0, rightdomain=3, bc="rotor_right", copy= rotor_l)
                geo.Append(["line", pnums[7], pnums[1]], leftdomain=1, rightdomain=0, bc="air_right", copy= air_l)
                geo.Append(["line", pnums[4], pnums[3]], leftdomain=2, rightdomain=1, bc="magnet_left", maxh = 0.0005)
                geo.Append(["line", pnums[6], pnums[5]], leftdomain=2, rightdomain=1, bc="magnet_right", maxh = 0.0005)

        else:
                lines = [ (0, 1, 0, 1, "outer"), (5, 4, 2, 1, "magnet"),
                 (3, 6, 2, 3, "magnet"), (8, 9, 3, 0, "inner")]
                for p1, p2, left, right, bc in lines:
                        geo.Append( ["line", pnums[p1], pnums[p2]], leftdomain=left, rightdomain=right, bc = bc, maxh=0.001)
                rotor_l = geo.Append(["line", pnums[3], pnums[8]], leftdomain=3, rightdomain=0, bc="rotor_left")
                air_l = geo.Append(["line", pnums[0], pnums[4]], leftdomain=1, rightdomain=0, bc="air_left", maxh=0.0005)
                magnet_l = geo.Append(["line", pnums[4], pnums[3]], leftdomain=2, rightdomain=0, bc="magnet_left", maxh = 0.0005)
                geo.Append(["line", pnums[6], pnums[9]], leftdomain=0, rightdomain=3, bc="rotor_right", copy=rotor_l)
                geo.Append(["line", pnums[1], pnums[5]], leftdomain=0, rightdomain=1, bc="air_right", copy=air_l, maxh = 0.0005)
                geo.Append(["line", pnums[5], pnums[6]], leftdomain=0, rightdomain=2, bc="magnet_right", copy = magnet_l, maxh = 0.0005)


        geo.SetMaterial( 1, "air")
        geo.SetMaterial( 2, "magnet")
        geo.SetMaterial( 3, "rotor")

        geo.SetDomainMaxH(1, 0.001)
        geo.SetDomainMaxH(2, maxh*delta_mag)
        geo.SetDomainMaxH(3, maxh*delta_rot)
        return geo

#b = np.pi * 38.19e-3      #Breite b entspricht halben Umfang, also
                #Pi*r mit r = r_rot + h_magnet/2, halber Umfang ist pole pitch bei Schmid
                #liefert
PZ = 2
#A = 0.00010798
tau = 1
#tau = A/((6e-3)*tau) #liefert 0.1199773 für PZ 2 und 0.02399555 für PZ 10
#print("b lautet ", b)
d_L = 2e-3
d_M = 3*d_L
tau_p = 60*d_L

nu =7

f = 50
omega = 2*np.pi*f
K0 = 10000
mu_air = 4e-7*np.pi
mu_magnet = 1.05*mu_air
mu_rotor = mu_air*5e2
mu = {"air":mu_air, "rotor": mu_rotor, "magnet": mu_magnet}

sigma_rotor =1e-12
sigma_magnet = 1.86e6
sigma = {"air":0, "rotor": sigma_rotor, "magnet": sigma_magnet}

sigma_v = {"air":None, "rotor": None, "magnet": sigma_magnet}

delta = lambda omega, sigma, mu : sqrt(2/(omega*sigma*mu))

#maxh ist der Koeffizient der max Elementgröße im Rotor. d.h. maxh_rotor = maxh*delta_rot
maxh = 0.3

geo = MakeGeometry(d_M=d_M, delta_rot = 5e-3, delta_mag=1e-3, d_L=2e-3, tau_p=tau_p, tau=tau, maxh = maxh)
print(geo.GetNSplines())
mesh = geo.GenerateMesh(maxh = 0.4)
mesh = Mesh(mesh)

def Phi(x):
    return 2*pi/(PZ*tau_p)

def K(x, nu=1):
      K = K0*exp(-1j*np.pi*x*nu/(tau_p))
      return K

#Randwert Probleme u. (Bi-)Linearform

#V = Periodic(H1(mesh, order = 2, dirichlet = "inner", complex=True), phase=[0,0])
if tau <1:
        phase = [-1,-1]
else:
        phase = [-1,-1,-1]
x_val = np.logspace(0, 6, 80)

p_values=[]
print(f'sweep_flat_onlymag_{tau}_{nu}.csv')
with open(f'sweep_flat_onlymag_{tau}_{nu}_PZ{PZ}.csv', 'w') as file:
        for freq in x_val:
                omega = 2*np.pi*freq
                if(sigma_rotor>1e-12):
                        delta_rot = delta(omega, sigma_rotor, mu_rotor)
                else:
                        delta_rot = 5e-3

                delta_mag = delta(omega, sigma_magnet, mu_magnet)
                print("delta_mag = ", delta_mag)
                print("f = ", freq)

                if delta_mag<1e-3:
                        geo = MakeGeometry(d_M=d_M, delta_rot = delta_rot, delta_mag = delta_mag, d_L=2e-3, tau_p=tau_p, tau=tau, maxh = maxh)
                        print(geo.GetNSplines())
                        mesh = geo.GenerateMesh(maxh = 0.4)
                        mesh = Mesh(mesh)

                muCF = mesh.MaterialCF(mu, default=mu_air)

                sigmaCF = mesh.MaterialCF(sigma, default=0)

                sigma_visual = mesh.MaterialCF(sigma_v)

                V = Periodic(H1(mesh, order = 2, dirichlet = "inner", complex=True), phase = phase)
                trial = V.TrialFunction()
                test = V.TestFunction()


                u = GridFunction(V)

                a = BilinearForm(V, symmetric = True)
                a +=  1/muCF*grad(trial)*grad(test) * dx
                a += 1j*omega*sigmaCF*test * trial * dx#("rotor|magnet|air") 1j*

                c = Preconditioner(a, type="direct", inverse = "sparsecholesky")

                f = LinearForm(V)
                f += K(x, nu=nu)*test.Trace()*ds("outer")

                with TaskManager():
                        a.Assemble()
                        f.Assemble()
                        bvp = BVP(bf=a, lf=f, gf=u, pre=c)
                        bvp.Do()
                E = -1j * omega * u
                J = sigmaCF * E
                p = E*Conj(J)/2
                losses = Integrate(p, mesh, definedon=mesh.Materials("magnet"))*PZ
                file.write(f'{freq},{losses.real}\n')
                p_values.append(losses.real)
#np.savetxt(f'sweep_flat_both_{nu}.csv', np.column_stack((x_val, p_values)), delimiter=',', fmt='%s')
exit()
if(False):
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





