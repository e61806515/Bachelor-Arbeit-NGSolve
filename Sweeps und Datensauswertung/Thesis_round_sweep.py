from netgen.occ import *
from ngsolve import *
from contextlib import nullcontext
import time

#from netgen.webgui import Draw as DrawGeo
#from ngsolve.webgui import Draw
from netgen.meshing import MeshingParameters, IdentificationType
import numpy as np
import math as math
import netgen.gui
import matplotlib.pyplot as plot
import argparse
import sys

#ngsglobals.msg_level = 5



#              Geometrie
#
#
#

def MakeGeometry(H_L, H_M, delta_rot, delta_mag, r_Fe, tau, PZ, maxh):

    print("tau = ", tau)
    d_Fe = 5*delta_rot
    r_i = r_Fe - d_Fe
    r_L = r_Fe + H_L
    rect = WorkPlane(Axes((-r_L,-r_L,0), n=Z, h=X)).Rectangle(2*r_L, r_L).Face()

    inner = WorkPlane().Circle(r_i).Face() #- rect
    inner.edges.name="inner"

    rotor = WorkPlane().Circle(r_Fe).Face() #-rect
    rotor.maxh = maxh*delta_rot
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
                magnets[i].maxh = maxh*delta_mag
                #magnets[i].vertices.maxh = 0.00001 #geht nicht??
                magnets[i].name = f"magnets_{i}" #{i:"Anweisungen zur Formatierung"} Formated String
                magnets[i].col = (0, 0, 1)


    else:
        magnets = WorkPlane().Circle(r_Fe+H_M).Face() #-rect
        magnets.name = "magnets"
        magnets.maxh = maxh*delta_mag
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

mu_air = 4e-7*np.pi
mu_magnet = 1.05*mu_air
mu_rotor = mu_air*5e2

sigma_magnet = 1.86e6
sigma_rotor =  1e-12

order0 = 2

savetime = 0

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
tau = 0.999999
Br = 1
K0= 10000
f = 100
omega = 2*np.pi*f

delta = lambda omega, sigma, mu : sqrt(2/(omega*sigma*mu))

#maxh ist der Koeffizient der max Elementgröße im Rotor. d.h. maxh_rotor = maxh*delta_rot
maxh = 0.3
mp = MeshingParameters(maxh=0.4)
mesh = Mesh(OCCGeometry(MakeGeometry(H_L=8e-3, H_M=6e-3, delta_rot = 5e-3, delta_mag = 1e-3, r_Fe=28e-3, tau=tau, PZ=PZ, maxh = maxh), dim = 2).GenerateMesh(mp=mp))
mesh.Curve(3)


#print(mesh.GetMaterials())
#print(set(mesh.GetBoundaries()))
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



sigma = {"air":0, "rotor": sigma_rotor, "magnets.*": sigma_magnet}

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




# Phi berechnung
def Phi(x,y):
    return atan2(y,x)

def K(x,y, nu=1):
     return K0*exp(-1j*nu*Phi(x,y))


x_val = np.logspace(0, 6, 80)
p_values=[]
p_flat =np.loadtxt('p_flat.csv', delimiter=',')
print(len(p_flat))
with (open(f'sweep_round_time_{tau}_{nu}.csv', 'w') if tau is 1 and savetime is 1 else nullcontext()) as time_file:
    with open(f'sweep_round_onlymag_{tau}_{nu}_PZ{PZ}.csv', 'w') as file:
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
                mesh = Mesh(OCCGeometry(MakeGeometry(H_L=8e-3, H_M=6e-3, delta_rot = 5e-3, delta_mag = delta_mag, r_Fe=28e-3, tau=tau, PZ=PZ, maxh = maxh), dim = 2).GenerateMesh(mp=mp))
                mesh.Curve(3)
            muCF = mesh.MaterialCF(mu, default=mu_air)
            sigmaCF = mesh.MaterialCF(sigma, default=0)

            V = H1(mesh, order = 2, dirichlet = "inner", complex=True)
            trial = V.TrialFunction()
            test = V.TestFunction()

            u = GridFunction(V)

            a = BilinearForm(V, symmetric = True)
            a +=  1/muCF*grad(trial)*grad(test) * dx
            a += 1j*omega*sigmaCF*test * trial * dx#("rotor|magnet|air") 1j*

            c = Preconditioner(a, type="direct", inverse = "sparsecholesky")

            f = LinearForm(V)
            f += K(x,y, nu)*test.Trace()*ds("outer")
            start_time = time.time()
            with TaskManager():
                    a.Assemble()
                    f.Assemble()
                    bvp = BVP(bf=a, lf=f, gf=u, pre=c)
                    bvp.Do()
            end_time = time.time()
            elapsed_time = end_time - start_time
            E = -1j * omega * u
            J = sigmaCF * E
            p = E*Conj(J)/2
            losses = Integrate(p, mesh, definedon=mesh.Materials("magnets.*"))
            file.write(f'{freq},{losses.real}\n')
            if(savetime is 1):
                time_file.write(f'{mesh.ne},{elapsed_time}\n')
            p_values.append(losses.real)

print(len(p_values), "vs ", len(x_val))
#write x_val and p_val to one .txt file
#np.savetxt(f'sweep_round_onlymag_{tau}_{nu}.csv', np.column_stack((x_val, p_values)), delimiter=',', fmt='%s')
exit()

if(False):
    plot.plot(x_val, p_values)
    plot.title('Eddy Current 1st O. Losses Half-Round, tau=3/4, PZ=2')
    plot.xscale('log')
    plot.yscale('log')
    plot.xlabel('Frequency (Hz)')
    plot.ylabel('Power Losses (W/m)')
    plot.grid(True, which='both')
    plot.savefig('Losses_Deg_PZ_10_3Phase.pdf', i=PZ, format='pdf')
    plot.show()


    r_vs_f=[]
    for i in np.arange(0,len(x_val)):
        d = p_values[i]/p_flat[i]
        r_vs_f.append(d)

    plot.semilogx(x_val, r_vs_f)
    plot.title('Eddy Current 1st Order Losses Ratio Flat vs Half-Round, tau=3/4, PZ=2')
    plot.xlabel('Frequency (Hz)')
    plot.ylabel('Ratio Power Losses (W/m)')
    plot.grid(True, which='both')
    plot.savefig('Ratio_1st_Order_180_Deg_3Phase.pdf', format='pdf')
    plot.show()

print(mesh.GetMaterials())
print("A = ", Integrate(1, mesh, region_wise=True))
