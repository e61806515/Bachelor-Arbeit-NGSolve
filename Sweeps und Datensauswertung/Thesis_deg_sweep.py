from netgen.occ import *
from ngsolve import *
from contextlib import nullcontext
import time

#from netgen.webgui import Draw as DrawGeo
#from ngsolve.webgui import Draw
from netgen.meshing import MeshingParameters, IdentificationType
import numpy as np
import math as math
import matplotlib.pyplot as plot


############################################
def drawBndAll(mesh, block=True, old=None, drawFunc=None, useBBnd=False):

    if old == None:
        old = True if mesh.dim == 2 else False
    [drawBnd(mesh, x, block, old, drawFunc=drawFunc, useBBnd=useBBnd) for x in set(mesh.GetBoundaries() if not useBBnd else mesh.GetBBoundaries())]

def drawBnd(mesh, name="bottom|right|top|left|ibot|itop|interface|ileft|iright|natural", block=False, old=False, drawFunc=None, useBBnd=False):
    from ngsolve import CF, H1, GridFunction, BND, BBND, Integrate

    if drawFunc == None:
        from ngsolve import Draw
        drawFunc = Draw


    val = {bnd:1 for bnd in name.split("|")}


    if old:
        if useBBnd:
            fes = H1(mesh, dirichlet_bbnd=name)
            sol = GridFunction(fes, "bbnd")
            sol.Set(CF([val[bbnd] if bbnd in val.keys() else 0 for bbnd in mesh.GetBBoundaries()]), VOL_or_BND=BBND)
        else:
            fes = H1(mesh, dirichlet=name)
            sol = GridFunction(fes, "bnd")
            sol.Set(CF([val[bnd] if bnd in val.keys() else 0 for bnd in mesh.GetBoundaries()]), VOL_or_BND=BND)
        drawFunc(sol)
        print("-----", name, sum(sol.vec))

    else:
        bnd_CF = CF([val[bnd] if bnd in val.keys() else 0 for bnd in (mesh.GetBoundaries() if not useBBnd else mesh.GetBBoundaries())])
        drawFunc(bnd_CF, mesh, "bnd", draw_vol=False)
        print("-----", name, Integrate(bnd_CF, mesh, VOL_or_BND=(BND if not useBBnd else BBND)))

    if block:
        #cmdInput(locals(), globals())
        input()

#########################################
#ngsglobals.msg_level = 5



#              Geometrie
#
#
#


def MakeGeometry(H_L, H_M, delta_rot, delta_mag, r_Fe, tau, PZ, maxh, faktor_d_rotor):

    print("tau = ", tau)
    print("maxh =", maxh)
    print("maxh_mag = ", maxh * delta_mag)
    d_Fe = faktor_d_rotor*delta_rot
    r_i = r_Fe - d_Fe
    r_L = r_Fe + H_L
    d_phi = 360/PZ
    phi_M = d_phi*tau/2
    subtract = WorkPlane(Axes((0,0,0), n=Z, h=Y)).Rotate(d_phi/2)
    subtract.Line(r_L).Rotate(90)
    subtract.Arc(r_L, (360-d_phi)).Rotate(90)
    subtract.Line(r_L)
    subtract = subtract.Face()

    inner = WorkPlane().Circle(r_i).Face() - subtract
    inner.edges.name="inner"

    rotor = WorkPlane().Circle(r_Fe).Face() #-rect
    rotor.edges.maxh = 2*maxh*delta_rot
    rotor.name = "rotor"
    rotor = rotor - subtract
    rotor.col = (1,0,0)

    magnets =  []#[WorkPlane(Axes((0,0,0), n=Z, h=X))] * PZ
    right = 360/(PZ*2)
    left = 360*(1-1/(PZ*2))
    for i in range(PZ):
        if (i*d_phi < right) or (i*d_phi > left):
            magnets = magnets + [WorkPlane(Axes((0,0,0), n=Z, h=X))]
            magnets[i].MoveTo(r_Fe*sin(i*d_phi*pi/180), r_Fe*cos(i*d_phi*pi/180))
            magnets[i].Arc(r_Fe, -phi_M).Rotate(90)
            magnets[i].Line(H_M).Rotate(90)
            magnets[i].Arc((r_Fe + H_M), 2*phi_M).Rotate(90)
            magnets[i].Line(H_M).Rotate(90)
            magnets[i].Arc(r_Fe, -phi_M).Rotate(-d_phi)
            magnets[i] = magnets[i].Face() #- rect
            magnets[i].maxh = delta_mag*maxh
            magnets[i].edges.name = "magnet_edge"
            magnets[i].name = f"magnets_{i}" #{i:"Anweisungen zur Formatierung"} Formated String
            magnets[i].col = (0, 0, 1)

    """ else:
        magnets= WorkPlane().Circle(r_Fe+H_M).Face()
        magnets.maxh = maxh*delta_mag
        magnets = magnets - subtract
        magnets.name = "magnets"
        magnets.edges.name = "magnet_edge"
        magnets.col = (0,0,1)"""

    outer = WorkPlane().Circle(r_L).Face()
    outer.name = "air"
    outer = outer - subtract
    outer.edges.name = "outer"
    outer.edges

    outer.col = (0,1,0)


    geo = Glue(([outer-rotor - sum(magnets), rotor- inner] + magnets))


    #PERIODIC EDGES
    trafo = Rotation(Axis((0,0,0), Z), -d_phi)

    if tau<1 :
        a = geo.edges.Nearest((-(r_Fe+H_L/2)*sin(d_phi/2*pi/180),(r_Fe+H_L/2)*cos(d_phi/2*pi/180),0))
        a.name = "left_o"
        b = geo.edges.Nearest((-(r_Fe-delta_rot)*sin(d_phi/2*pi/180),(r_Fe-delta_rot)*cos(d_phi/2*pi/180),0))
        b.name = "left_i"
        d = geo.edges.Nearest(((r_Fe+H_L/2)*sin(d_phi/2*pi/180),(r_Fe+H_L/2)*cos(d_phi/2*pi/180),0))
        d.name = "right_o"
        c = geo.edges.Nearest(((r_Fe/2)*sin(d_phi/2*pi/180),(r_Fe/2)*cos(d_phi/2*pi/180),0))
        c.name = "right_i"
    else:
        a = geo.edges.Nearest((-(r_Fe+H_L - (H_L-H_M)/2)*sin(d_phi/2*pi/180),(r_Fe+H_L - (H_L-H_M)/2)*cos(d_phi/2*pi/180),0))
        a.name = "left_o"
        b = geo.edges.Nearest((-(r_Fe-delta_rot)*sin(d_phi/2*pi/180),(r_Fe-delta_rot)*cos(d_phi/2*pi/180),0))
        b.name = "left_i"
        e = geo.edges.Nearest((-(r_Fe+H_M/2)*sin(d_phi/2*pi/180),(r_Fe+H_M/2)*cos(d_phi/2*pi/180),0))
        e.name = "left_m"
        f = geo.edges.Nearest(((r_Fe+H_M/2)*sin(d_phi/2*pi/180),(r_Fe+H_M/2)*cos(d_phi/2*pi/180),0))
        f.name = "right_m"
        d = geo.edges.Nearest(((r_Fe+H_L - (H_L-H_M)/2)*sin(d_phi/2*pi/180),(r_Fe+H_L - (H_L-H_M)/2)*cos(d_phi/2*pi/180),0))
        d.name = "right_o"
        c = geo.edges.Nearest(((r_Fe/2*sin(d_phi/2*pi/180)),(r_Fe/2*cos(d_phi/2*pi/180)),0))
        c.name = "right_i"

        e.Identify(f, "magnet_edge", type = IdentificationType.PERIODIC, trafo=trafo)

    a.Identify(d, "air_edge",type = IdentificationType.PERIODIC, trafo=trafo)
    b.Identify(c, "rotor_edge",type = IdentificationType.PERIODIC, trafo=trafo)

    return geo

sweep = True

# -----------------------------------------------------------------------------
# --- Parameters
# -----------------------------------------------------------------------------

##### WERTE WOLFSCHLUCKNER ######
#PZ = 8
#mu_air = 4e-7*np.pi
#mu_magnet = 1.05*mu_air
#mu_rotor = mu_air*5e2

#sigma_magnet = 8e5
#sigma_rotor =  0
#################################
#Dicke des Rotors ist faktor_d_rotor*delta_rot
f_dr = 8

mu_air = 4e-7*np.pi
mu_magnet = 1.05*mu_air
mu_rotor = mu_air*5e2
sigma_magnet = 8e5
sigma_rotor =  0 #1.86e6

order0 = 3
tau = 1
nu=9
PZ = 8
savetime = 0

K0= 10000
f = 1e6
omega = 2*np.pi*f

delta = lambda omega, sigma, mu : sqrt(2/(nu*omega*sigma*mu))

omega = 2*np.pi*f

delta_rot = delta(omega, 1.86e6, mu_rotor)
delta_mag = delta(omega, sigma_magnet, mu_magnet)
#maxh ist der Koeffizient der max Elementgröße im Rotor. d.h. maxh_rotor = maxh*delta_rot


mu = {"air":mu_air, "rotor": mu_rotor, "magnets.*": mu_magnet }


sigma = {"air":0, "rotor": sigma_rotor, "magnets.*": sigma_magnet}

#A_mag INSGESAMT (d.h. *PZ=8)
A_mags = 0.011519981499672

# Phi berechnung
def Phi(x,y):
    return atan2(y,x)

def K(x,y, nu=1):
     return K0*exp(-1j*nu*Phi(x,y)*PZ/2)

#              Finite Elemente Raum
#
#
#
if(tau<1):
    phase = [-1,-1]
else:
    phase = [-1,-1,-1]

#Frequenz zwischen 0 und 5e5, weil sonst das meshing zu schwierig wird und sinnlos.
x_val = np.logspace(0, 6, 5)
p_values=[]
#p_flat =np.loadtxt(f'sweep_flat_onlymag_{tau}_{nu}_PZ{PZ}.csv', delimiter=',')
#print(len(p_flat))
i=0
with (open(f'sweep_deg_time_{tau}_{nu}.csv', 'w') if tau is 1 and savetime is 1 else nullcontext()) as time_file:
    with open(f'sweep_deg_onlymag_{tau}_{nu}_PZ{PZ}.csv', 'w') as file:
        for freq in x_val:
            i=i+1
            print(f"Starting {i}th Simulation at f = {freq} and nu = {nu}\n")
            omega = 2*np.pi*freq

            delta_rot = delta(omega, 1.86e6, mu_rotor)
            delta_mag = delta(omega, sigma_magnet, mu_magnet)
            print("delta_mag = ", delta_mag)
            print(f"delta_rot = {delta_rot}")

            mp = MeshingParameters(maxh=0.1)
            maxh = 0.5
            if(freq > 1e5):
                maxh = 0.6
            if(freq > 7e5):
                maxh = 0.7
            #Adaptive Meshing
            mesh = Mesh(OCCGeometry(MakeGeometry(H_L=8e-3, H_M=6e-3, delta_rot = delta_rot, delta_mag = delta_mag, r_Fe=302.577e-3, tau=tau, PZ=PZ, maxh = sqrt(nu)*maxh, faktor_d_rotor = f_dr), dim = 2).GenerateMesh(mp=mp))
            mesh.Curve(3)

            muCF = mesh.MaterialCF(mu, default=mu_air)

            sigmaCF = mesh.MaterialCF(sigma, default=0)

            V = Periodic(H1(mesh, order = order0, dirichlet = "inner", complex=True), phase=phase)
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
            print("HERE")
            end_time = time.time()
            elapsed_time = end_time - start_time
            E = -1j * omega * u
            J = sigmaCF * E
            p = E*Conj(J)/2
            losses = Integrate(p, mesh, definedon=mesh.Materials("magnets.*"))*PZ
            print(f"losses are {losses.real} at freq {freq}")
            file.write(f'{freq},{losses.real},{losses.real/A_mags}\n')
            if savetime:
                time_file.write(f'{mesh.ne},{elapsed_time}\n')
            p_values.append(losses.real)

print(len(p_values), "vs ", len(x_val))
#write x_val and p_val to one .txt file
#np.savetxt(f'sweep_flat_onlymag_{tau}_{nu}.csv', np.column_stack((x_val, p_values)), delimiter=',', fmt='%s')
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
