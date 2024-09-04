from netgen.occ import *
from ngsolve import *

#from netgen.webgui import Draw as DrawGeo
#from ngsolve.webgui import Draw
from netgen.meshing import MeshingParameters, IdentificationType
import numpy as np
import math as math
import netgen.gui


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

def MakeGeometry(H_L, H_M, maxh_rotor, maxh_mag, r_Fe, tau, PZ, maxh):
    print("maxh = ", maxh)
    print("tau = ", tau)
    print("maxh_mag = ", maxh * delta_mag)
    d_Fe = 8*maxh_rotor
    r_i = r_Fe - d_Fe
    r_L = r_Fe + H_L
    d_phi = 360/PZ
    phi_M = d_phi*tau/2
    print("r_Fe = ", r_Fe)
    print("r_L = ", r_L)
    print("r_M = ", r_Fe + H_M)
    print("r_i = ", r_i)
    subtract = WorkPlane(Axes((0,0,0), n=Z, h=Y)).Rotate(d_phi/2)
    subtract.Line(r_L).Rotate(90)
    subtract.Arc(r_L, (360-d_phi)).Rotate(90)
    subtract.Line(r_L)
    subtract = subtract.Face()

    inner = WorkPlane().Circle(r_i).Face() - subtract
    inner.edges.name="inner"

    rotor = WorkPlane().Circle(r_Fe).Face() #-rect
    rotor.maxh = maxh_rotor
    rotor.name = "rotor"
    rotor = rotor - subtract
    rotor.col = (1,0,0)

    print("ri = ", r_i)
    print("rFe = ", r_Fe)
    print("rL = ", r_L)
    print("r_mag", r_Fe+H_M)
    magnets =  []
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
            magnets[i].maxh = maxh_mag
            magnets[i].edges.name = "magnet_edge"
            magnets[i].name = f"magnets_{i}"
            magnets[i].col = (0, 0, 1)


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


# -----------------------------------------------------------------------------
# ---------------------Parameters
# -----------------------------------------------------------------------------

mu_air = 4e-7*np.pi
mu_magnet = 1.05*mu_air
mu_rotor = mu_air*5e2
sigma_magnet = 8e5
sigma_rotor =  0

order0 = 3
tau = 5/6
nu=1
PZ = 8

K0= 10000
f = 1e6
omega = 2*np.pi*f

delta = lambda omega, sigma, mu : sqrt(2/(nu*omega*sigma*mu))


delta_rot = delta(omega, 1.86e6, mu_rotor)
delta_mag = delta(omega, sigma_magnet, mu_magnet)
print("Delta_mag ist ", delta_mag)

print("delta_rot", delta_rot)
print(f"Frequenz = {f} und u = {nu}\n")
r_Fe = 149.7887453682e-3
maxh =  0.5
maxh_rotor = max(maxh*delta_rot, 1e-4)
maxh_mag = 4*maxh_rotor
mp = MeshingParameters(maxh = maxh_mag)
mesh = Mesh(OCCGeometry(MakeGeometry(H_L=8e-3, H_M=6e-3, maxh_rotor=maxh_rotor, maxh_mag=maxh_mag, r_Fe = r_Fe, tau=tau, PZ=PZ, maxh =  maxh), dim = 2).GenerateMesh())
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
    return atan2(y,x)+np.pi/2

def K(x,y, nu=1):
     return K0*exp(-1j*nu*PZ*Phi(x,y)/2)

# -----------------------------------------------------------------------------
# ---------------------Finite Elemente Raum
# -----------------------------------------------------------------------------

if(tau<1):
    phase = [-1,-1]
else:
    phase = [-1, -1, -1]
V = Periodic(H1(mesh, order = order0, dirichlet = "inner", complex=True), phase=phase)
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
Draw(sigmaCF, mesh, 'materials')
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

try:
    A = Integrate(1, mesh, definedon=mesh.Materials("magnets.*"))
    print("Fläche = ", PZ*A)
    p = E*Conj(J)/2
    #losses = Integrate(p, mesh, definedon=mesh.Materials("magnets.*"))
    losses = Integrate(p, mesh, definedon=mesh.Materials("magnets.*"))

    print("P(u, u) = ", PZ*losses)
    print("P/omega = ", PZ*losses/omega)

    with open("simulations.txt", "a") as file:
        file.write(f"Periodic Sim: PZ = {PZ}, Tau = {tau}, omega = {omega}, A_magnet = {A*PZ}: losses = {PZ*losses.real}, delta_magnet = {delta_mag}\n")
        print("written to file")
except Exception as e:
    print("An error occurred: ", e)