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
PZ = 2

def MakeGeometry(H_L, H_M, delta, r_Fe, tau, PZ, maxh):

    print("tau = ", tau)
    print("maxh = ", maxh * delta)
    d_Fe = 5*delta
    r_i = r_Fe - d_Fe
    r_L = r_Fe + H_L
    d_phi = 360/PZ
    phi_M = d_phi*tau/2
    #rect = WorkPlane(Axes((-r_L,-r_L,0), n=Z, h=X)).Rectangle(2*r_L, r_L).Face()
    subtract = WorkPlane(Axes((0,0,0), n=Z, h=Y)).Rotate(d_phi/2)
    #subtract.MoveTo((r_Fe+H_M)*sin(d_phi/2*pi/180), (r_Fe+H_M)*cos(d_phi/2*pi/180))
    subtract.Line(r_L).Rotate(90)
    subtract.Arc(r_L, (360-d_phi)).Rotate(90)
    subtract.Line(r_L)
    #subtract.Line(r_L).Rotate(180+d_phi)
    #subtract.Line(r_L)
    subtract = subtract.Face()

    inner = WorkPlane().Circle(r_i).Face() - subtract
    inner.edges.name="inner"

    rotor = WorkPlane().Circle(r_Fe).Face() #-rect
    rotor.maxh = maxh*delta
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
            magnets[i].maxh = 0.2e-3
            magnets[i].edges.name = "magnet_edge"
            #magnets[i].vertices.maxh = 0.00001 #geht nicht??
            magnets[i].name = f"magnets_{i}" #{i:"Anweisungen zur Formatierung"} Formated String
            magnets[i].col = (0, 0, 1)


    outer = WorkPlane().Circle(r_L).Face()
    outer.name = "air"
    outer = outer - subtract
    outer.edges.name = "outer"
    outer.edges

    outer.col = (0,1,0)





    geo = Glue(([outer-rotor - sum(magnets), rotor- inner] + magnets))

    #geo = Glue([geo - rect])
    #for i in range(PZ):
    #    geo = Glue([geo, magnets[i]])
    #PERIODIC EDGES
    trafo = Rotation(Axis((0,0,0), Z), -d_phi)

    if tau<1 :
        a = geo.edges.Nearest((-(r_Fe+H_L/2)*sin(d_phi/2*pi/180),(r_Fe+H_L/2)*cos(d_phi/2*pi/180),0))
        a.name = "left_o"
        b = geo.edges.Nearest((-(r_Fe/2)*sin(d_phi/2*pi/180),(r_Fe/2)*cos(d_phi/2*pi/180),0))
        b.name = "left_i"
        d = geo.edges.Nearest(((r_Fe+H_L/2)*sin(d_phi/2*pi/180),(r_Fe+H_L/2)*cos(d_phi/2*pi/180),0))
        d.name = "right_o"
        c = geo.edges.Nearest(((r_Fe/2)*sin(d_phi/2*pi/180),(r_Fe/2)*cos(d_phi/2*pi/180),0))
        c.name = "right_i"
    else:
        a = geo.edges.Nearest((-(r_Fe+H_L - (H_L-H_M)/2)*sin(d_phi/2*pi/180),(r_Fe+H_L - (H_L-H_M)/2)*cos(d_phi/2*pi/180),0))
        a.name = "left_o"
        b = geo.edges.Nearest((-(r_Fe/2)*sin(d_phi/2*pi/180),(r_Fe/2)*cos(d_phi/2*pi/180),0))
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
# --- Parameters
# -----------------------------------------------------------------------------

PZ = 10

mu_air = 4e-7*np.pi
mu_magnet = 1.05*mu_air
mu_rotor = mu_air*5e2

sigma_magnet = 8e5
sigma_rotor =  1.86e6

order0 = 2
tau = 3/4

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
#print(MakeGeometry(H_L=8e-3, H_M=6e-3, H_Fe=20e-3, H_W=30e-3, Tau=1, PZ=PZ))
mp = MeshingParameters(maxh=0.4)
#mp.RestrictH(x=0, y=0, z=1, h=0.0025)
mesh = Mesh(OCCGeometry(MakeGeometry(H_L=8e-3, H_M=6e-3, delta = 0.3e-3, r_Fe = r_Fe, tau=tau, PZ=PZ, maxh = 0.5), dim = 2).GenerateMesh(mp=mp))
mesh.Curve(3)

#Materials
#('air', 'rotor', 'magnets_0', 'magnets_1', 'magnets_2', 'magnets_3', 'magnets_4', 'magnets_5', 'magnets_6', 'magnets_7', 'magnets_8', 'magnets_9', 'magnets_10', 'magnets_11')
#Boundaries
#{'default', 'inner', 'outer'}

print(mesh.GetMaterials())
print(set(mesh.GetBoundaries()))
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
muCF = mesh.MaterialCF(mu, default=mu_air)



sigma = {"air":0, "rotor": sigma_rotor, "magnets.*": sigma_magnet}
sigmaCF = mesh.MaterialCF(sigma, default=0)
versions = [muCF, sigmaCF]
[Draw(versions[i], mesh, str(i)) for i in range(len(versions))]
Draw(sigmaCF, mesh, "sigma_v")

regions = mesh.MaterialCF({"air":1, "rotor": 2, "magnets.*": 3})
Draw(regions, mesh, "regions")
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

#              Finite Elemente Raum
#
#
#
if(tau<1):
    phase = [-1,-1]
else:
    phase = [-1,-1,-1]



# -----------------------------------------------------------------------------
#  WIRBELSTROMVERLUSTE
# -----------------------------------------------------------------------------
