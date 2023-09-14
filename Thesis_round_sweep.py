from netgen.occ import *
from ngsolve import *

#from netgen.webgui import Draw as DrawGeo
#from ngsolve.webgui import Draw
from netgen.meshing import MeshingParameters, IdentificationType
import numpy as np
import math as math
import netgen.gui
import matplotlib.pyplot as plot

#ngsglobals.msg_level = 5



#              Geometrie
#
#
#
def MakeGeometry(H_L, H_M, H_Fe, H_W, Tau, PZ):
    
    print("Tau = ", Tau)

    r_Fe = H_Fe + H_W
    r_L = r_Fe + H_L
    rect = WorkPlane(Axes((-r_L,-r_L,0), n=Z, h=X)).Rectangle(2*r_L, r_L).Face()

    inner = WorkPlane().Circle(H_W).Face() - rect
    inner.edges.name="inner"

    rotor = WorkPlane().Circle(r_Fe).Face() #-rect
    rotor.name = "rotor"
    rotor = rotor - rect
    rotor.col = (1,0,0)

    magnets =  []#[WorkPlane(Axes((0,0,0), n=Z, h=X))] * PZ
    d_phi = 360/PZ
    phi_M = d_phi*Tau/2
    for i in range(PZ):
        if (i*d_phi < 90) or (i*d_phi > 270):
            magnets = magnets + [WorkPlane(Axes((0,0,0), n=Z, h=X))]
            magnets[i].MoveTo(r_Fe*sin(i*d_phi*pi/180), r_Fe*cos(i*d_phi*pi/180))
            magnets[i].Arc(r_Fe, -phi_M).Rotate(90)
            magnets[i].Line(H_M).Rotate(90)
            magnets[i].Arc((r_Fe + H_M), 2*phi_M).Rotate(90)
            magnets[i].Line(H_M).Rotate(90)
            magnets[i].Arc(r_Fe, -phi_M).Rotate(-d_phi)
            magnets[i] = magnets[i].Face() #- rect
            magnets[i].maxh = 0.001
            magnets[i].edges.maxh = 0.0005
            #magnets[i].vertices.maxh = 0.00001 #geht nicht??
            magnets[i].name = f"magnets_{i}" #{i:"Anweisungen zur Formatierung"} Formated String
            magnets[i].col = (0, 0, 1)

       
    outer = WorkPlane().Circle(r_L).Face() 
    outer.name = "air"
    outer = outer - rect
    outer.edges.name = "outer"
    outer.edges

    outer.col = (0,1,0)

    
    
    

    geo = Glue(([outer-rotor - sum(magnets), rotor- inner] + magnets))
    #geo = Glue([geo - rect])
    #for i in range(PZ):
    #    geo = Glue([geo, magnets[i]])
    #PERIODIC EDGES
    trafo = Rotation(Axis((0,0,0), Z), 180)

    if Tau<1 :
        a = geo.edges.Nearest((-(r_Fe+H_L/2),0,0))
        a.maxh = 0.0005
        a.name = "left_o"
        b = geo.edges.Nearest((-(r_Fe/2),0,0)) 
        b.name = "left_i"
        b.maxh = 0.0005
        d = geo.edges.Nearest(((r_Fe+H_L/2),0,0))
        d.name = "right_o"
        d.maxh = 0.0005
        c = geo.edges.Nearest(((r_Fe/2),0,0))
        c.name = "right_i"
        c.maxh = 0.0005
    else:
        a = geo.edges.Nearest((-(r_Fe+H_L - (H_L-H_M)/2),0,0))
        a.name = "left_o"
        b = geo.edges.Nearest((-(r_Fe/2),0,0)) 
        b.name = "left_i"
        e = geo.edges.Nearest((-(r_Fe+H_M/2),0,0))
        e.name = "left_m"
        f = geo.edges.Nearest(((r_Fe+H_M/2),0,0))
        f.name = "right_m"
        d = geo.edges.Nearest(((r_Fe+H_L/2 - (H_L-H_M)/2),0,0))
        d.name = "right_o"
        c = geo.edges.Nearest(((r_Fe/2),0,0))
        c.name = "right_i"

        e.Identify(f, "magnet_edge", type = IdentificationType.PERIODIC, trafo=trafo)
    
    a.Identify(d, "air_edge",type = IdentificationType.PERIODIC, trafo=trafo)
    b.Identify(c, "rotor_edge",type = IdentificationType.PERIODIC, trafo=trafo)
    
    return geo

PZ = 2 
mp = MeshingParameters(maxh=0.4)
#mp.RestrictH(x=0, y=0, z=1, h=0.0025)
phase = [-1,-1]

#radius Rotor ist bei Schmid 35.19mm, r_rot = H_Fe+H_W
geo = MakeGeometry(H_L=8e-3, H_M=6e-3, H_Fe=26.19e-3, H_W=9e-3, Tau=3/4, PZ=PZ)
mesh = Mesh(OCCGeometry(geo, dim = 2).GenerateMesh(mp=mp))
exit()
mesh.Curve(3)


mu_air = 4e-7*np.pi
mu_magnet = 1.08*mu_air
mu_rotor = mu_air*5e3
mu = {"air":mu_air, "rotor": mu_rotor}
[mu.update({f"magnets_{i}": mu_magnet}) for i in range(PZ)]
#version_2 = CF([val[mat] if mat in val.keys() else 0 for mat in mesh.GetMaterials()])
muCF = mesh.MaterialCF(mu, default=mu_air)

sigma = {"air":1e-3, "rotor": 1.86e6}
[sigma.update({f"magnets_{i}": 0.667e6}) for i in range(PZ)]
#sigma = {"air":0, "rotor": 0}
#[sigma.update({f"magnets_{i}": 0}) for i in range(PZ)]
sigmaCF = mesh.MaterialCF(sigma, default=0)
sigma_v = {"air":None, "rotor": None}
[sigma_v.update({f"magnets_{i}": 0.667e6}) for i in range(PZ)]
sigma_visual = mesh.MaterialCF(sigma_v)

versions = [muCF, sigmaCF]

Br = 1
K0= 10000
omega = 50

def K(x,y):
        K1 = -K0*cos(Phi(x,y))
        K2 = -K0*cos(Phi(x,y) + 2*pi/3)*exp(1j*2*pi/3)
        K3 = -K0*cos(Phi(x,y) + 4*pi/3)*exp(1j*4*pi/3)
        return K1 + K2 + K3   

# Phi berechnung
def Phi(x,y):
    return atan2(y,x)


V = Periodic(H1(mesh, order = 3, dirichlet = "inner", complex=True), phase=phase)
trial = V.TrialFunction()
test = V.TestFunction()

sweep = False

x_val = np.logspace(0, 6, 100)
p_values=[]
for omega in x_val:

    if(sweep==True):
         continue
    else:
         omega = 50
    u = GridFunction(V) 

    a = BilinearForm(V, symmetric = True)
    a +=  1/muCF*grad(trial)*grad(test) * dx
    a += 1j*omega*sigmaCF*test * trial * dx#("rotor|magnet|air") 1j*

    c = Preconditioner(a, type="direct", inverse = "umfpack")    

    f = LinearForm(V)
    f += K(x,y)*test.Trace()*ds("outer")
    with TaskManager():
            a.Assemble()
            f.Assemble()
            bvp = BVP(bf=a, lf=f, gf=u, pre=c)
            bvp.Do()
    
    p = sigmaCF*omega*omega*u*Conj(u)/2
    energy = Integrate(p, mesh, definedon=mesh.Materials("magnets_0"))*2
    if(sweep!=True):
        print("Verluste = ", energy)
        break

    p_values.append(energy.real)



if(sweep==True):
    p_flat =np.loadtxt('p_flat.csv', delimiter=',')
    plot.plot(x_val, p_values)
    plot.title('Eddy Current 1st O. Losses Half-Round, Tau=3/4, PZ=2')
    plot.xscale('log')
    plot.yscale('log')
    plot.xlabel('Frequency (Hz)')
    plot.ylabel('Power Losses (W/m)')
    plot.grid(True, which='both')
    plot.savefig('Losses_Round_180_Deg.jpg', format='jpg')
    plot.show()

    r_vs_f=[]
    for i in np.arange(0,len(x_val)):
        d = p_values[i].real/p_flat[i]
        r_vs_f.append(d)

    plot.semilogx(x_val, r_vs_f)
    plot.title('Eddy Current 1st Order Losses Ratio Flat vs Half-Round, Tau=3/4, PZ=2')
    plot.xlabel('Frequency (Hz)')
    plot.ylabel('Ratio Power Losses (W/m)')
    plot.grid(True, which='both')
    plot.savefig('Ratio_1st_Order_180_Deg.jpg', format='jpg')
    plot.show()
