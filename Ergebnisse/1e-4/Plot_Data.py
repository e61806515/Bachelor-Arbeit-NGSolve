import matplotlib.pyplot as plot
import numpy as np
import matplotlib.font_manager as fm
from matplotlib.ticker import LogLocator

fm.fontManager.addfont(r'C:\Users\malbu\Desktop\Computer Modern\Serif\cmunrm.ttf')
font_prop = fm.FontProperties(fname=r'C:\Users\malbu\Desktop\Computer Modern\Serif\cmunrm.ttf')
print(font_prop.get_name())
deg1_1_80 = "sweep_deg_onlymag_1_1_PZ8_80samples.csv"
deg3_1_80 = "sweep_deg_onlymag_1_3_PZ8_80samples.csv"
deg5_1_80 = "sweep_deg_onlymag_1_5_PZ8_80samples.csv"
deg7_1_80 = "sweep_deg_onlymag_1_7_PZ8_80samples.csv"
deg9_1_80 = "sweep_deg_onlymag_1_9_PZ8_80samples.csv"

deg1_1_MAGROT = "sweep_deg_MAGROT_1_1_PZ8_60samples.csv"
deg3_1_MAGROT = "sweep_deg_MAGROT_1_3_PZ8_60samples.csv"
deg5_1_MAGROT = "sweep_deg_MAGROT_1_5_PZ8_60samples.csv"
deg7_1_MAGROT = "sweep_deg_MAGROT_1_7_PZ8_60samples.csv"
deg9_1_MAGROT = "sweep_deg_MAGROT_1_9_PZ8_60samples.csv"

deg1_1 = "sweep_deg_onlymag_1_1_PZ8_60samples_149mm.csv"
deg3_1 = "sweep_deg_onlymag_1_3_PZ8_60samples_149mm.csv"
deg5_1 = "sweep_deg_onlymag_1_5_PZ8_60samples_149mm.csv"
deg7_1 = "sweep_deg_onlymag_1_7_PZ8_60samples_149mm.csv"
deg9_1 = "sweep_deg_onlymag_1_9_PZ8_60samples_149mm.csv"

flatA1_1 = "sweep_flat_A_onlymag_1_1_PZ8_60samples_149mm.csv"
flatA3_1 = "sweep_flat_A_onlymag_1_3_PZ8_60samples_149mm.csv"
flatA5_1 = "sweep_flat_A_onlymag_1_5_PZ8_60samples_149mm.csv"
flatA7_1 = "sweep_flat_A_onlymag_1_7_PZ8_60samples_149mm.csv"
flatA9_1 = "sweep_flat_A_onlymag_1_9_PZ8_60samples_149mm.csv"

flatT1_1 = "sweep_flat_tau_onlymag_1_1_PZ8_60samples.csv"
flatT3_1 = "sweep_flat_tau_onlymag_1_3_PZ8_60samples.csv"
flatT5_1 = "sweep_flat_tau_onlymag_1_5_PZ8_60samples.csv"
flatT7_1 = "sweep_flat_tau_onlymag_1_7_PZ8_60samples.csv"
flatT9_1 = "sweep_flat_tau_onlymag_1_9_PZ8_60samples.csv"

deg1_56_MAGROT = "sweep_deg_MAGROT_0.8333333333333334_1_PZ8_60samples.csv"
deg3_56_MAGROT = "sweep_deg_MAGROT_0.8333333333333334_3_PZ8_60samples.csv"
deg5_56_MAGROT = "sweep_deg_MAGROT_0.8333333333333334_5_PZ8_60samples.csv"
deg7_56_MAGROT = "sweep_deg_MAGROT_0.8333333333333334_7_PZ8_60samples.csv"
deg9_56_MAGROT = "sweep_deg_MAGROT_0.8333333333333334_9_PZ8_60samples.csv"

deg1_56 = "sweep_deg_onlymag_0.8333333333333334_1_PZ8_60samples_149mm.csv"
deg3_56 = "sweep_deg_onlymag_0.8333333333333334_3_PZ8_60samples_149mm.csv"
deg5_56 = "sweep_deg_onlymag_0.8333333333333334_5_PZ8_60samples_149mm.csv"
deg7_56 = "sweep_deg_onlymag_0.8333333333333334_7_PZ8_60samples_149mm.csv"
deg9_56 = "sweep_deg_onlymag_0.8333333333333334_9_PZ8_60samples_149mm.csv"

flatA1_56 = "sweep_flat_A_onlymag_0.8333333333333334_1_PZ8_60samples_149mm.csv"
flatA3_56 = "sweep_flat_A_onlymag_0.8333333333333334_3_PZ8_60samples_149mm.csv"
flatA5_56 = "sweep_flat_A_onlymag_0.8333333333333334_5_PZ8_60samples_149mm.csv"
flatA7_56 = "sweep_flat_A_onlymag_0.8333333333333334_7_PZ8_60samples_149mm.csv"
flatA9_56 = "sweep_flat_A_onlymag_0.8333333333333334_9_PZ8_60samples_149mm.csv"

flatT1_56 = "sweep_flat_tau_onlymag_0.8333333333333334_1_PZ8_60samples.csv"
flatT3_56 = "sweep_flat_tau_onlymag_0.8333333333333334_3_PZ8_60samples.csv"
flatT5_56 = "sweep_flat_tau_onlymag_0.8333333333333334_5_PZ8_60samples.csv"
flatT7_56 = "sweep_flat_tau_onlymag_0.8333333333333334_7_PZ8_60samples.csv"
flatT9_56 = "sweep_flat_tau_onlymag_0.8333333333333334_9_PZ8_60samples.csv"

deg1_1_144mm = "sweep_deg_onlymag_1_1_PZ8_60samples_144.5mm.csv"

flatA1_1_144mm = "sweep_flat_A_onlymag_1_1_PZ8_60samples_144.5mm.csv"

flatT1_1_144mm = "sweep_flat_tau_onlymag_1_1_PZ8_60samples_144.5mm.csv"

#deg1_rot_1e3 = "sweep_deg_rotor_1_1_PZ8_40samples_1e-3.csv"
#deg1_rot_1e4 = "sweep_deg_rotor_1_1_PZ8_40samples_1e-4.csv"

#load frequencies
freqs = np.loadtxt(deg1_1, delimiter=',', usecols=(0,))
freqs80 = np.loadtxt(deg1_1_80, delimiter=',', usecols=(0,))
#freqs40 = np.loadtxt(deg1_rot_1e3, delimiter=',', usecols=(0,))

deg_1_1 = np.loadtxt(deg1_1, delimiter=',', usecols=(1,))
deg_1_3 = np.loadtxt(deg3_1, delimiter=',', usecols=(1,))
deg_1_5 = np.loadtxt(deg5_1, delimiter=',', usecols=(1,))
deg_1_7 = np.loadtxt(deg7_1, delimiter=',', usecols=(1,))
deg_1_9 = np.loadtxt(deg9_1, delimiter=',', usecols=(1,))

deg_1_1_MAGROT = np.loadtxt(deg1_1_MAGROT, delimiter=',', usecols=(1,))
deg_1_3_MAGROT = np.loadtxt(deg3_1_MAGROT, delimiter=',', usecols=(1,))
deg_1_5_MAGROT = np.loadtxt(deg5_1_MAGROT, delimiter=',', usecols=(1,))
deg_1_7_MAGROT = np.loadtxt(deg7_1_MAGROT, delimiter=',', usecols=(1,))
deg_1_9_MAGROT = np.loadtxt(deg9_1_MAGROT, delimiter=',', usecols=(1,))

deg_1_1_80 = np.loadtxt(deg1_1_80, delimiter=',', usecols=(1,))
deg_1_3_80 = np.loadtxt(deg3_1_80, delimiter=',', usecols=(1,))
deg_1_5_80 = np.loadtxt(deg5_1_80, delimiter=',', usecols=(1,))
deg_1_7_80 = np.loadtxt(deg7_1_80, delimiter=',', usecols=(1,))
deg_1_9_80 = np.loadtxt(deg9_1_80, delimiter=',', usecols=(1,))

flatA_1_1 = np.loadtxt(flatA1_1, delimiter=',', usecols=(1,))
flatA_1_3 = np.loadtxt(flatA3_1, delimiter=',', usecols=(1,))
flatA_1_5 = np.loadtxt(flatA5_1, delimiter=',', usecols=(1,))
flatA_1_7 = np.loadtxt(flatA7_1, delimiter=',', usecols=(1,))
flatA_1_9 = np.loadtxt(flatA9_1, delimiter=',', usecols=(1,))

flatT_1_1 = np.loadtxt(flatT1_1, delimiter=',', usecols=(1,))
flatT_1_3 = np.loadtxt(flatT3_1, delimiter=',', usecols=(1,))
flatT_1_5 = np.loadtxt(flatT5_1, delimiter=',', usecols=(1,))
flatT_1_7 = np.loadtxt(flatT7_1, delimiter=',', usecols=(1,))
flatT_1_9 = np.loadtxt(flatT9_1, delimiter=',', usecols=(1,))

deg_56_1 = np.loadtxt(deg1_56, delimiter=',', usecols=(1,))
deg_56_3 = np.loadtxt(deg3_56, delimiter=',', usecols=(1,))
deg_56_5 = np.loadtxt(deg5_56, delimiter=',', usecols=(1,))
deg_56_7 = np.loadtxt(deg7_56, delimiter=',', usecols=(1,))
deg_56_9 = np.loadtxt(deg9_56, delimiter=',', usecols=(1,))

deg_56_1_MAGROT = np.loadtxt(deg1_56_MAGROT, delimiter=',', usecols=(1,))
deg_56_3_MAGROT = np.loadtxt(deg3_56_MAGROT, delimiter=',', usecols=(1,))
deg_56_5_MAGROT = np.loadtxt(deg5_56_MAGROT, delimiter=',', usecols=(1,))
deg_56_7_MAGROT = np.loadtxt(deg7_56_MAGROT, delimiter=',', usecols=(1,))
deg_56_9_MAGROT = np.loadtxt(deg9_56_MAGROT, delimiter=',', usecols=(1,))

flatA_56_1 = np.loadtxt(flatA1_56, delimiter=',', usecols=(1,))
flatA_56_3 = np.loadtxt(flatA3_56, delimiter=',', usecols=(1,))
flatA_56_5 = np.loadtxt(flatA5_56, delimiter=',', usecols=(1,))
flatA_56_7 = np.loadtxt(flatA7_56, delimiter=',', usecols=(1,))
flatA_56_9 = np.loadtxt(flatA9_56, delimiter=',', usecols=(1,))

flatT_56_1 = np.loadtxt(flatT1_56, delimiter=',', usecols=(1,))
flatT_56_3 = np.loadtxt(flatT3_56, delimiter=',', usecols=(1,))
flatT_56_5 = np.loadtxt(flatT5_56, delimiter=',', usecols=(1,))
flatT_56_7 = np.loadtxt(flatT7_56, delimiter=',', usecols=(1,))
flatT_56_9 = np.loadtxt(flatT9_56, delimiter=',', usecols=(1,))

deg_1_1_144mm = np.loadtxt(deg1_1_144mm, delimiter=',', usecols=(1,))

flatA_1_1_144mm = np.loadtxt(flatA1_1_144mm, delimiter=',', usecols=(1,))

flatT_1_1_144mm = np.loadtxt(flatT1_1_144mm, delimiter=',', usecols=(1,))

deg_1_1_144mm_A = np.loadtxt(deg1_1_144mm, delimiter=',', usecols=(2,))

flatT_1_1_144mm_A = np.loadtxt(flatT1_1_144mm, delimiter=',', usecols=(2,))

with open('Ratio_deg_flatA_PZ8_1_1_149mm.csv', 'w') as f:
    for i in range(len(freqs)):
        f.write(str(freqs[i]) + ',' + str(8*deg_1_1[i]/flatA_1_1[i]) + '\n')
    f.close()
with open('Ratio_deg_flatA_PZ8_1_3_149mm.csv', 'w') as f:
    for i in range(len(freqs)):
        f.write(str(freqs[i]) + ',' + str(8*deg_1_3[i]/flatA_1_3[i]) + '\n')
    f.close()
with open('Ratio_deg_flatA_PZ8_1_5_149mm.csv', 'w') as f:
    for i in range(len(freqs)):
        f.write(str(freqs[i]) + ',' + str(8*deg_1_5[i]/flatA_1_5[i]) + '\n')
    f.close()
with open('Ratio_deg_flatA_PZ8_1_7_149mm.csv', 'w') as f:
    for i in range(len(freqs)):
        f.write(str(freqs[i]) + ',' + str(8*deg_1_7[i]/flatA_1_7[i]) + '\n')
    f.close()
with open('Ratio_deg_flatA_PZ8_1_9_149mm.csv', 'w') as f:
    for i in range(len(freqs)):
        f.write(str(freqs[i]) + ',' + str(8*deg_1_9[i]/flatA_1_9[i]) + '\n')
    f.close()
with open('Ratio_deg_flatA_PZ8_0.8333333333333334_1_149mm.csv', 'w') as f:
    for i in range(len(freqs)):
        f.write(str(freqs[i]) + ',' + str(8*deg_56_1[i]/flatA_56_1[i]) + '\n')
    f.close()
with open('Ratio_deg_flatA_PZ8_0.8333333333333334_3_149mm.csv', 'w') as f:
    for i in range(len(freqs)):
        f.write(str(freqs[i]) + ',' + str(8*deg_56_3[i]/flatA_56_3[i]) + '\n')
    f.close()
with open('Ratio_deg_flatA_PZ8_0.8333333333333334_5_149mm.csv', 'w') as f:
    for i in range(len(freqs)):
        f.write(str(freqs[i]) + ',' + str(8*deg_56_5[i]/flatA_56_5[i]) + '\n')
    f.close()
with open('Ratio_deg_flatA_PZ8_0.8333333333333334_7_149mm.csv', 'w') as f:
    for i in range(len(freqs)):
        f.write(str(freqs[i]) + ',' + str(8*deg_56_7[i]/flatA_56_7[i]) + '\n')
    f.close()
with open('Ratio_deg_flatA_PZ8_0.8333333333333334_9_149mm.csv', 'w') as f:
    for i in range(len(freqs)):
        f.write(str(freqs[i]) + ',' + str(8*deg_56_9[i]/flatA_56_9[i]) + '\n')
    f.close()
#Relative Fehler
with open('RelError_deg_flatA_PZ8_1_1_144mm.csv', 'w') as f:
    for i in range(len(freqs)):
        f.write(str(freqs[i]) + ',' + str(abs((flatA_1_1_144mm[i]-deg_1_1_144mm[i])/deg_1_1_144mm[i])) + '\n')
    f.close()

with open('RelError_deg_flatT_PZ8_1_1_144mm.csv', 'w') as f:
    for i in range(len(freqs)):
        f.write(str(freqs[i]) + ',' + str((flatT_1_1_144mm[i]/8-deg_1_1_144mm[i])/deg_1_1_144mm[i]) + '\n')
    f.close()

with open('RelErrorA_deg_flatT_PZ8_1_1_144mm.csv', 'w') as f:
    for i in range(len(freqs)):
        f.write(str(freqs[i]) + ',' + str((flatT_1_1_144mm_A[i]-deg_1_1_144mm_A[i])/deg_1_1_144mm_A[i]) + '\n')
    f.close()

#Berechnung der Flächen
PZ = 8
r_Fe = 144.5e-3
d_M = 6e-3
d_L = 2e-3
A_mag = np.pi*((r_Fe+d_M)**2 - r_Fe**2)
print(f"Magnetfläche des gekrümmten Problems, über r_Fe hergeleitet ist {A_mag/8}")

tau_p = 60*d_L
A_flat = tau_p*d_M*PZ
print(f"Fläche des flachen Problems, über tau_p hergeleitet ist {A_flat/8}")

#alternativ tau_p herleiten über Kreisumfang
U = (r_Fe+d_M/2)*2*np.pi
tau_p = U/PZ
A_flat = tau_p*d_M*PZ
print(f"Fläche des flachen Problems, über U hergeleitet ist {A_flat/8}")
exit()
#deg_rot_o1_1e3 = np.loadtxt(deg1_rot_1e3, delimiter=',', usecols=(1,))
#deg_rot_o1_1e4 = np.loadtxt(deg1_rot_1e4, delimiter=',', usecols=(1,))

#flat_o1 = np.loadtxt(flat1, delimiter=',', usecols=(1,))
#flat_o3 = np.loadtxt(flat3, delimiter=',', usecols=(1,))
print(fm.findfont("CMU Serif"))
plot.rcParams.update({
    "font.family": "serif",          # Setze die Schriftartfamilie auf Serif
    "font.serif": ["CMU Serif", "DejaVu Serif"],     # Beispiel für den Namen der Schriftart
    "text.usetex": False             # Mathtext verwenden
})

plot.plot(freqs, deg_1_1, color='black', label=r'$\nu = 1$')
plot.plot(freqs, deg_1_3, color='red', label=r'$\nu = 3$')
plot.plot(freqs, deg_1_5, color='green', label=r'$\nu = 5$')
plot.plot(freqs, deg_1_7, color='orange', label=r'$\nu = 7$')
plot.plot(freqs, deg_1_9, color='purple', label=r'$\nu = 9$')
plot.legend(loc = 'lower right')
plot.xscale('log')
plot.yscale('log')
plot.grid(True, which="both", linestyle=':', linewidth=0.75)
plot.minorticks_on()
plot.savefig('Losses_tau=1_deg_onlymags_PZ8_.pdf', format='pdf')
plot.show()

plot.plot(freqs, flatA_1_1/8, color='black', label=r'$\nu = 1$')
plot.plot(freqs, flatA_1_3/8, color='red', label=r'$\nu = 3$')
plot.plot(freqs, flatA_1_5/8, color='green', label=r'$\nu = 5$')
plot.plot(freqs, flatA_1_7/8, color='orange', label=r'$\nu = 7$')
plot.plot(freqs, flatA_1_9/8, color='purple', label=r'$\nu = 9$')
plot.legend(loc = 'lower right')
plot.xscale('log')
plot.yscale('log')
plot.grid(True, which="both", linestyle=':', linewidth=0.75)
plot.minorticks_on()
plot.savefig('Losses_tau=1_flatA_onlymags_PZ8_.pdf', format='pdf')
plot.show()

plot.plot(freqs, flatT_1_1/8, color='black', label=r'$\nu = 1$')
plot.plot(freqs, flatT_1_3/8, color='red', label=r'$\nu = 3$')
plot.plot(freqs, flatT_1_5/8, color='green', label=r'$\nu = 5$')
plot.plot(freqs, flatT_1_7/8, color='orange', label=r'$\nu = 7$')
plot.plot(freqs, flatT_1_9/8, color='purple', label=r'$\nu = 9$')
plot.legend(loc = 'lower right')
plot.xscale('log')
plot.yscale('log')
plot.grid(True, which="both", linestyle=':', linewidth=0.75)
plot.minorticks_on()
plot.savefig('Losses_tau=1_flatT_onlymags_PZ8_.pdf', format='pdf')
plot.show()


plot.plot(freqs80, deg_1_1_80/8, color='black', label=r'$\nu = 1$')
plot.plot(freqs80, deg_1_3_80/8, color='red', label=r'$\nu = 3$')
plot.plot(freqs80, deg_1_5_80/8, color='green', label=r'$\nu = 5$')
plot.plot(freqs80, deg_1_7_80/8, color='orange', label=r'$\nu = 7$')
plot.plot(freqs80, deg_1_9_80/8, color='purple', label=r'$\nu = 9$')
plot.legend(loc = 'lower right')
plot.xscale('log')
plot.yscale('log')
#plot.xlabel('Frequenz (Hz)')
#plot.ylabel('Verluste (W/m)')
plot.grid(True, which="both", linestyle=':', linewidth=0.75)
plot.minorticks_on()
plot.savefig('Losses_tau=1_deg_onlymags_PZ8_80SAMPLES.pdf', format='pdf')
plot.show()

plot.plot(freqs, deg_56_1, color='black', label=r'$\nu = 1$')
plot.plot(freqs, deg_56_3, color='red', label=r'$\nu = 3$')
plot.plot(freqs, deg_56_5, color='green', label=r'$\nu = 5$')
plot.plot(freqs, deg_56_7, color='orange', label=r'$\nu = 7$')
plot.plot(freqs, deg_56_9, color='purple', label=r'$\nu = 9$')
plot.legend(loc = 'lower right')
plot.xscale('log')
plot.yscale('log')
plot.grid(True, which="both", linestyle=':', linewidth=0.75)
plot.minorticks_on()
plot.savefig('Losses_tau=0.833_deg_onlymags_PZ8_.pdf', format='pdf')
plot.show()

plot.plot(freqs, flatA_56_1/8, color='black', label=r'$\nu = 1$')
plot.plot(freqs, flatA_56_3/8, color='red', label=r'$\nu = 3$')
plot.plot(freqs, flatA_56_5/8, color='green', label=r'$\nu = 5$')
plot.plot(freqs, flatA_56_7/8, color='orange', label=r'$\nu = 7$')
plot.plot(freqs, flatA_56_9/8, color='purple', label=r'$\nu = 9$')
plot.legend(loc = 'lower right')
plot.xscale('log')
plot.yscale('log')
plot.grid(True, which="both", linestyle=':', linewidth=0.75)
# Haupt- und Nebengitter-Lokatoren festlegen
ax = plot.gca()  # aktuelle Achse abrufen
ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs='auto', numticks=100))
ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
ax.yaxis.set_minor_locator(LogLocator(base=10.0,subs=[2, 3, 4, 5, 6, 7, 8, 9], numticks=100))
plot.minorticks_on()
plot.savefig('Losses_tau=0.833_flatA_onlymags_PZ8_.pdf', format='pdf')
plot.show()

plot.plot(freqs, flatT_56_1/8, color='black', label=r'$\nu = 1$')
plot.plot(freqs, flatT_56_3/8, color='red', label=r'$\nu = 3$')
plot.plot(freqs, flatT_56_5/8, color='green', label=r'$\nu = 5$')
plot.plot(freqs, flatT_56_7/8, color='orange', label=r'$\nu = 7$')
plot.plot(freqs, flatT_56_9/8, color='purple', label=r'$\nu = 9$')
plot.legend(loc = 'lower right')
plot.xscale('log')
plot.yscale('log')
plot.grid(True, which="both", linestyle=':', linewidth=0.75)
plot.minorticks_on()
plot.savefig('Losses_tau=0.833_flatT_onlymags_PZ8_.pdf', format='pdf')
plot.show()

#Ratio between flatA and deg Tau=1
plot.plot(freqs, flatA_1_1/(8*deg_1_1), color='black', label=r'$\nu = 1$')
plot.plot(freqs, flatA_1_3/(8*deg_1_3), color='red', label=r'$\nu = 3$')
plot.plot(freqs, flatA_1_5/(8*deg_1_5), color='green', label=r'$\nu = 5$')
plot.plot(freqs, flatA_1_7/(8*deg_1_7), color='orange', label=r'$\nu = 7$')
plot.plot(freqs, flatA_1_9/(8*deg_1_9), color='purple', label=r'$\nu = 9$')
plot.legend(loc = 'lower right')
plot.xscale('log')
plot.ylim(0.9, 1.1)
plot.grid(True, which="both", linestyle=':', linewidth=0.75)
plot.minorticks_on()
plot.savefig('Ratio_flatA_deg_PZ8_Tau=1.pdf', format='pdf')
plot.show()

'''
#Ratio between flatT and deg Tau=1
plot.plot(freqs, flatT_1_1/deg_1_1, color='black', label=r'$\nu = 1$')
plot.plot(freqs, flatT_1_3/deg_1_3, color='red', label=r'$\nu = 3$')
plot.plot(freqs, flatT_1_5/deg_1_5, color='green', label=r'$\nu = 5$')
plot.plot(freqs, flatT_1_7/(8*deg_1_7), color='orange', label=r'$\nu = 7$')
plot.plot(freqs, flatT_1_9/deg_1_9, color='purple', label=r'$\nu = 9$')
plot.legend(loc = 'lower right')
plot.xscale('log')
plot.ylim(0.5, 0.8)
plot.grid(True, which="both", linestyle=':', linewidth=0.75)
plot.minorticks_on()
plot.savefig('Ratio_flatT_deg_PZ8_Tau=1.pdf', format='pdf')
plot.show()
'''

#Ratio between flatA and deg Tau=0.833
plot.plot(freqs, flatA_56_1/(8*deg_56_1), color='black', label=r'$\nu = 1$')
plot.plot(freqs, flatA_56_3/(8*deg_56_3), color='red', label=r'$\nu = 3$')
plot.plot(freqs, flatA_56_5/(8*deg_56_5), color='green', label=r'$\nu = 5$')
plot.plot(freqs, flatA_56_7/(8*deg_56_7), color='orange', label=r'$\nu = 7$')
plot.plot(freqs, flatA_56_9/(8*deg_56_9), color='purple', label=r'$\nu = 9$')
plot.legend(loc = 'lower right')
plot.xscale('log')
plot.ylim(0.9, 1.1)
plot.grid(True, which="both", linestyle=':', linewidth=0.75)
plot.minorticks_on()
plot.savefig('Ratio_flatA_deg_PZ8_Tau=0.833.pdf', format='pdf')
plot.show()

#Ratio between flatT and deg Tau=0.833
plot.plot(freqs, flatT_56_1/deg_56_1, color='black', label=r'$\nu = 1$')
plot.plot(freqs, flatT_56_3/deg_56_3, color='red', label=r'$\nu = 3$')
plot.plot(freqs, flatT_56_5/deg_56_5, color='green', label=r'$\nu = 5$')
plot.plot(freqs, flatT_56_7/deg_56_7, color='orange', label=r'$\nu = 7$')
plot.plot(freqs, flatT_56_9/deg_56_9, color='purple', label=r'$\nu = 9$')
plot.legend(loc = 'lower right')
plot.xscale('log')
plot.ylim(0.9,1.1)
plot.grid(True, which="both", linestyle=':', linewidth=0.75)
plot.minorticks_on()
plot.savefig('Ratio_flatT_deg_PZ8_Tau=0.833.pdf', format='pdf')
plot.show()


'''
plot.plot(freqs40, deg_rot_o1_1e3, color='black')
plot.plot(freqs40, deg_rot_o1_1e4, color='red')

plot.xscale('log')
plot.yscale('log')
#plot.xlabel('Frequenz (Hz)')
#plot.ylabel('Verluste (W/m)')
plot.grid(True, which='both')
plot.savefig('Losses_tau=1_deg_rotor_PZ8_o1_2grids.pdf', format='pdf')
plot.show()
'''
