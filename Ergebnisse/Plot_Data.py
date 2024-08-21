import matplotlib.pyplot as plot
import numpy as np

deg1 = "sweep_deg_onlymag_1_1_PZ8_80samples.csv"
deg3 = "sweep_deg_onlymag_1_3_PZ8_80samples.csv"
deg5 = "sweep_deg_onlymag_1_5_PZ8_80samples.csv"
deg7 = "sweep_deg_onlymag_1_7_PZ8_80samples.csv"
deg9 = "sweep_deg_onlymag_1_9_PZ8_80samples.csv"

flat1 = "sweep_flat_A_onlymag_1_1_PZ8_80samples.csv"
flat3 = "sweep_flat_A_onlymag_1_3_PZ8_80samples.csv"

deg1_rot_1e3 = "sweep_deg_rotor_1_1_PZ8_40samples_1e-3.csv"
deg1_rot_1e4 = "sweep_deg_rotor_1_1_PZ8_40samples_1e-4.csv"

#load frequencies
freqs = np.loadtxt(deg1, delimiter=',', usecols=(0,))
freqs40 = np.loadtxt(deg1_rot_1e3, delimiter=',', usecols=(0,))

deg_o1 = np.loadtxt(deg1, delimiter=',', usecols=(1,))
deg_o3 = np.loadtxt(deg3, delimiter=',', usecols=(1,))
deg_o5 = np.loadtxt(deg5, delimiter=',', usecols=(1,))
deg_o7 = np.loadtxt(deg7, delimiter=',', usecols=(1,))
deg_o9 = np.loadtxt(deg9, delimiter=',', usecols=(1,))

deg_rot_o1_1e3 = np.loadtxt(deg1_rot_1e3, delimiter=',', usecols=(1,))
deg_rot_o1_1e4 = np.loadtxt(deg1_rot_1e4, delimiter=',', usecols=(1,))

flat_o1 = np.loadtxt(flat1, delimiter=',', usecols=(1,))
flat_o3 = np.loadtxt(flat3, delimiter=',', usecols=(1,))

plot.plot(freqs, deg_o1, color='black')
plot.plot(freqs, deg_o3, color='red')
plot.plot(freqs, deg_o5, color='green')
plot.plot(freqs, deg_o7, color='orange')
plot.plot(freqs, deg_o7, color='purple')
plot.plot(freqs, deg_o9, color='blue')

plot.xscale('log')
plot.yscale('log')
#plot.xlabel('Frequenz (Hz)')
#plot.ylabel('Verluste (W/m)')
plot.grid(True, which='both')
plot.savefig('Losses_tau=1_deg_onlymags_PZ8_.pdf', format='pdf')
plot.show()

plot.plot(freqs40, deg_rot_o1_1e3, color='black')
plot.plot(freqs40, deg_rot_o1_1e4, color='red')

plot.xscale('log')
plot.yscale('log')
#plot.xlabel('Frequenz (Hz)')
#plot.ylabel('Verluste (W/m)')
plot.grid(True, which='both')
plot.savefig('Losses_tau=1_deg_rotor_PZ8_o1_2grids.pdf', format='pdf')
plot.show()

