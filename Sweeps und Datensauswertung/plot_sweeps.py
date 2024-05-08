import matplotlib.pyplot as plot
import numpy as np

time = 0

round1 = "sweep_round_onlymag_1_1_PZ2.csv"
round3 = "sweep_round_onlymag_1_3_PZ2.csv"
round5 = "sweep_round_onlymag_1_5_PZ2.csv"
round7 = "sweep_round_onlymag_1_7_PZ2.csv"
round9 = "sweep_round_onlymag_1_9_PZ2.csv"
flat1 = "sweep_flat_onlymag_1_1_PZ2.csv"
flat3 = "sweep_flat_onlymag_1_3_PZ2.csv"
flat5 = "sweep_flat_onlymag_1_5_PZ2.csv"
flat7 = "sweep_flat_onlymag_1_7_PZ2.csv"
flat9 = "sweep_flat_onlymag_1_9_PZ2.csv"
deg1 = "sweep_deg_onlymag_1_1_PZ2.csv"
deg3 = "sweep_deg_onlymag_1_3_PZ2.csv"
deg5 = "sweep_deg_onlymag_1_5_PZ2.csv"
deg7 = "sweep_deg_onlymag_1_7_PZ2.csv"
deg9 = "sweep_deg_onlymag_1_9_PZ2.csv"
if time is 1:
    degtime = "sweep_deg_time_1_1.csv"
    roundtime = "sweep_round_time_1_1.csv"
#load data from file
#load frequencies
freqs = np.loadtxt(round1, delimiter=',', usecols=(0,))
#load p_values
deg_o1 = np.loadtxt(deg1, delimiter=',', usecols=(1,))

#load p_flat
flat_o1 = np.loadtxt(flat1, delimiter=',', usecols=(1,))

flat_o3 = np.loadtxt(flat3, delimiter=',', usecols=(1,))

flat_o5 = np.loadtxt(flat5, delimiter=',', usecols=(1,))

flat_o7 = np.loadtxt(flat7, delimiter=',', usecols=(1,))

flat_o9 = np.loadtxt(flat9, delimiter=',', usecols=(1,))

round_o1 = np.loadtxt(round1, delimiter=',', usecols=(1,))

round_o3 = np.loadtxt(round3, delimiter=',', usecols=(1,))

round_o5 = np.loadtxt(round5, delimiter=',', usecols=(1,))

round_o7 = np.loadtxt(round7, delimiter=',', usecols=(1,))

round_o9 = np.loadtxt(round9, delimiter=',', usecols=(1,))

deg_o3 = np.loadtxt(deg3, delimiter=',', usecols=(1,))

deg_o5 = np.loadtxt(deg5, delimiter=',', usecols=(1,))

deg_o7 = np.loadtxt(deg7, delimiter=',', usecols=(1,))

deg_o9 = np.loadtxt(deg9, delimiter=',', usecols=(1,))

if time is 1:
    deg_t = np.loadtxt(degtime, delimiter=',', usecols=(1,))
    deg_se = np.loadtxt(degtime, delimiter=',', usecols=(0,))

    round_t = np.loadtxt(roundtime, delimiter=',', usecols=(1,))
    round_se = np.loadtxt(roundtime, delimiter=',', usecols=(0,))

PZ =2

plot.plot(freqs, deg_o1, color='black')
plot.plot(freqs, deg_o3, color='red')
plot.plot(freqs, deg_o5, color='green')
plot.plot(freqs, deg_o7, color='orange')
plot.plot(freqs, deg_o7, color='purple')
plot.plot(freqs, deg_o9, color='blue')

plot.xscale('log')
plot.yscale('log')
plot.xlabel('Frequenz (Hz)')
plot.ylabel('Verluste (W/m)')
plot.grid(True, which='both')
plot.savefig('Losses_tau=0.75_deg_1Magnet.pdf', format='pdf')
plot.show()

plot.plot(freqs, flat_o1, color='black')
plot.plot(freqs, flat_o3, color='red')
plot.plot(freqs, flat_o5, color='green')
plot.plot(freqs, flat_o7, color='orange')
plot.plot(freqs, flat_o7, color='purple')
plot.plot(freqs, flat_o9, color='blue')
plot.xscale('log')
plot.yscale('log')
plot.xlabel('Frequenz (Hz)')
plot.ylabel('Verluste (W/m)')
plot.grid(True, which='both')
plot.savefig('Losses_tau=0.75_flat_1Magnet.pdf', format='pdf')
plot.show()

plot.plot(freqs, round_o1, color='black')
plot.plot(freqs, round_o3, color='red')
plot.plot(freqs, round_o5, color='green')
plot.plot(freqs, round_o7, color='orange')
plot.plot(freqs, round_o7, color='purple')
plot.plot(freqs, round_o9, color='blue')
plot.xscale('log')
plot.yscale('log')
plot.xlabel('Frequenz (Hz)')
plot.ylabel('Verluste (W/m)')
plot.grid(True, which='both')
plot.savefig('Losses_tau=0.75_round_1Magnet_onlymag.pdf', format='pdf')
plot.show()
'''
plot.plot(freqs, p_round)
plot.xscale('log')
plot.yscale('log')
plot.xlabel('Frequency (Hz)')
plot.ylabel('Losses (W/m)')
plot.grid(True, which='onlymag')
plot.savefig('Losses_tau=1_round_1Magnet.pdf', format='pdf')
plot.show()
'''

r_vs_f=[]
for i in np.arange(0,len(freqs)):
    d = deg_o1[i]/flat_o1[i]
    r_vs_f.append(d)

plot.semilogx(freqs, r_vs_f)
plot.xlabel('Frequenz (Hz)')
plot.ylabel('Verhältnis Modell 1/Modell 2')
plot.grid(True, which='both')
plot.savefig('Ratio_half_vs_flat_tau=1.pdf', format='pdf')
plot.show()

r_vs_f=[]
for i in np.arange(0,len(freqs)):
    d = deg_o1[i]/round_o1[i]
    r_vs_f.append(d)

plot.semilogx(freqs, r_vs_f)
plot.xlabel('Frequenz (Hz)')
plot.ylabel('Verhältnis Modell 2/Modell 3')
plot.grid(True, which='both')
plot.savefig('Ratio_half_vs_round_tau=1.pdf', format='pdf')
plot.show()

if time is 1:
    time=[]
    for i in np.arange(0,len(freqs)):
        d = round_t[i]/deg_t[i]
        time.append(d)

    plot.plot(round_se, time)
    plot.xlabel('Number of surface elements round geometry')
    plot.ylabel('Time from matrix assembly till solution')
    plot.grid(True, which='both')
    plot.savefig('Time_round_vs_deg.pdf', format='pdf')
    plot.show()



