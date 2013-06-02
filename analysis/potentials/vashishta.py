# -*- coding: utf-8 -*-
"""
Created on Mon May 27 18:25:11 2013

@author: svenni
"""

from pylab import *

r = linspace(0.05, 7, 1000)

r4s = 2.6

limit = 70

Z = {"Si": 8.0, "O": -0.8}
a = {"Si": 0.0, "O": 2.40}

eta = {"SiSi": 11, "SiO": 9, "OO": 7}
H = {"SiSi": 0.057, "SiO": 11.387, "OO": 51.692}

def ewald(qi, qj, rij):
    alpha = 1 / 4.43
    result = zeros(len(rij))
    for i in range(len(rij)):
        result[i] = 0.5 * qi * qj * math.erfc(alpha * rij[i]) / rij[i]
    return result

figure("SiSi")
VH = H["SiSi"] / r**eta["SiSi"]
Vcol = Z["Si"] * Z["Si"] / r
Vewald = ewald(Z["Si"], Z["Si"], r)
Vexp = -(a["Si"] * Z["Si"]**2) * exp(-r / r4s) / r**4
Vtot = VH + Vcol + Vexp
plot(r, VH, label="VH")
plot(r, Vcol, label="Vcol")
plot(r, Vexp, label="Vexp")
plot(r, Vtot, label="Vtot")
ylim(-limit,limit)
grid()
legend()

figure("SiSiEwald")
Vewald = ewald(Z["Si"], Z["Si"], r)
plot(r, Vcol, label="Vcol")
plot(r, Vewald, label="Vewald")
ylim(-limit,limit)
grid()
legend()

figure("SiSiEwald2")
Vewald = ewald(Z["Si"], Z["Si"], r)
Vtot = VH + Vewald + Vexp
plot(r, Vewald, label="Vewald")
plot(r, Vexp, label="Vexp")
plot(r, Vtot, label="Vtot")
ylim(-limit,limit)
grid()
legend()

figure("SiO")
VH = H["SiO"] / r**eta["SiO"]
Vcol = Z["Si"] * Z["O"] / r
Vewald = ewald(Z["Si"], Z["O"], r)
Vexp = -0.5*(a["Si"] * Z["O"]**2 + a["O"] * Z["Si"]**2) * exp(-r / r4s) / r**4
Vtot = VH + Vcol + Vexp
plot(r, VH, label="VH")
plot(r, Vcol, label="Vcol")
plot(r, Vexp, label="Vexp")
plot(r, Vtot, label="Vtot")
ylim(-limit,limit)
grid()
legend()

figure("SiOEwald")
Vewald = ewald(Z["Si"], Z["O"], r)
Vtot = VH + Vewald + Vexp
plot(r, VH, label="VH")
plot(r, Vewald, label="Vewald")
plot(r, Vexp, label="Vexp")
plot(r, Vtot, label="Vtot")
ylim(-limit,limit)
grid()
legend()

figure("OO")
VH = H["OO"] / r**eta["OO"]
Vcol = Z["O"] * Z["O"] / r
Vewald = ewald(Z["O"], Z["O"], r)
Vexp = -(a["O"] * Z["O"]**2) * exp(-r / r4s) / r**4
Vtot = VH + Vcol + Vexp
plot(r, VH, label="VH")
plot(r, Vcol, label="Vcol")
plot(r, Vexp, label="Vexp")
plot(r, Vtot, label="Vtot")
ylim(-1,1)
legend()
grid()


figure("OOEwald")
Vewald = ewald(Z["O"], Z["O"], r)
Vtot = VH + Vewald + Vexp
plot(r, VH, label="VH")
plot(r, Vewald, label="Vewald")
plot(r, Vexp, label="Vexp")
plot(r, Vtot, label="Vtot")
ylim(-limit,limit)
grid()
legend()

show()
