# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:06:13 2013

@author: svenni
"""

from pylab import *

lij = linspace(0,10,1000)
r0 = 2.6
lik = 1
l = 1.0
y = exp(l / (lij - r0) + l / (lik - r0)) * (lij < r0)
figure()
plot(lij,y)

figure()
dydx = diff(y) / diff(lij)
plot(lij[:-1],dydx)

figure()
dfdrij = -l*exp(l/(-r0 + lik) + l/(-r0 + lij))/pow(-r0 + lij, 2) * (lij < r0)
plot(lij, dfdrij)

show()