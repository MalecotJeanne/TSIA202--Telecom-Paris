# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 10:32:40 2022

@author: fontaine
"""

import numpy as np
from scipy import signal as sig
import matplotlib.pyplot as plt
import matplotlib.widgets as wid
import randproc as rp

# %% Initializations
Nit = 500     # number of paths
H = 20      # number of points for covariance
n = 2*H - 1  # number of process' samples
displayRatio = 4  # 1/displayRatio graphs are displayed
# increase displayRatio to speed up,
# in particular if Nit/displayRatio>100

EX = np.zeros(n)          # Expectations are stored here
gamma = np.zeros(n)          # ACF's are stored here
tc = np.arange(-(H-1), H)  # temporal axis for ACF
t = np.arange(0, n)       # temporal axis for the signal


# %% Create a figure for interactive display of results
class btnAction(object):
    stopnow = 0

    def stop(self, event):
        self.stopnow = 1
        if k == Nit:
            plt.close(1)

# Choose different signals here
# signal = 'SIN' # 'AR' 'WN' 'SIN'


# Create the figure where the results will be shown
plt.figure(1, [H, 8])
plt.subplot(211, label='mean')
plt.grid()
plt.title('Sample averages')
plt.subplot(212, label='ACF')
plt.grid()
plt.title('Sample covariance')
axstop = plt.axes([0.81, 0.01, 0.08, 0.05])
btn = wid.Button(axstop, 'Stop')
callback = btnAction()
btn.on_clicked(callback.stop)


# %% Iterations
for k in range(Nit+1):
    # k-th realization of the process

    if signal == 'AR':

        # Generate an AR(1) signal
        Z = np.random.normal(0, 1, n)  # White noise
        # Build the rational function P/Q
        phi1 = 0.6
        Pcoeffs = np.array([1.])
        Qcoeffs = np.poly((phi1,))
        # Use P/Q for filtering
        X = sig.lfilter(Pcoeffs, Qcoeffs, Z)

    elif signal == 'SIN':
        #  Armonic signal
        # A=2.;
        #A_0= -A+2*A*np.random.random_sample();
        A_0 = 1.
        omega = np.pi/3.
        phi = np.pi + 2*np.pi*np.random.random_sample()
        X = A_0 * np.cos(omega*tc+phi) + np.random.normal(0, 1, n)

    else:  # white noise
        sigma = 1.
        X = np.random.normal(0, sigma, n)

    # Compute empirical means
    EX = EX + X

    # Compute exmpirical covariances
    for h in range(H):
        # compiute cov for lag = h
        XX = X[0:n-H] * X[h:h+n-H]
        EXX = XX.sum() / float(n-H)      # computation of the unbiased ACF estim
        gamma[h+H-1] = gamma[h+H-1] + EXX  # gamma[H-1] = cov(X_n, X_n)
        gamma[H-h-1] = gamma[h+H-1].conj()

    # Displaying results (only once every "displayRatio" iterations)
    if np.mod(k, displayRatio) == 0:
        plt.figure(1)

        # Dislpay sample averages
        plt.subplot(211, label='mean')
        if k == 0:  # Create the line at iteration 0
            line1, = plt.plot(t, EX, '-o')
        else:  # update the average EX/(k+1)
            line1.set_ydata(EX/(k+1))
            plt.ylim([-1.5, 1.5])
            plt.xticks(t)

        # Display the unbiased sample covariance
        plt.subplot(212, label='ACF')
        if k == 0:  # Create the line at iteration 0
            line2, = plt.plot(tc, gamma, '-o')
        else:  # Update the sample covariance plot
            line2.set_ydata(gamma/(k+1))

        plt.pause(0.001)
        plt.ylim([-1.2, 2])
        plt.show()
        plt.xticks(tc)
        plt.text(4, -1, 'Iteration %5d' %
                 k, bbox=dict(facecolor='white', alpha=1))

    if callback.stopnow == 1:
        plt.close(1)
        break


# %%Exercice 1


A0 = 1
w = np.pi/3
sigma = 1


# SIN
x1 = [A0**2/2*np.cos(w*k)for k in range(-19, 20)]
x1[len(x1)//2] += sigma*sigma

# AR
x2 = [0.6**abs(k)*sigma*sigma/(1-0.6*0.6) for k in range(-19, 20)]

# WN
x3 = [0 for k in range(-19, 20)]
x3[(len(x3))//2] = 1

# Echelle des temps
t = [k for k in range(-19, 20)]

plt.clf()
plt.figure(2, [H, 8])
plt.subplot(311, label='SIN')
plt.grid()
plt.plot(t, x1)
plt.title('sinusoidal signal')
plt.subplot(312, label='AR')
plt.grid()
plt.plot(t, x2)
plt.title('AR signal')
plt.subplot(313, label='WN')
plt.grid()
plt.plot(t, x3)
plt.title('white noise')
plt.show()

# p.drawZ_DTFT_AR(rp.genAR(4,n))

plt.close('all')
n = 1000
p = 4

X, phi = rp.genAR(p, n)

rp.drawZ_DTFT_AR(X, phi)
# %%Exercice 2


# 1)

# 2)

# We define DTFT
def DTFT(y, lambd, m):
    res = []
    for i in range(len(y)):
        res.append(y[i]*np.exp(-1j*lambd*i))
    res = sum(res)/(2*np.pi)
    return res

# Now we compute In


def I(y, n, m, k):
    return 2*np.pi/n*abs(DTFT(y, 2*np.pi*k/m, m))**2


X1, X2, X3 = [], [], []
Nit = 100
n = 39


# AR
Z = np.random.normal(0, 1, n)
phi1 = 0.6
Pcoeffs = np.array([1.])
Qcoeffs = np.poly((phi1,))
X2 = sig.lfilter(Pcoeffs, Qcoeffs, Z)

# SIN
A_0 = 1.
omega = np.pi/3.
phi = np.pi + 2*np.pi*np.random.random_sample()
X1 = A_0 * np.cos(omega*tc+phi) + np.random.normal(0, 1, n)

# WN
sigma = 1.
X3 = np.random.normal(0, sigma, n)
WN = np.random.normal(0, sigma, 100)  # signal de white noise plus long


# We compute In for each process
I1 = [I(X1, 20, 39, k) for k in range(39)][1:]
I2 = [I(X2, 20, 39, k) for k in range(39)][1:]
I3 = [I(X3, 20, 39, k) for k in range(39)][1:]

tbis = [k+0.5 for k in range(-len(I1)//2, len(I1)//2)]


plt.clf()
plt.figure(3, [H, 8])
plt.subplot(311, label='AR')
plt.grid()
plt.title('AR')
plt.plot(tbis, I2)
plt.subplot(312, label='WN')
plt.grid()
plt.title('WN')
plt.plot(tbis, I3)
plt.subplot(313, label='SIN')
plt.grid()
plt.title('SIN')
plt.plot(tbis, I1)
plt.show()


V = np.zeros(6)
n = [10, 20, 50, 100, 500, 1000]
for i in range(6):
    V[i] = [I(X3, n[i], 39, k) for k in range(39)][0]

plt.figure()
plt.plot(n, V)
plt.xlabel("n")
plt.ylabel("var")
plt.show()
