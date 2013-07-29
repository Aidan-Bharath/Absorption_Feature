# written by Aidan Bharath
# Fabry-Perot Etalon Absorption Features


import numpy as np
from matplotlib import pyplot as plt

# Variable Definitions

dg = 1       # HWHM of the incident light source
dh = 1         # HWHM of flatness defect  Gaussian
df =8         # Free spectral range of the Etalon
dt = 4       # HWHM of the top hat function

sigma = np.arange(-10,10,0.1)
sigma_o = np.arange(-10,10,0.1)
# range of wavenumbers
# np.arange(start,stop,stepsize)



k = .01	       # absorption dip lorentzian HWHM
l = 0.1	       # convolution lorentzin HWHM
A = 0          # Etalon absorption coefficient
R = 0.99       # Etalon reflection Coefficient
hr = 1         # decimal value of absorption
a = 50     # temporary fill in soon

###### Code Functionality
# 1 = true , 0 = false

abs_dip = 0
etalon_trans =1
gauss_conv = 0
lorentz_conv = 0
tophat_conv = 0

####### Definitions ############

#Gaussian

def G(dg,df):
    return (2*np.pi*dg)/(df*np.sqrt(np.log(2)))

def S(df,sigma,sigma_o):
    return ((2*np.pi)/df)*(sigma-sigma_o)

def Gauss(dg,df,sigma,sigma_o):
    return (1/(np.sqrt(np.pi)*(G(dg,df))))*np.exp(-(S(df,sigma,sigma_o)/G(dg,df))**2)

#Lorentzian

def L(l,df):
    return (2*np.pi*l)/df

def T(df,sigma,sigma_o):
    return ((2*np.pi)/df)*(sigma-sigma_o)

def Lorentz(l,df,sigma,sigma_o):
    return (L(l,df)/(np.pi))*(1/(L(l,df)**2+T(df,sigma,sigma_o)**2))


#TopHat Function

def TH(df,sigma,sigma_o):
    return ((2*np.pi)/df)*(sigma-sigma_o)

def tophat(dt,df,sigma,sigma_o):

    #toha = np.zeros(len(sigma))
    #toha[abs(TH(df,sigma,sigma_o))<dt] = 1

    if abs(TH(df,sigma,sigma_o))<dt:
        return 1
    else:
        return 0


#Etalon Transmission Function

def S(df,sigma,sigma_o):
    return ((2*np.pi)/df)*(sigma-sigma_o)

def Etalon(R,A,df,sigma,sigma_o):
    return ((1-(A/(1-R)))**2)*((1-R**2)/(2*np.pi))*(1+R**2-2*R*np.cos(S(df,sigma,sigma_o)))**(-1)


#Incident Lineshape

def incident(dg,df,sigma,sigma_o,L=None):
    if abs_dip == 1:
        return Gauss(dg,df,sigma,sigma_o)*np.exp(-a*Lorentz(L,df,sigma,sigma_o))

    else:
        return Gauss(dg,df,sigma,sigma_o)


#Gaussian Convolution

def conv_gauss(sigma,sigma_o):
    total = 0
    for i in xrange(len(sigma_o)):
        total += Gauss(dh,df,sigma_o[i],0)*incident(dg,df,sigma,sigma_o[i],l)
    return total


#Lorentzian Convolution

def conv_lorentz(sigma,sigma_o):
    if gauss_conv == 1:
        total = 0
        for i in xrange(len(sigma_o)):
            total += lineshape1[i]*Lorentz(k,df,sigma,sigma_o[i])
        return total
    else:
        total = 0
        for i in xrange(len(sigma_o)):
            total += Lorentz(k,df,sigma_o[i],0)*incident(dg,df,sigma,sigma_o[i],l)
        return total


#Tophat Convolution

def conv_tophat(sigma,sigma_o):
    if gauss_conv == 1 and lorentz_conv == 1:
        total = 0
        for i in xrange(len(sigma_o)):
            total += lineshape2[i]*tophat(dt,df,sigma_o[i],sigma)
        return total

    elif gauss_conv == 1 and lorentz_conv == 0:
        total = 0
        for i in xrange(len(sigma_o)):
            total += lineshape1[i]*tophat(dt,df,sigma_o[i],sigma)
        return total


    elif gauss_conv == 0 and lorentz_conv == 1:
        total = 0
        for i in xrange(len(sigma_o)):
            total += lineshape2[i]*tophat(dt,df,sigma_o[i],sigma)
        return total

    else:
        total = 0
        for i in xrange(len(sigma_o)):
            total += tophat(dt,df,sigma_o[i],0)*incident(dg,df,sigma,sigma_o[i],l)
        return total

# Etalon Transmission Function Convolution

def conv_etalon(sigma,sigma_o):
    if gauss_conv == 1 and lorentz_conv == 1 and tophat_conv == 1:
        total = 0
        for i in xrange(len(sigma_o)):
            total += lineshape3[i]*Etalon(R,A,df,sigma_o[i],sigma)
        return total

    elif gauss_conv == 1 and lorentz_conv == 1 and tophat_conv == 0:
        total = 0
        for i in xrange(len(sigma_o)):
            total += lineshape2[i]*Etalon(R,A,df,sigma_o[i],sigma)
        return total

    elif gauss_conv == 1 and lorentz_conv == 0 and tophat_conv == 1:
        total = 0
        for i in xrange(len(sigma_o)):
            total += lineshape3[i]*Etalon(R,A,df,sigma_o[i],sigma)
        return total

    elif gauss_conv == 0 and lorentz_conv == 1 and tophat_conv == 1:
        total = 0
        for i in xrange(len(sigma_o)):
            total += lineshape3[i]*Etalon(R,A,df,sigma_o[i],sigma)
        return total

    elif gauss_conv == 1 and lorentz_conv == 0 and tophat_conv == 0:
        total = 0
        for i in xrange(len(sigma_o)):
            total += lineshape1[i]*Etalon(R,A,df,sigma,sigma_o[i])
        return total

    elif gauss_conv == 0 and lorentz_conv == 1 and tophat_conv == 0:
        total = 0
        for i in xrange(len(sigma_o)):
            total += lineshape2[i]*Etalon(R,A,df,sigma_o[i],sigma)
        return total

    elif gauss_conv ==0 and lorentz_conv == 0 and tophat_conv == 1:
        total = 0
        for i in xrange(len(sigma_o)):
            total += lineshape3[i]*Etalon(R,A,df,sigma_o[i],sigma)
        return total

    else:
        total = 0
        for i in xrange(len(sigma_o)):
            total += Etalon(R,A,df,sigma_o[i],0)*incident(dg,df,sigma,sigma_o[i],l)
        return total



####### Calculations ####################



# Gaussian Linshape Convolution

if gauss_conv == 1:

    lineshape1 = np.zeros(len(sigma))

    for i in xrange(len(sigma)):
        lineshape1[i] = conv_gauss(sigma[i],sigma_o)


# Lorentzian Lineshape Convolution

if lorentz_conv == 1:

    lineshape2 = np.zeros(len(sigma))

    for i in xrange(len(sigma)):
        lineshape2[i] = conv_lorentz(sigma[i],sigma_o)


# Top Hat Lineshape Convolution

if tophat_conv == 1:

    lineshape3 = np.zeros(len(sigma))

    for i in xrange(len(sigma)):
        lineshape3[i] = conv_tophat(sigma[i],sigma_o)


# Convolution with Ideal Etalon

if etalon_trans == 1:

    lineshape4 = np.zeros(len(sigma))

    for i in xrange(len(sigma)):
        lineshape4[i] = conv_etalon(sigma[i],sigma_o)


### Plotting functions
hat = np.zeros(len(sigma))
for i in xrange(len(sigma)):
    hat[i] = Etalon(R,A,df,sigma[i],sigma_o[10])
plt.plot(sigma,lineshape4)
plt.show()




