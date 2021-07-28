"""
Created on ...

@author: Guilherme Samuel

1º: O que é o z ? Como define z ? 
2º: O que é o que em cada for (p, period, adivR, w): range(len(array)) or for i in array

"""
# TEM QUE DAR UM CONTROL -, TIRAR ZOOM
# Libraries
from math import acos, log, pi, sqrt
import numpy as np

# Help functions
def gaussian_sum(X, A):
    return


def ellk(k):
    """
    Computes polynomial approximation for the complete
    elliptic integral of the first kind (Hasting's approximation)

    https://doi.org/10.2307/2004103
    Table II, coefficients values from n=4

    :param FLOAT? k: 

    :return: The complete elliptical integral of the first kind
    :rtype: FLOAT?
    """
    m1 = 1 - k**2

    # Coefficients for K*
    a0 = log(4)
    a1 = 0.09666344259
    a2 = 0.03590092383
    a3 = 0.03742563713
    a4 = 0.01451196212
    b0 = 0.5
    b1 = 0.12498593597
    b2 = 0.06880248576
    b3 = 0.03328355346
    b4 = 0.00441787012

    ek1 = a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
    ek2 = (b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*log(m1)

    return ek1 - ek2


def ellec(k):
    """
    Computes polynomial approximation for the complete
    elliptic integral of the second kind (Hasting's approximation)

    https://doi.org/10.2307/2004103
    Table III, coefficients values from n=4

    :param float k:

    :return: The complete elliptical integral of the second kind
    :rtype: 
    """
    m1 = 1 - k**2

    # Coefficients for E*
    c1 = 0.44325141463
    c2 = 0.06260601220
    c3 = 0.04757383546
    c4 = 0.01736506451
    d1 = 0.24998368310
    d2 = 0.09200180037
    d3 = 0.04069697526
    d4 = 0.00526449639

    ee1 = 1+m1*(c1+m1*(c2+m1*(c3+m1*c4)))
    ee2 = m1*(d1+m1*(d2+m1*(d3+m1*d4)))*log(1/m1)

    return ee1 + ee2


def ellpic_bulirsch(n, k):
    """
    Computes the complete elliptical integral of the third kind 
    using the algorithm of Bulirsch (1965)

    https://doi.org/10.1007/BF02165405

    :param FLOAT? n: 
    :param FLOAT? k: 

    :return: The complete elliptical integral of the third kind
    :rtype: 
    """
    kc = sqrt(1 - k**2)
    p = n + 1

    # if min(p) < 0:
    if p < 0:
        print('Negative p')

    m0 = 1
    c = 1
    p = sqrt(p)
    d = 1/p
    e = kc
    d = 1/p
    e = kc

    iter = 0
    while iter < 20:
        f = c
        c = d/p+c
        g = e/p
        d = 2*(f*g + d)

        p = g + p
        g = m0
        m0 = kc + m0

        # if max(abs(1 - kc/g)) > 1-8:
        if (abs(1 - kc/g)) > 1-8:
            kc = 2 * sqrt(e)
            e = kc * m0

        else:
            return 0.5*pi*(c*m0+d)/(m0*(m0+p))

        iter += 1

    return 0.5*pi*(c*m0+d)/(m0*(m0+p))

## Functions from Table 1, Mandel & Agol (2008)
def lambda_1(a, b, k, p, q, w, z):
    return (1/(9*pi*sqrt(p*z[w]))) * (((1-b)*(2*b+a-3)-3*q*(b-2))*ellk(k)+4*p*z*(z[w]**2+7*p**2-4)*ellec(k)-(3*q/a)*ellpic_bulirsch(abs((a-1)/a), k))

def lambda_2():
    pass

def lambda_3():
    pass

def lambda_4():
    pass

def lambda_5():
    pass

def lambda_6():
    pass

def calculate_eta_2(p, w, z):
    return ((p**2)/2)*(p**2+2*z[w]**2)

def calculate_eta_1(a, b, k0, k1, p, w, z):
    return (2*pi)**(-1)*(k1+2*calculate_eta_2(p, w, z)*k0-0.25*(1+5*p**2+z[w]**2)*sqrt((1-a)*(b-1)))

def lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z):
    if (p <= z[w]):
        return 1-(4*Omega)**(-1)*((1-c2)*lambda_e+c2*(lambda_d-c4*eta_d))

    if (p > z[w]):
        return 1-(4*Omega)**(-1)*((1-c2)*lambda_e+c2*(lambda_d+(2/3)-c4*eta_d))



"""
Simula_Curvas_Luz_Plus
"""

# Sampling intervals of parameter's values
delta_b = 0.01
delta_adivR = 0.1
delta_periodo = 0.01
delta_p = 0.01

# Coefficientes of limb darkening
gamma1 = 0.44
gamma2 = 0.23

# Lightcurves simulations (with limb darkening) and calculation of chi squares
v = 0

#######################################################
b_values = [1, 2]
x_values = [1, 2]
p_values = [1, 2]
period_values = [1, 2]
adivR_values = [1, 2]
#######################################################

F = []

for i in b_values:
    b_impact = b_values[i]
    # z = sqrt(x_values**2 + b_impact**2)
    z = np.sqrt(np.power(x_values, 2) + np.power(b_values, 2))

    for j in p_values:
        p = p_values[j]

        # for l in period_values:
        for l in range(len(period_values)):
            period = period_values[l]

            # for m in adivR_values:
            for m in range(len(adivR_values)):
                adivR = adivR_values[m]

                # for w in x_values:
                for w in range(len(x_values)):

                    # Application of basic modeling
                    if ((1+p) < z[w]):
                        lambda_e = 0

                    elif (abs(1-p) < z[w] and (z[w] <= (1+p))):
                        k0 = acos((p**2 + z[w]**2 - 1) / (2*p*z[w]))
                        k1 = acos((1 - p**2 + z[w]**2) / (2*z[w]))
                        lambda_e = (1/pi) * (p**2*k0 + k1 -
                                             sqrt((4*z[w]**2 - (1+z[w]**2-p**2)**2)/(4)))

                    elif (z[w] <= 1-p):
                        lambda_e = p**2

                    elif (z[w] <= p-1):
                        lambda_e = 1

                    # Application of limb darkening
                    a = (z[w]-p)**2
                    b = (z[w]+p)**2
                    k = sqrt((1-a)/(4*z[w]*p))
                    q = p**2 - z[w]**2

                    c1 = 0
                    c2 = gamma1 + 2*gamma2
                    c3 = 0
                    c4 = -1*gamma2
                    c0 = 1 - c1 - c2 - c3 - c4

                    # Omega = sum_{n=0}^{4} c_n/(n+4)
                    Omega = (c0/4) + (c1/5) + (c2/6) + (c3/7) + (c4/8)

                    ## Case I
                    if (p > 0) and (z[w] >= 1+p):
                        print('Case I')
                        lambda_d = 0
                        eta_d = 0
                        # F[w] = lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z)
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    if (p == 0) and (z[w] >= 0):
                        print('Case I')
                        lambda_d = 0
                        eta_d = 0
                        # F[w] = lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z)       
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    ## Case II
                    elif (p > 0) and (z[w] > 0.5+abs(p-0.5)) and (z[w] < 1+p):
                        print('Case II')
                        lamda1 = lambda_1(a, b, k, p, q, w, z)
                        eta_2 = calculate_eta_2(p, w, z)
                        eta_1 = calculate_eta_1(a, b, k0, k1, p, w, z)

                        lambda_d = lamda1
                        eta_d = eta_1
                        # F[w] = lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z)
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    ## Case III
                    elif (p > 0) and (p < 0.5) and (z[w] > p) and (z[w] < 1-p):
                        print('Case III')
                        lambda2 = lambda_2() #TODO
                        eta_2 = calculate_eta_2(p, w, z)

                        lambda_d = lambda2
                        eta_d = eta_2
                        # F[w] = lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z)
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    ## Case IV
                    elif (p > 0) and (p < 0.5) and (z[w] == 1-p):
                        print('Case IV')
                        lambda5 = lambda_5() #TODO
                        eta_2 = calculate_eta_2(p, w, z)

                        lambda_d = lambda5
                        eta_d = eta_2
                        # F[w] = lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z)
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
 
                    ## Case V
                    elif (p > 0) and (p < 0.5) and (z[w] == p):
                        print('Case V')
                        lambda4 = lambda_4() #TODO
                        eta_2 = calculate_eta_2(p, w, z)

                        lamda_d = lambda4
                        eta_d = eta_2
                        F[w] = lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z)

                    ## Case VI
                    elif (p == 0.5) and (z[w] == 0.5):
                        print('Case VI')
                        lambda_d = (1/3) - (4/(9*pi))
                        eta_d = 3/32
                        # F[w] = lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z)
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    ## Case VII
                    elif (p > 0.5) and (z[w] == p):
                        print('Case VII')
                        lambda3 = lambda_3() #TODO
                        eta_2 = calculate_eta_2(p, w, z)
                        eta_1 = calculate_eta_1(a, b, k0, k1, p, w, z)

                        lambda_d = lambda3
                        eta_d = eta_1
                        # F[w] = lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z)
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    ## Case VIII
                    elif (p > 0.5) and (z[w] >= abs(1-p)) and (z[w] < p):
                        print('Case VIII')
                        lambda1 = lambda_1(a, b, k, p, q, w, z)
                        eta_2 = calculate_eta_2(p, w, z)
                        eta_1 = calculate_eta_1(a, b, k0, k1, p, w, z)

                        lambda_d = lambda1
                        eta_d = eta_1
                        # F[w] = lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z)
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    ## Case IX
                    elif (p > 0) and (p < 1) and (z[w] > 0) and (z[w] < 0.5-abs(p-0.5)):
                        print('Case IX')
                        lambda2 = lambda_2() #TODO
                        eta_2 = calculate_eta_2(p, w, z)

                        lambda_d = lambda2
                        eta_d = eta_2
                        # F[w] = lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z)
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    ## Case X
                    elif (p > 0) and (p < 1) and (z[w] == 0):
                        print('Case X')
                        lambda6 = lambda_6() #TODO
                        eta_2 = calculate_eta_2(p, w, z)
                        eta_1 = calculate_eta_1(a, b, k0, k1, p, w, z)

                        lambda_d = lambda6
                        eta_d = eta_2
                        # F[w] = lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z)
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    ## Case XI
                    elif (p > 1) and (z[w] >= 0) and (z[w] < p-1):
                        print('Case XI')
                        lambda_d = 1
                        eta_d = 1
                        # F[w] = lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z)
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
            
            
            # Conversao de x para tempo
            
            # Reamostragem

            # Cálculo do qui^2





