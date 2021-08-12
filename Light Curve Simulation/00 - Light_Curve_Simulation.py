"""
Started in: July 28, 2021
Finished in: 

@author: Guilherme Samuel
"""

from os import system
system("cls")

# Libraries
from math import acos, log, pi, sqrt
import numpy as np
from scipy import interpolate
from numpy import loadtxt


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
def calculate_lambda_1(a, b, k, p, q, w, z):
    return (1/(9*pi*sqrt(p*z[w]))) * (((1-b)*(2*b+a-3)-3*q*(b-2))*ellk(k)+4*p*z*(z[w]**2+7*p**2-4)*ellec(k)-(3*q/a)*ellpic_bulirsch(abs((a-1)/a), k))

def calculate_lambda_2(a, b, k, p, q, w, z):
    return (2/(9*pi*sqrt(1-a))) * ((1-5*z[w]**2+p**2 + q**2)*ellk(k**(-1))+(1-a)*(z[w]**2+7*p**2-4)*ellec(k**(-1))-(3*q/a)*ellpic_bulirsch(abs((a-b)/a), k**(-1)))

def calculate_lambda_3():
    pass

def calculate_lambda_4():
    pass

def calculate_lambda_5():
    pass

def calculate_lambda_6():
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

# Complete path to the lightcurve phase folded
observed_curve = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
# observed_curve_path = 'files\curva_luz_eclipse_medio_ID100725706_Butterworth_n2_f02_autocalibrada.txt'
# observed_curve = loadtxt(observed_curve_path, dtype='float', delimiter='\n')


# Sampling intervals of parameter's values
delta_b = 0.01
delta_adivR = 0.1
delta_periodo = 0.01
delta_p = 0.01


# Planet coordinate, along the x-axis, as a function of the start's radius
x_values_path = 'files\Valores_b_simulacao.txt'
x_values = loadtxt(x_values_path, dtype='float', delimiter='\n')

# Transit impact parameter
b_values_path = 'files\Valores_b_simulacao.txt'
b_values = loadtxt(b_values_path, dtype='float', delimiter='\n')

# Radius values of the planet compared to the star
p_values_path = 'files\Valores_p_simulacao.txt'
p_values = loadtxt(p_values_path, dtype='float', delimiter='\n')

# Orbital period values to be considered
period_values_path = 'files\Valores_periodo_simulacao.txt'
period_values = loadtxt(period_values_path, dtype='float', delimiter='\n')

# Orbital radius values compared to star radius
adivR_values_path = 'files\Valores_adivR_simulacao.txt'
adivR_values = loadtxt(adivR_values_path, dtype='float', delimiter='\n')


# Coefficientes of limb darkening
gamma1 = 0.44
gamma2 = 0.23

# Lightcurves simulations (with limb darkening) and calculation of chi squares


simulated_curve = np.ones((2, len(x_values)))
simulated_curve_resampled = np.ones((2, len(observed_curve)))
observed_curve_eclipse = np.ones((3, len(observed_curve)))

F = []
z = []

v = 0
for b_impact in b_values:
    z = np.sqrt((np.power(x_values, 2) + b_impact))

    for p in p_values:

        for period in period_values:

            for adivR in adivR_values:

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
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    elif (p == 0) and (z[w] >= 0):
                        # print('Case I')
                        lambda_d = 0
                        eta_d = 0    
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    ## Case II
                    elif (p > 0) and (z[w] > 0.5+abs(p-0.5)) and (z[w] < 1+p):
                        # print('Case II')
                        lambda_d = calculate_lambda_1(a, b, k, p, q, w, z)
                        eta_d = calculate_eta_1(a, b, k0, k1, p, w, z)

                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    ## Case III
                    elif (p > 0) and (p < 0.5) and (z[w] > p) and (z[w] < 1-p):
                        # print('Case III')
                        lambda_d = calculate_lambda_2(a, b, k, p, q, w, z)
                        eta_d = calculate_eta_2(p, w, z) 

                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    ## Case IV
                    elif (p > 0) and (p < 0.5) and (z[w] == 1-p):
                        # print('Case IV')
                        lambda_d = calculate_lambda_5() #TODO
                        eta_d = calculate_eta_2(p, w, z)
                        
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
 
                    ## Case V
                    elif (p > 0) and (p < 0.5) and (z[w] == p):
                        # print('Case V')
                        lamda_d = calculate_lambda_4() #TODO
                        eta_d = calculate_eta_2(p, w, z)

                        F[w] = lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z)

                    ## Case VI
                    elif (p == 0.5) and (z[w] == 0.5):
                        # print('Case VI')
                        lambda_d = (1/3) - (4/(9*pi))
                        eta_d = 3/32
                        
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    ## Case VII
                    elif (p > 0.5) and (z[w] == p):
                        # print('Case VII')
                        lambda_d = calculate_lambda_3() #TODO
                        eta_d = calculate_eta_1(a, b, k0, k1, p, w, z)

                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    ## Case VIII
                    elif (p > 0.5) and (z[w] >= abs(1-p)) and (z[w] < p):
                        # print('Case VIII')
                        lambda_d = calculate_lambda_1(a, b, k, p, q, w, z)
                        eta_d = calculate_eta_1(a, b, k0, k1, p, w, z)
                        
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    ## Case IX
                    elif (p > 0) and (p < 1) and (z[w] > 0) and (z[w] < 0.5-abs(p-0.5)):
                        # print('Case IX')
                        lambda_d = calculate_lambda_2(a, b, k, p, q, w, z)
                        eta_d = calculate_eta_2(p, w, z) 

                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    ## Case X
                    elif (p > 0) and (p < 1) and (z[w] == 0):
                        # print('Case X')
                        lambda_d = calculate_lambda_6() #TODO
                        eta_d = calculate_eta_2(p, w, z)
                        
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

                    ## Case XI
                    elif (p > 1) and (z[w] >= 0) and (z[w] < p-1):
                        print('Case XI')
                        lambda_d = 1
                        eta_d = 1
                        
                        F.insert(w, lightcurve_F(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
            
#############################################################################################            
                # Conversao de x para tempo
                ttrans = period/(pi*adivR)
                mid = (len(x_values)/2) - 0.5
                for u in range(len(x_values)):
                    if (u < mid):
                        simulated_curve[0, u] = (-1)*x_values[u]*ttrans/2
                    elif (u >= mid):
                        simulated_curve[0, u] = x_values[u]*ttrans/2
                
                # simulated_curve[1, *] = F[*]          
                    
                # Reamostragem
                # simulated_curve_resampled[0, *] = observed_curve_eclipse[0, *]
                # simulated_curve_resampled[1, *] = interpolate.interp1d(simulated_curve[1, *], simulated_curve[0, *], simulated_curve_resampled[0, *])
                
                # Cálculo do qui^2
                qui2 = 0
                for u in range(len(observed_curve)):
                    qui2 += ((observed_curve_eclipse[1, u]-simulated_curve_resampled[1, u])**2/observed_curve_eclipse[2, u]**2)

                final_table = np.ones(5)
                final_table[0] = b_impact
                final_table[1] = p
                final_table[2] = period
                final_table[3] = adivR
                final_table[4] = qui2

                # string_v = str(v+1)
                # string_v = strtrim(string_v, 1)
                # v += 1

                # print("b_impact:", final_table[0])
                # print("p:", final_table[1])
                # print("period:", final_table[2])
                # print("adivR:", final_table[3])
                # print("qui2:", final_table[4])
    break
# 481 line

# Ordenamento da tabela de quis quadrados
print(qui2)




