"""

    This module defines a Simulate

"""

# Imports
from utils.models.phasefolded import PhaseFoldedLightCurve
from utils.simulation.help_functions import *

from dataclasses import dataclass, field
import pandas as pd
import numpy as np
import scipy.signal as scs
from tqdm import tqdm

from os import system
system('cls')

@dataclass(frozen=True)
class Simulate():
    """
    Class that applies the modeling created by Mandel & Agol (2008).
    """
    # Constructor
    ## Coefficientes of limb darkening
    gamma1: int = field(default=0.44)
    gamma2: int = field(default=0.23)

    ## Sampling intervals of parameter's values
    delta_b: float = field(default=0.01, repr=False)
    delta_adivR: float = field(default=0.1, repr=False)
    delta_period: float = field(default=0.01, repr=False)
    delta_p: float = field(default=0.01, repr=False)

    ## Best fitting parameters
    b_impact_best: float = field(default=None, repr=False)
    p_best: float = field(default=None, repr=False)
    period_best: float = field(default=None, repr=False)
    adivR_best: float = field(default=None, repr=False)
    chi2_best: float = field(default=None, repr=False)

    ## Final table, computes after the `simulate_values` method run
    simulation_table: pd.DataFrame = field(default=pd.DataFrame(columns=['b_impact', 'p', 'period', 'adivR', 'chi2']), repr=False)

    # Class Methods 
    def simulate_values(self, observed_curve, b_values, p_values, period_values, adivR_values, x_values, set_best_values=True, results_to_csv=False) -> pd.DataFrame:
        """
        Parameters
        ----------
        observed_curve : numpy ndarray
            Lightcurve phase-folded

        b_values : numpy ndarray
            Transit impact parameter

        p_values : numpy ndarray
            Radius values of the planet compared to the star

        period_values : numpy ndarray
            Orbital period values to be considered

        adivR_values : numpy ndarray
            Orbital radius values compared to star radius                

        x_values : numpy ndarray
            Planet coordinate, along the x-axis, as a function of the start's radius

        Returns
        -------
        A `pd.DataFrame`, where the columns are the input parameters of the simulation modeling (b_impact, p, period, adivR) together with the error (Chi-squared) attributed to each combination of these parameters
        """    
        print('Starting simulation...')
        list_b_values = []
        list_p_values = []
        list_period_values = []
        list_adivR_values = []
        list_chi2_values = []

        total = len(b_values) * len(p_values) * len(period_values) * len(adivR_values)

        with tqdm(range(total), colour='blue', desc='Simulating') as pbar:
            for b_impact in b_values:
                for p in p_values:
                    for period in period_values:
                        for adivR in adivR_values:
                            try:
                                simulated_curve = self.simulate(observed_curve=observed_curve, b_impact=b_impact, p=p, period=period, adivR=adivR, x_values=x_values)
                                chi2 = self.calculate_chi2(observed_curve=observed_curve, simulated_curve=simulated_curve)
                            except KeyboardInterrupt:
                                print('Caught KeyboardInterrupt')
                                raise StopIteration

                            except ValueError:
                                print(f"Problem with: \nb_impact = {b_impact}\np = {p}\nperiod = {period}\nadivR = {adivR}")
                                raise 
                                

                            list_b_values.append(b_impact)
                            list_p_values.append(p)
                            list_period_values.append(period)
                            list_adivR_values.append(adivR)
                            list_chi2_values.append(chi2)
                            pbar.update(1)

        
        self.simulation_table['b_impact'] = list_b_values
        self.simulation_table['p'] = list_p_values
        self.simulation_table['period'] = list_period_values
        self.simulation_table['adivR'] = list_adivR_values
        self.simulation_table['chi2'] = list_chi2_values

    
        sorted_table = self.simulation_table.sort_values(by='chi2')
        if set_best_values:
            object.__setattr__(self, 'b_impact_best', sorted_table.loc[sorted_table.index[0]][0])
            object.__setattr__(self, 'p_best', sorted_table.loc[sorted_table.index[0]][1])
            object.__setattr__(self, 'period_best', sorted_table.loc[sorted_table.index[0]][2])
            object.__setattr__(self, 'adivR_best', sorted_table.loc[sorted_table.index[0]][3])
            object.__setattr__(self, 'chi2_best', sorted_table.loc[sorted_table.index[0]][4])
        if results_to_csv:
            sorted_table.to_csv('final_table.csv', index=False)

        print("\nBest parameters are:")
        print('-> Best b_impact =', self.b_impact_best)
        print('-> Best p =', self.p_best)
        print('-> Best period =', self.period_best)
        print('-> Best adivR =', self.adivR_best)
        print('-> Best chi2 =', self.chi2_best, end='\n\n')

        return sorted_table

    def simulate_lightcurve(self, observed_curve=None, b_impact=None, p=None, period=None, adivR=None, x_values=None) -> PhaseFoldedLightCurve:
        """
        Parameters
        ----------
        observed_curve : numpy ndarray
            Lightcurve phase-folded

        b_impact : float
            Transit impact parameter

        p : float
            Radius values of the planet compared to the star

        period : float
            Orbital period values to be considered

        adivR : float
            Orbital radius values compared to star radius                

        x_values : float
            Planet coordinate, along the x-axis, as a function of the start's radius

        Returns
        -------
        A 2-D `np.ndarray`. First column is the time values, second colums is the flux values
        """
        print('\nBuilding the light curve...')
        # If no parameters input was given, the default values are the best ones, computed by the `simulate_values` method
        if b_impact is None:
            print('Using the best b_impact, computed earlier')
            b_impact = self.b_impact_best

        if p is None:
            print('Using the best p, computed earlier')
            p = self.p_best

        if period is None:
            print('Using the best period, computed earlier')
            period = self.period_best

        if adivR is None:
            print('Using the best adivR, computed earlier')
            adivR = self.adivR_best

        # Simulate value
        simulated_curve =  self.simulate(observed_curve, b_impact, p, period, adivR, x_values)
        chi2 = self.calculate_chi2(observed_curve, simulated_curve)

        phasefoldedobject = PhaseFoldedLightCurve(original_time=observed_curve[:, 0], original_flux=observed_curve[:, 1], simulated_time=simulated_curve[:, 0], simulated_flux=simulated_curve[:, 1], chi2=chi2)

        return phasefoldedobject

    def simulate(self, observed_curve, b_impact, p, period, adivR, x_values) -> np.ndarray:
        flux = []

        simulated_curve = np.zeros((len(x_values), 2))
        # First column -> simulated_curve[: 0] => x_values [time]
        # Second column -> simulated_curve[: 1] => flux values

        resampled_simulated_curve = np.zeros((observed_curve.shape[0], 2))
        # First column -> resampled_simulated_curve[: 0] => x_values [time]
        # Second column -> resampled_simulated_curve[: 1] => flux values


        z = np.sqrt((np.power(x_values, 2) + b_impact**2))

        for w in range(len(x_values)):
            # Application of basic modeling
            if ((1+p) < z[w]):
                lambda_e = 0

            elif (abs(1-p) < z[w] and (z[w] <= (1+p))):
                k0 = acos((p**2 + z[w]**2 - 1) / (2*p*z[w]))
                k1 = acos((1 - p**2 + z[w]**2) / (2*z[w]))
                lambda_e = (1/pi) * (p**2*k0 + k1 - sqrt((4*z[w]**2 - (1+z[w]**2-p**2)**2)/(4)))

            elif (z[w] <= 1-p):
                lambda_e = p**2

            elif (z[w] <= p-1):
                lambda_e = 1

            # Application of limb darkening
            a = (z[w]-p)**2
            b = (z[w]+p)**2
            try: 
                k = sqrt((1-a)/(4*z[w]*p))
            except ValueError:
                # print('Math domain error: ')
                # print('Trying to calculate sqrt of', (1-a)/(4*z[w]*p))
                # raise
                pass
            q = (p**2) - (z[w]**2)
            
            c1 = 0
            c2 = self.gamma1
            c3 = 0
            c4 = -1*self.gamma2
            c0 = 1 - c1 - c2 - c3 - c4

            Omega = (c0/4) + (c1/5) + (c2/6) + (c3/7) + (c4/8)

            ## Case I
            # print('Case I') 
            if (p > 0) and (z[w] >= 1+p):
                lambda_d = 0
                eta_d = 0
                
                # print(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
                flux.append(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

            elif (p == 0) and (z[w] >= 0):
                # print('Case I.I')
                lambda_d = 0
                eta_d = 0 

                # print(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
                flux.append(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

            ## Case II
            elif (p > 0) and (z[w] > 0.5+abs(p-0.5)) and (z[w] < 1+p):
                # print('Case II')
                lambda_d = calculate_lambda_1(a, b, k, p, q, w, z)
                eta_d = calculate_eta_1(a, b, k0, k1, p, w, z)
                
                # print(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
                flux.append(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

            ## Case III
            elif (p > 0) and (p < 0.5) and (z[w] > p) and (z[w] < 1-p):
                # print('Case III')
                lambda_d = calculate_lambda_2(a, b, k, p, q, w, z)
                eta_d = calculate_eta_2(p, w, z)

                # print(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
                flux.append(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

            ## Case IV  
            elif (p > 0) and (p < 0.5) and (z[w] == 1-p):
                # print('Case IV')
                lambda_d = calculate_lambda_5(p)
                eta_d = calculate_eta_2(p, w, z)
                
                # print(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
                flux.append(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

            ## Case V
            elif (p > 0) and (p < 0.5) and (z[w] == p):
                # print('Case V')
                lambda_d = calculate_lambda_4(p)
                eta_d = calculate_eta_2(p, w, z)

                # print(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
                flux.append(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

            ## Case VI
            elif (p == 0.5) and (z[w] == 0.5):
                # print('Case VI')
                lambda_d = (1/3) - (4/(9*pi))
                eta_d = 3/32

                # print(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
                flux.append(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

            ## Case VII
            elif (p > 0.5) and (z[w] == p):
                # print('Case VII')
                lambda_d = calculate_lambda_3(p)
                eta_d = calculate_eta_1(a, b, k0, k1, p, w, z)

                # print(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
                flux.append(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

            ## Case VIII
            elif (p > 0.5) and (z[w] >= abs(1-p)) and (z[w] < p):
                # print('Case VIII')
                lambda_d = calculate_lambda_1(a, b, k, p, q, w, z)
                eta_d = calculate_eta_1(a, b, k0, k1, p, w, z)

                # print(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
                flux.append(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

            ## Case IX
            elif (p > 0) and (p < 1) and (z[w] > 0) and (z[w] < 0.5-abs(p-0.5)):
                # print('Case IX')
                lambda_d = calculate_lambda_2(a, b, k, p, q, w, z)
                eta_d = calculate_eta_2(p, w, z)

                # print(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
                flux.append(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

            ## Case X
            elif (p > 0) and (p < 1) and (z[w] == 0):
                # print('Case X')
                lambda_d = calculate_lambda_6(p)
                eta_d = calculate_eta_2(p, w, z)

                # print(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
                flux.append(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))

            ## Case XI
            elif (p > 1) and (z[w] >= 0) and (z[w] < p-1):
                # print('Case XI')
                lambda_d = 1
                eta_d = 1

                # print(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
                flux.append(calculate_flux(c2, c4, Omega, lambda_e, lambda_d, eta_d, p, w, z))
            
        # x -> time
        ttrans = period/(pi * adivR)
        mid = (len(x_values)/2) - 0.5

        for u in range(len(x_values)):
            if (u < mid):
                simulated_curve[u, 0] = (-1)*x_values[u]*ttrans/2
            elif (u >= mid):
                simulated_curve[u, 0] = x_values[u]*ttrans/2

        simulated_curve[:, 1] = flux[:]

        # Resampling   
        resampled_simulated_curve[:, 0] = observed_curve[:, 0]
        resampled_simulated_curve[:, 1] = scs.resample(simulated_curve[:, 1], observed_curve.shape[0])


        return resampled_simulated_curve

    def calculate_chi2(self, observed_curve, simulated_curve) -> float:
        chi2 = 0
        chi2 = sum(((observed_curve[:, 1] - simulated_curve[:, 1])**2)/(observed_curve[:, 2]**2) )
        return chi2




