# %% [markdown]
# 

# %%
# Basic Code for CLL-761 Project by Debasis

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import math
import unittest
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize as opt
from scipy.optimize import fsolve
from scipy.optimize import minimize
from scipy.optimize import root
from scipy.optimize import brentq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy
import inspect

from Solver import Solver
from Database import COMP_DB, R
#

# %%


# %%
# Create Streaam Class
# create a class for the stream
class Stream():
    """
    A class to represent a chemical stream with its properties.
    Attributes:
    ----------
    name : str
        The name of the stream.
    composition : dict, optional
        A dictionary containing the components as keys and their mole fractions as values.
    temperature : float, optional
        The temperature of the stream in Kelvin.
    pressure : float, optional
        The pressure of the stream in Pascals.
    flowrate : float, optional
        The flowrate of the stream.
    Methods:
    -------
    get_mole_fraction():
        Returns the mole fractions of the components in the stream.
    get_temperature():
        Returns the temperature of the stream.
    get_pressure():
        Returns the pressure of the stream.
    get_flowrate():
        Returns the flowrate of the stream.

    """
    def __init__(self, name=None, composition=None, temperature=None, pressure=None, flowrate=None):
        self.name = name
        self.components = list(composition.keys())
        self.mole_fractions = list(composition.values())
        self.temperature = temperature
        self.pressure = pressure
        self.flowrate = flowrate

    def get_mole_fraction(self):
        return self.mole_fractions

    def get_temperature(self):
        return self.temperature

    def get_pressure(self):
        return self.pressure

    def get_flowrate(self):
        return self.flowrate
    

    
# %%
# Create  a Stream object
# stream1 = Stream(composition={"benzene": 0.2, "toluene": 0.3, "xylene": 0.5}, temperature=350, pressure=101325, flowrate=1000)

# %%
# Create a empty Class as EOS
class EOS:
    def __init__(self):
        pass
class UOP:
    def __init__(self):
        pass

# %%
# Class EOS for Peng-Robinson Equation of State
class PENG_ROBINSON (EOS):
    def __init__(self, stream: Stream, compData=COMP_DB, solver=Solver()):
        # attaching the class args
        self.stream = stream
        self.compData = compData
        self.solver = solver
        #f loading of stream data 
        self.components = stream.components
        self._z_i = stream.get_mole_fraction()
        self._T = stream.get_temperature()
        self._P = stream.get_pressure()
        self.flowrate = stream.get_flowrate()
        self.N_component = len(self.components)
        # self.volumeFlowRate = self.calculate_total_volume_flowrate()
        # loading of critical properties and calculation of parameters in init method
        self.load_critical_data()
        self.calc_parameters()
        # Phase composition: intialize with toatal mle raction then update with Flash calculation
        
        self.x_i, self.y_i, self.K, self.vf = self.estimate_xi_yi_K_beta() # initial guess for liquid and vapor compositions
        self.A, self.B = None, None
        self.Z_l = None
        self.Z_v = None
        self.phi_l = None
        self.phi_v = None

        self.h_l = None
        self.h_v = None
        self.h = None

        self.adiabatic_flash()
        self.update_stream_enthalpy()
        # self.phase = self.phaseCheck(self.K)
        # print(self.x_i, self.y_i, self.vf)
        # self.x_i, self.y_i, self.vf = self.z_i, self.z_i, 0.5 # initial guess for liquid and vapor compositions
        # number of components
        # run calc_parametrs on change of value of T,P, composition 
    
    @property
    def T(self):
        return self._T

    @T.setter
    def T(self, new_T):
        self._T = new_T
        self.on_state_change()

    @property
    def P(self):
        return self._P

    @P.setter
    def P(self, new_P):
        self._P = new_P
        self.on_state_change()

    @property
    def z_i(self):
        return self._z_i

    @z_i.setter
    def z_i(self, new_z):
        self._z_i = new_z
        self.on_state_change()

    def on_state_change(self):
        """Update parameters when temperature changes."""
        self.calc_parameters()
        self.adiabatic_flash()
        self.update_stream_enthalpy()    
    


    def load_critical_data(self):
        '''Loading of critical properties from data bank at imitialization'''
        compCriticalProp = self.compData["critical"].loc[self.components]

        self.omega = compCriticalProp['omega'].to_numpy()
        self.tc = compCriticalProp['tc'].to_numpy()
        self.pc = compCriticalProp['pc'].to_numpy()
        # self.zc = self.compCriticalProp['zc'].to_numpy()
        # self.vc = self.compCriticalProp['vc'].to_numpy()
        self.bp = compCriticalProp['bp'].to_numpy()
        self.mw = compCriticalProp['mw'].to_numpy()
        # self.k_ij = calculate_kij(self.tc) # binary interaction parameter
        self.k_ij = np.array(self.compData['kij'][self.components].loc[self.components])
        # self.k_ij = np.array([[0.0, 0.0, 0],
        #                       [0.0, 0.0, 0],
        #                       [0.0, 0.0, 0]])
        # self.tc = np.array([562.2, 591.8])  # K
        # self.pc = np.array([4890000.0, 4108000.0])  # Pa
        # self.omega = np.array([0.212, 0.263])
        # print(self.k_ij)
        # loading Cp ideal values
        self.idealCp_coeff = self.compData["Cp_ideal"].loc[self.components].to_numpy()
        self.compHf = self.compData["heatFormation"].loc[self.components].to_numpy()
        # Calculate A, B, parameters

    def calc_parameters(self):
        """Calculate EOS parameters for Peng-Robinson."""
        self.tr = self.T / self.tc
        self.pr = self.P / self.pc
        self.m = (0.37464 + 1.54226*self.omega - 0.26992*self.omega**2)
        self.alpha = (1 + self.m * (1 - np.sqrt(self.tr)))**2
        self.a = 0.45724 * (R**2 * self.tc**2 * self.alpha) / self.pc
        self.b = 0.07780 * R * self.tc / self.pc
        self.a_ij = np.sqrt(np.outer(self.a, self.a)) * (1 - self.k_ij)
        print('##### Parameters Calculation:')
        print(f'Temperature: {self.T}, Pressure: {self.P}')
        print(f'Tr: {self.tr}, alpha: {self.alpha}, a: {self.a}, b: {self.b}')
        print('#####')
    
    def basic_info(self):
        '''print Basic data of a stream'''
        print(f'### Post Flash Results for Stream:{self.stream.name} ###')
        print(f'T:, {round(self.T, 2)} K, P: {round(self.P, 2)} Pa, F: {round(self.flowrate, 2)} mol/hr') 
        print(f'Z_l: {round(self.Z_l,4)}, Z_v: {round(self.Z_v,4)}, phi_l: {np.round(self.phi_l,2)}, phi_v: {np.round(self.phi_v,2)}')
        print(f'K: {np.round(self.K, 4)}')
        print(f'VF: {round(self.vf,2)}, xi: {np.round(self.x_i,4)}, yi: {np.round(self.y_i,4)}')
        print(f'h_l: {np.round(self.h_l,2)} J/mol, h_v: {np.round(self.h_v,2)} J/mol, h: {np.round(self.h,2)} J/mol')
        print('###')


    def residual_enthalpy_mixture(self, x_y, Z):
        """Calculate residual enthalpy for the mixture.
        Args:
            z (array): Mole fractions of components in the mixture.
            Z (float): Compressibility factor.
        Returns:
        float: Residual enthalpy of the mixture.
            """
        z = x_y
        T = self.T
        P = self.P
        Tc = self.tc
        Pc = self.pc
        omega = self.omega
        kij = self.k_ij
        Tr = T / Tc
        m = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
        alpha = (1 + m * (1 - np.sqrt(Tr)))**2
        a = 0.45724 * R**2 * Tc**2 * alpha / Pc
        b = 0.07780 * R * Tc / Pc

        # Mixing rules
        a_ij = np.sqrt(np.outer(a, a)) * (1 - kij)
        a_mix = np.dot(z, np.dot(a_ij, z))
        b_mix = np.dot(z, b)

        # da/dT for each component
        dalpha_dT = -m * (1 / np.sqrt(Tr)) * (1 / Tc) * (1 + m * (1 - np.sqrt(Tr)))
        da_dT = 0.45724 * R**2 / Pc * (2 * Tc * alpha * (-1 / T) + Tc**2 * dalpha_dT)

        # Mixture da/dT
        da_ij = np.outer(da_dT, a) + np.outer(a, da_dT)
        da_ij_dT = 0.5 * da_ij / np.sqrt(np.outer(a, a)) * (1 - kij)
        da_mix_dT = np.dot(z, np.dot(da_ij_dT, z))

        # Residual enthalpy [J/mol]
        term1 = (Z - 1) * R * T
        B = b_mix * P / (R * T)
        eps = 0
        term2 = ((T * da_mix_dT - a_mix) / (2 * np.sqrt(2) * b_mix + eps)) 
        # term2 = 2.078 * (1 + m) * np.sqrt(alpha) 
        term3 = np.log((Z + (1 + np.sqrt(2)) * B+ eps) / (Z + (1 - np.sqrt(2)) * B + eps))
        H_res = term1 + term2 * term3
        print('##### Residual Enthalpy Calculation:')
        print(f'Term1: {term1}, Term2: {term2}')
        print((f"Residual Enthalpy: {H_res} J/mol"))
        print('#####')
        return H_res  # J/mol

    def ideal_enthalpy_mixture(self, x_y):
        """Calculate ideal enthalpy for the mixture."""
        z = x_y
        coeff = self.idealCp_coeff
        T_ref = 273+ 25  # Reference temperature in K (25 °C)    
        T = self.T -T_ref
        # ideal_H = coeff[:, 0] + coeff[:, 1] * T + coeff[:, 2] * T**2  + coeff[:, 3] * T**3  + coeff[:, 4] * T**4 + coeff[:, 5] * T**5  
        ideal_H = coeff[:, 0]*T + coeff[:, 1] * T**2/2 + coeff[:, 2] * T**3/3  + coeff[:, 3] * T**4/4  + coeff[:, 4] * T**5/5 + coeff[:, 5] * T**6/6 

        ideal_H = ideal_H * 3  # Convert to J/mol (HYSYS calculation gives in Cal/mol)
        print('##### Ideal Enthalpy Calculation:')
        print(f'Ideal Enthalpy: {ideal_H} J/mol')
        print('#####')
        return np.dot(z, ideal_H)  # Mole fraction-weighted sum
    
    def heatOfformation_mixture(self, x_y):
        """Calculate heat of formation for the mixture."""
        z = x_y
        # Hf = self.compHf['Hf'].to_numpy()
        hf = self.compHf
        return np.dot(z, hf)  # J/mol
    
    def update_stream_enthalpy(self):
        """Calculate stream enthalpy."""

        # h_l = self.residual_enthalpy_mixture(self.x_i, self.Z_l) + self.ideal_enthalpy_mixture(self.x_i) if not np.all(self.x_i==0) else 0
        # h_v = self.residual_enthalpy_mixture(self.y_i, self.Z_v) + self.ideal_enthalpy_mixture(self.y_i)if not np.all(self.y_i==0) else 0
        h_l = self.residual_enthalpy_mixture(self.x_i, self.Z_l) + self.heatOfformation_mixture(self.x_i) + self.ideal_enthalpy_mixture(self.x_i) if not np.all(self.x_i==0) else 0
        h_v = self.residual_enthalpy_mixture(self.y_i, self.Z_v) + self.heatOfformation_mixture(self.y_i) + self.ideal_enthalpy_mixture(self.x_i) if not np.all(self.y_i==0) else 0
        h = h_v * self.vf + h_l * (1 - self.vf)
        self.h_l = h_l
        self.h_v = h_v
        self.h = h
        print('##### Enthalpy Calculation:')
        print(f'Liquid Enthalpy: {h_l} J/mol')
        print(f'Vapor Enthalpy: {h_v} J/mol')
        print(f'Mixture Enthalpy: {h} J/mol')
        print('#####')
        return h_l, h_v, h  # J/mol
    
    def estimate_xi_yi_K_beta(self):
        """Estimate initial guess for liquid and vapor compositions using Wilson's method."""
        # Calculate K-values using Wilson's method
        ln_K =np.log(1/self.pr) + 5.37 * (1 + self.omega) * (1 - 1 / self.tr)
        K = np.exp(ln_K)  # K-values for each component
        # Check if the system is single phase or two phase
        if np.all(K < 1):
            print("Single phase: Liquid")
            # x_i = self.z_i
            # y_i = np.zeros_like(self.z_i)
            # if single phase keep z as liquid and vapor phase composition to test in PR
            x_i = self.z_i
            y_i = self.z_i
            return x_i, y_i, K, 0
        
        if np.all(K >= 1):
            print("Single phase: Vapor")
            # x_i = np.zeros_like(self.z_i)
            # y_i = self.z_i
            # if single phase keep z as liquid and vapor phase composition to test in PR
            x_i = self.z_i
            y_i = self.z_i
            return x_i, y_i, K, 1
        
        else:
            print("Two phase")
            beta_min = max(-1 / (np.max(K) - 1), 0)  # Ensure beta ≥ 0
            beta_max = min(1 / (1 - np.min(K)), 1)  # Ensure beta ≤ 1
            # Two phase system
            beta_solver = self.solver.beta_solver
            beta = beta_solver(z=self.z_i, K=K, tol=1e-8, max_iter=100)
            x_i = self.z_i / (1 + beta * (K - 1))
            y_i = K * x_i

            # Normalize to sum to 1
            x_i /= np.sum(x_i)
            y_i /= np.sum(y_i)
        print(f'Willson: xi:{x_i}, yi:{y_i}, K: {K}, beta:{beta}')
        return x_i, y_i, K, beta
            


    def calc_AB(self, z):
        """Compute mixture parameters A and B using mixing rules."""
        a_mix = np.sum(np.outer(z, z) * self.a_ij)
        b_mix = np.sum(z * self.b)
        A = a_mix * self.P / (R**2 * self.T**2)
        B = b_mix * self.P / (R * self.T)
        self.A = A
        self.B = B
        print(f'A: {A}, B: {B}, a_mix: {a_mix}, b_mix: {b_mix}')
        return A, B, a_mix, b_mix

    def calc_Z(self, A, B):
        """Solve cubic Peng-Robinson EOS for Z-factor."""
        Z_l, Z_v = self.solver.Z_solver(A, B, solverName='np_roots')
        print('##### Z-factor Calculation:')
        print(f'Z_l: {Z_l}, Z_v: {Z_v}')
        print('#####')
        return Z_l, Z_v
        

    def calc_fugacity_coeff(self, Z, x_or_y, A, B, a_mix, b_mix):
        """Calculate fugacity coefficient φ for each component in given phase composition x_or_y."""
        eps = 0  # Small positive value to avoid division by zero or log of non-positive values
        sqrt2 = np.sqrt(2)
        n = len(x_or_y)
        phi = np.zeros(n)

        term1 = self.b / (b_mix + eps) * (Z - 1)
        term2 = - np.log(Z - B + eps)
        term3 = - (A / (2 * sqrt2 * B + eps))
        term4 = ((2 * np.dot(x_or_y, self.a_ij) / (a_mix + eps)) - (self.b / (b_mix + eps)))
        term5 = np.log((Z + (1 + sqrt2) * B + eps) / (Z + (1 - sqrt2) * B + eps))

        ln_phi = term1 + term2 + term3 * term4 * term5
        phi = np.exp(ln_phi)
        print('##### Fugacity Coefficient Calculation:')
        print(f'Z: {Z}, x_or_y: {x_or_y}, A: {A}, B: {B}, a_mix: {a_mix}, b_mix: {b_mix}')
        # print(f'term1: {term1}, term2: {term2}, term3: {term3}, term4: {term4}, term5: {term5}')
        print(f'ln_phi: {ln_phi}')
        print(f'phi: {phi}')    
        print('#####')
        
        return phi
    # run update K on change of x,y,z,T,P
    def update_K(self):
        """Compute equilibrium K-values from fugacity coefficients."""
        # Compute A, B, and Z for liquid phase
        A_l, B_l, a_l, b_l = self.calc_AB(self.x_i)
        Z_l, _ = self.calc_Z(A_l, B_l)
        phi_l = self.calc_fugacity_coeff(Z_l, self.x_i, A_l, B_l, a_l, b_l)

        # Compute A, B, and Z for vapor phase
        A_v, B_v, a_v, b_v = self.calc_AB(self.y_i)
        _, Z_v = self.calc_Z(A_v, B_v)
        phi_v = self.calc_fugacity_coeff(Z_v, self.y_i, A_v, B_v, a_v, b_v)

        # Calculate K-values
        K = phi_l / phi_v
        print('##### K-value Calculation:')
        print(f"phi_l:{phi_l}, phi_v:{phi_v}")
        print(f'K:{K}')
        print('#####')
        print('####fugacity balance:')
        print('Vapor fugacity:',  self.y_i * phi_v)
        print('Liquid fugacity:', self.x_i * phi_l)
        print('####')
        self.Z_l = Z_l
        self.Z_v = Z_v
        self.phi_l = phi_l
        self.phi_v = phi_v
        self.K = K
        


    # detection of phase
    def phaseCheck(self):
        K = self.K
        # Check if the system is single phase or two phase
        if np.all(K < 1):
            print("Single phase: Liquid")
            return "LIQUID"
        elif np.all(K >= 1):
            print("Single phase: Vapor")
            return "VAPOR"
        else:
            print("Two phase")
            return "TWOPHASE"    
    
    def adiabatic_flash(self):      
        self.x_i, self.y_i, self.K, self.vf = self.estimate_xi_yi_K_beta()
        # self.update_K()

        # if self.phaseCheck() == "LIQUID":
        #     print("Single phase: Liquid")
        #     self.x_i = self.z_i
        #     self.y_i = np.zeros_like(self.z_i)
        #     return
        
        # if self.phaseCheck() == "VAPOR":
        #     print("Single phase: Vapor")
        #     self.x_i = np.zeros_like(self.z_i)
        #     self.y_i = self.z_i
        #     return
        
        # if self.phaseCheck() == "TWOPHASE":
        #     # Two phase flash 
        #     print("Two phase flash") 
        beta  = self.vf # initial guess for beta
        iter = 0
        error = 1
        beta_solver = self.solver.beta_solver

        while (error > 1e-6) and (iter < 100):  
            self.update_K()
            if self.phaseCheck() == "LIQUID":
                print("Single phase: Liquid")
                self.x_i = self.z_i
                self.y_i = np.zeros_like(self.z_i)
                self.vf = 0
                break
            if self.phaseCheck() == "VAPOR":
                print("Single phase: Vapor")
                self.x_i = np.zeros_like(self.z_i)
                self.y_i = self.z_i
                self.vf = 1
                break
            if self.phaseCheck() == "TWOPHASE":
                print("Two phase flash") 
                z = self.z_i
                K = self.K
                beta = beta_solver(z=z, K=K, tol=1e-8, max_iter=100)
                x = self.z_i / (1 + beta * (K - 1))
                y = K * x
                # normalize and limit the values to 0-1
                x = x / np.sum(x)
                y = y / np.sum(y)
                # Calculate error
                error = np.sum(np.abs(x - self.x_i))  + np.sum(np.abs(y - self.y_i)) + np.abs(beta - self.vf)
                
                # cases for beta values
                if beta <=0: # liquid only
                    self.x_i = self.z_i
                    self.y_i = np.zeros_like(self.z_i)
                    self.vf = 0
                    print(f'Two phase to Liquid Phase: iter:{iter} x:{self.x_i}, y:{self.y_i}, beta:{self.vf}, error:{error}')
                    break
                if beta >=1: # vapor only
                    self.x_i = np.zeros_like(self.z_i)
                    self.y_i = self.z_i
                    self.vf = 1
                    print(f'Two phase to Vapor Phase: iter:{iter} x:{self.x_i}, y:{self.y_i}, beta:{self.vf}, error:{error}')
                    break
                if beta >0 and beta <1: # two phase
                    self.x_i = x
                    self.y_i = y
                    self.vf = beta
                    print(f'two-phase flash iter:{iter} x:{self.x_i}, y:{self.y_i}, beta:{self.vf}, error:{error}')
                    iter += 1
                    continue

        return
        
class Q_Flash(UOP):
    """
    This class implements the non-adiabatic flash calculation for a given mixture.
    It uses the Peng-Robinson equation of state to calculate the phase equilibrium.
    """

    def __init__(self, fs: PENG_ROBINSON, dT=None, dP=None, beta=None, dQ=None):
        """"""
        self.fs = copy.copy(fs) # feed stream object
        self.ps = None # product stream object

        if dT is not None and dP is not None :
            # calculate dQ from dT and dP
            self.ps = self.calc_on_dT_dP(dT, dP)
        elif dP is not None and beta is not None:
            # calculate dQ from beta and dP
            self.ps = self.calc_on_beta_dP(beta, dP)
        elif dT is not None and beta is not None:
            # calculate dQ from beta and dP
            self.ps = self.calc_on_beta_dT(beta, dT)

        elif dQ is not None and dP is not None:
            # calculate dT from dQ and dP
            self.ps = self.calc_on_dQ_dP(dQ, dP)
        elif dQ is not None and dT is not None:
            # calculate dT from dQ and dP
            self.ps = self.calc_on_dQ_dT(dQ, dT)
        elif dQ is not None and beta is not None:
            # calculate dT from dQ and beeta
            self.ps = self.calc_on_dQ_beta(beta, dQ)

        else:
            raise ValueError("Invalid input parameters. Please provide proper combination of dQ, dT, and dP or beta.")
            

    def bubble_dew_temp_approx(self):
        P=self.fs.P
        x=self.fs.x_i
        Tc=self.fs.tc
        Pc=self.fs.pc
        omega=self.fs.omega

        A = 5.37 * (1 + omega)
        num = np.sum(x * Tc * A)
        denom_bubble = np.sum(x * (np.log(Pc / P) + A))
        denom_dew = np.sum(x * (np.log(Pc / P) + A))
        T_bubble = num / denom_bubble
        T_dew = num / denom_dew
        
        return T_bubble, T_dew

    def calc_on_dT_dP(self, dT, dP):
        """Calculate the heater duty."""
        fs = copy.deepcopy(self.fs) # copy the stream object
        
        fs_h = fs.h # get feed total enthalpy. as it is a float(immutable) object, it will not change with the original stream object
        fs_F = fs.flowrate # get the feed flowrate
        # update the stream temperature and pressure
        fs.P = fs.P - dP
        fs.T = fs.T + dT
        ps = fs
        # print('##### calc_dQ_on_dT_dP:')
        # print(f'Feed Stream: {fs_h} J/mol, Product Stream: {ps_h} J/mol')
        # print(f'Q: {dQ} J/hr')
        # print(f'Feed Stream Flowrate: {fs_F} mol/hr')
        # print('#####')
        return ps

    def calc_on_dQ_dP(self, dQ, dP):
        """Calculate the heater duty."""
        fs = copy.deepcopy(self.fs) # copy the stream object
        ps = copy.deepcopy(self.fs) # copy the stream object

        # update the stream temperature and pressure
        ps.P = fs.P - dP
        ps_h = fs.h + dQ / fs.flowrate # get the new stream enthalpy target
        # iteratefor T to match the enthalpy
        def enthalpy_gap(T):
            ps.T = T
            print(ps.h)
            ps_h_new = ps.h # get the new stream enthalpy
            return ps_h_new - ps_h # return the difference
        
        T_min, T_max = 10, 500
        brentq(enthalpy_gap, T_min, T_max, xtol=1e-6) # solve for T using brentq method
        
        dT = ps.T - fs.T
        # print('##### calc_dT_on_dQ_dP')
        # print(f'Inlet Temperature: {fs.T}')
        # print(f'Outlet Temperature: {ps.T}')
        return ps

    def calc_on_beta_dP(self, beta, dP):
        """Calculate the heater duty to meet a vapor fraction and pressure drop"""
        fs = copy.copy(self.fs)
        ps = copy.copy(self.fs)
        ps.P = fs.P - dP

        def beta_gap(T):
            ps.T = T
            return ps.vf - beta
        
        # T_min, T_max = self.bubble_dew_temp_approx()
        T_min, T_max = 350, 450
        T_guess = 300
        brentq(beta_gap, T_min, T_max, xtol=1e-6) # solve for T using brentq method
        # res = minimize(beta_gap, x0=T_guess, method='BFGS')
        # print('##### calc_dQ_on_beta_dP')
        # print(f'Inlet Temperature: {fs.T}')
        # print(f'Outlet Temperature: {ps.T}')
        return ps
    
   
    def calc_on_beta_dT(self, beta, dT):
        """Calculate the heater duty to meet a vapor fraction and temperature change"""
        fs = copy.deepcopy(self.fs)
        ps = copy.deepcopy(self.fs)
        ps.T = fs.T + dT

        def beta_gap(P):
            ps.P = P
            return ps.vf - beta
        
        P_min, P_max = 1e4, 1e7
        brentq(beta_gap, P_min, P_max, xtol=1e-6)
        return ps
    
    def calc_on_dQ_dT(self, dQ, dT):
        """Calculate the heater duty to meet a vapor fraction and temperature change"""
        fs = copy.deepcopy(self.fs)
        ps = copy.deepcopy(self.fs)
        ps.T = fs.T + dT

        def enthalpy_gap(P):
            ps.P = P
            ps_h_new = ps.h
            return ps_h_new - (fs.h + dQ / fs.flowrate)
        
        P_min, P_max = 10, 10000000
        brentq(enthalpy_gap, P_min, P_max, xtol=1e-6)
        return ps
    
    def calc_on_dQ_beta(self, beta, dQ):
        """Calculate the heater duty to meet a vapor fraction and temperature change"""
        fs = copy.deepcopy(self.fs)
        ps = copy.deepcopy(self.fs)

        def beta_enthalpy_gap(P, T):
            ps.P, ps.T = P, T
            beta_new = ps.vf
            ps_h_new = ps.h
            return [(beta_new - beta),
                    (ps_h_new - (fs.h + dQ / fs.flowrate))]
        
        P_min, P_max = 10, 10000000
        T_min, T_max = 273, 1000
        P_guess, T_guess =10325, 300
        root(beta_enthalpy_gap, [P_guess, T_guess], method='hybr', bounds=[(P_min, P_max), (T_min, T_max)], xtol=1e-6)
        
        return ps


