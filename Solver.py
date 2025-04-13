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

# %%
class Solver:
    def __init__(self):
        pass

    def plot_cubic_PR(self, A, B):
        cubic_PR = lambda Z: Z**3 - (1 - B) * Z**2 + (A- 3 * B**2 -2 * B) * Z - (A * B - B**2 - B**3)
        Z = np.linspace(-0.5, 2, 100)
        plt.plot(Z, cubic_PR(Z), label='Cubic PR')
        plt.axhline(0, color='r', linestyle='--', label='y=0')
        plt.xlabel('Z')
        plt.ylabel('Cubic PR')
        plt.title('Cubic PR Equation value vs Z')
        plt.legend()
        plt.grid()
        plt.show()

    def Z_solver(self, A, B, tol=1e-8, max_iter=100, solverName='fsolve'):
        """Solve cubic Peng-Robinson EOS for Z-factor."""
        def cubic_PR(Z):
            return Z**3 - (1 - B) * Z**2 + (A- 3 * B**2 -2 * B) * Z - (A * B - B**2 - B**3)
        
        
        


        

        if solverName == 'fsolve':
            Z_l = fsolve(lambda Z: cubic_PR(Z), x0=0, xtol=1e-7 )
            Z_v = fsolve(lambda Z: cubic_PR(Z), x0=1, xtol=1e-7 )
            plot_cubic_PR()
            

        elif solverName == 'brentq':
            Z_l = brentq(lambda Z: cubic_PR(Z), 0, 0.5)
            Z_v = brentq(lambda Z: cubic_PR(Z), 0.5, 2)

        elif solverName == 'root':
            result = root(lambda Z: cubic_PR(Z), [0, 1], method='broyden1') 
            Z_l = min(result.x)
            if len(result.x) == 0:
                raise ValueError("No real roots found for Z.")
            Z_v = max(result.x)

        elif solverName == 'np_roots':
            coeffs = [1, - (1 - B) , (A- 3 * B**2 -2 * B),  - (A * B - B**2 - B**3)]
            roots = np.roots(coeffs)
            # Filter to keep only real roots
            print(f'Z roots: {roots}')
            real_roots = roots[np.isreal(roots)].real
            # get the min root which is >0.01 making physical sence
            # real_roots = real_roots[real_roots > 0.01]
            Z_l = min(real_roots)
            Z_v = max(real_roots)

        elif solverName == 'minimize':
            result = minimize(lambda Z: 1 - cubic_PR(Z), [0, 1], method='BFGS') 
            Z_l = min(result.x)
            if len(Z_l) == 0:
                raise ValueError("No real roots found for Z.")
            Z_v = max(result.x)
        elif solverName == 'newton':
        # Newton's method iteration
            for guess in [0, 1]:
                Z = guess
                for _ in range(max_iter):
                    f = cubic_PR(Z)
                    df_dZ = 3 * Z**2 - 2 * (1 - B) * Z + A
                    if abs(df_dZ) < tol or abs(f) < tol:
                        break
                    Z -= f / df_dZ  # Newton step
                if abs(f) < tol and guess == 0:
                    Z_l = Z
                elif abs(f) < tol and guess == 1:
                    Z_v = Z
                else: 
                    raise ValueError("No real roots found for Z. or convergence failed.")        

        else:
            raise ValueError("Unsupported solver name. Use 'fsolve' or 'brentq'.")

        # print(f'Z_l: {Z_l}, Z_v: {Z_v}')
        return Z_l, Z_v  # Z_l (liquid) and Z_v (vapor)
    
    def beta_solver(self, z, K, tol=1e-6, max_iter=1000, solverName='newton'):

        def RR(beta, z, K):
            return np.sum(z * (K - 1) / (1 + beta * (K - 1)))
        
        def dRR_dbeta(beta, z, K):
            return -np.sum(z * (K - 1) ** 2 / (1 + beta * (K - 1)) ** 2)
        
        # check for feasible solution domain
        def beta_fasibility(z, K):
            '''If there is a beta  solution value between 0 and 1 then return True'''
            return RR(beta=0, z=z, K=K) >0 and RR(beta=1, z=z, K=K)<0
        
        def plot_RRvsBeta(z, K):
            beta = np.linspace(0,1,10)
            rr = np.array([])
            for _ in beta:
                rr = np.append(rr, RR(_, z, K))
                
            plt.plot(beta, rr)
            plt.show()

        # plot_RRvsBeta(z, K)

        if solverName == 'newton':
            """Solves for beta using Newton's method."""
            # Define beta bounds
            beta_min = max(-1 / (np.max(K) - 1), 0)  # Ensure beta ≥ 0
            beta_max = min(1 / (1 - np.min(K)), 1)  # Ensure beta ≤ 1

            # Initial guess: midpoint of bounds
            beta = (beta_min + beta_max) / 2

            # Newton's method iteration
            for _ in range(max_iter):
                f_beta = RR(beta, z, K)
                df_beta = dRR_dbeta(beta, z, K)

                if  abs(df_beta) < tol or abs(f_beta) < tol:
                    break  # Converged

                beta -= f_beta / df_beta  # Newton step
                
            print(f'Calculated beta:{beta}')
            beta = max(beta_min, min(beta_max, beta))# Ensure beta remains within valid bounds

            return beta
        
