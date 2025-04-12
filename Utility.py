import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import exp, log, sqrt, pi, atan, tan, cos, sin, sinh, cosh, asinh, acosh, atanh
import matplotlib.cm as cm

from Flash import Stream



def plot_TXY(stream: Stream,*, P, T_min, T_max, T_step):
    '''Plot Txy diagram for a binary system'''
    stream = copy.copy(stream)

    X = np.linspace(T_min, T_max, T_step)
    Y1, Y2 = [], []
    stream.P = P # Set the pressure for the stream
    for x in X:
        stream.T = x
        Y1.append(stream.x_i)  # Assuming stream.x_i is an array
        Y2.append(stream.y_i)  # Assuming stream.y_i is an array

    Y1 = np.array(Y1)  # Shape: (len(X), num_components)
    Y2 = np.array(Y2)

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot each component's liquid and vapor compositions vs temperature
    # get colur scale
    colors = cm.plasma(np.linspace(0, 1, Y1.shape[1]))

    for i in range(Y1.shape[1]):
        ax.plot(X, Y1[:, i], label=f'Liquid x[{i}]', linestyle='-', alpha=0.7, color=colors[i])
        ax.plot(X, Y2[:, i], label=f'Vapor y[{i}]', linestyle='--', alpha=0.7, color=colors[i])
        # ax.plot(Y1[:, i], X, label=f'Liquid x[{i}]', linestyle='-', alpha=0.7, color=colors[i])
        # ax.plot(Y2[:, i], X, label=f'Vapor y[{i}]', linestyle='--', alpha=0.7, color=colors[i])

    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Mole Fraction')
    ax.set_title(f'T-x-y Diagram at P: {P} Pa for {stream.stream.name}')
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.show()

def plot_PXY(stream: Stream,*, T, P_min, P_max, P_step):
    '''Plot Pxy diagram for a binary system'''
    stream = copy.copy(stream)

    X = np.linspace(P_min, P_max, P_step)
    Y1, Y2 = [], []
    stream.T = T # Set the temperature for the stream
    for x in X:
        stream.P = x
        Y1.append(stream.x_i)  # Assuming stream.x_i is an array
        Y2.append(stream.y_i)  # Assuming stream.y_i is an array

    Y1 = np.array(Y1)  # Shape: (len(X), num_components)
    Y2 = np.array(Y2)

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot each component's liquid and vapor compositions vs temperature
    # get colur scale
    colors = cm.plasma(np.linspace(0, 1, Y1.shape[1]))
    for i in range(Y1.shape[1]):
        ax.plot(X, Y1[:, i], label=f'Liquid x[{i}]', linestyle='-', alpha=0.7, color=colors[i])
        ax.plot(X, Y2[:, i], label=f'Vapor y[{i}]', linestyle='--', alpha=0.7, color=colors[i])
        # ax.plot(Y1[:, i], X, label=f'Liquid x[{i}]', linestyle='-', alpha=0.7, color=colors[i])
        # ax.plot(Y2[:, i], X, label=f'Vapor y[{i}]', linestyle='--', alpha=0.7, color=colors[i])

    ax.set_xlabel('Pressure, Pa')
    ax.set_ylabel('Mole Fraction')
    ax.set_title(f'P-x-y Diagram at T:{T} K for {stream.stream.name}')
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.show()