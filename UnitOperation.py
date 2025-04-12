from Flash import Q_Flash, PENG_ROBINSON, EOS

class Heater(Q_Flash):
    """
    A class representing a Heater unit operation, inheriting from Q_Flash.

    The Heater class models a heating process where the enthalpy change, 
    pressure change, and temperature change are calculated based on the 
    inlet and outlet stream properties.

    Attributes:
        fs (PENG_ROBINSON): The feed stream object, representing the inlet stream.
        duty (float): The heat duty of the heater, calculated as the product 
                      of the enthalpy difference and the flow rate of the feed stream.
        dP (float): The pressure difference between the outlet and inlet streams.
        dT (float): The temperature difference between the outlet and inlet streams.

    Methods:
        __init__(fs: PENG_ROBINSON, dT=None, dP=None, beta=None, dQ=None):
            Initializes the Heater object with the feed stream and optional parameters.
    """
    def __init__(self, fs: PENG_ROBINSON, dT=None, dP=None, beta=None, dQ=None):
        super().__init__(fs, dT=dT, dP=dP, beta=beta, dQ=dQ)
        # getting ps from super class
        self.fs = fs
        self.dH = (self.ps.h - self.fs.h) * self.fs.flowrate
        self.dP = self.ps.P - self.fs.P
        self.dT = self.ps.T - self.fs.T
        # Add any Heater-specific initialization below if needed

    def basic_info(self):
        """
        Print basic information about the Heater unit operation.
        """
        print(f"### Heater Unit Operation ###")
        print(f"Feed Stream: {self.fs.stream.name}")
        print(f"Outlet Stream: {self.ps.stream.name}")
        print(f"Heat Duty: {round(self.dH/1e3, 2)} kJ/h")
        print(f"Pressure Drop: {round(self.dP, 2)} Pa")
        print(f"Temperature Change: {round(self.dT, 2)} K")
        print("#############################")

class Mixer(Q_Flash):
    """
    A class representing a Mixer unit operation, inheriting from Q_Flash.

    The Mixer class models a mixing process where the enthalpy change, 
    pressure change, and temperature change are calculated based on the 
    inlet and outlet stream properties.

    Attributes:
        fs (PENG_ROBINSON): The feed stream object, representing the inlet stream.
        duty (float): The heat duty of the mixer, calculated as the product 
                      of the enthalpy difference and the flow rate of the feed stream.
        dP (float): The pressure difference between the outlet and inlet streams.
        dT (float): The temperature difference between the outlet and inlet streams.

    Methods:
        __init__(fs: PENG_ROBINSON, dT=None, dP=None, beta=None, dQ=None):
            Initializes the Mixer object with the feed stream and optional parameters.
    """
    def __init__(self, fs1: PENG_ROBINSON, fs1: PENG_ROBINSON):
        super().__init__(fs, dT=dT, dP=dP, beta=beta, dQ=dQ)
        # getting ps from super class
        self.fs = fs
        self.dH = (self.ps.h - self.fs.h) * self.fs.flowrate
        self.dP = self.ps.P - self.fs.P
        self.dT = self.ps.T - self.fs.T
        # Add any Mixer-specific initialization below if needed

    def basic_info(self):
        """
        Print basic information about the Mixer unit operation.
        """
        print(f"### Mixer Unit Operation ###")
        print(f"Feed Stream: {self.fs.stream.name}")
        print(f"Outlet Stream: {self.ps.stream.name}")
        print(f"Heat Duty: {round(self.dH/1e3, 2)} kJ/h")
        print(f"Pressure Drop: {round(self.dP, 2)} Pa")
        print(f"Temperature Change: {round(self.dT, 2)} K")
        print("#############################")