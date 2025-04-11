from Flash import Q_Flash, PENG_ROBINSON, EOS

class Heater(Q_Flash):
    def __init__(fs: PENG_ROBINSON, dT=None, dP=None, beta=None, dQ=None):
        
        return super().__init__(fs, dT=dT, dP=dP, beta=beta, dQ=dQ)