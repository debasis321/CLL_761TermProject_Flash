import pandas as pd
import numpy as np

excelPath = "HYSYS_Props.xlsx"

COMP_DB = {"critical": pd.DataFrame(), "kij": pd.DataFrame(), "Cp_ideal": pd.DataFrame(), "heatFormation": pd.DataFrame()}

## HYSYS Property Data
COMP_DB["critical"] = pd.read_excel(excelPath, sheet_name="criticalProp")
COMP_DB["critical"].head()


# fix cpmpound name and Column name
COMP_DB["critical"].rename(columns={"CompName": "component"}, inplace=True)
# make it as index
COMP_DB["critical"].set_index("component", inplace=True)
# make comp name as lower case
COMP_DB["critical"].index = COMP_DB["critical"].index.str.lower()
# Fix accentricity factor and Column name
COMP_DB["critical"].rename(columns={"Acentricity": "omega"}, inplace=True)
# Fix BP value and Column name
COMP_DB["critical"].rename(columns={"B.P_C": "bp"}, inplace=True)
COMP_DB["critical"]["bp"] = COMP_DB["critical"]["bp"] + 273.14
# Fix Pc value and Column name
COMP_DB["critical"].rename(columns={"Pc_kg/cm2g": "pc"}, inplace=True)
COMP_DB["critical"]["pc"] = COMP_DB["critical"]["pc"] * 98066.5 + 101325  # convert to Pa Absolute
# fix Tc value and Column name
COMP_DB["critical"].rename(columns={"Tc_C": "tc"}, inplace=True)
COMP_DB["critical"]["tc"] = COMP_DB["critical"]["tc"] + 273.14   # to convert to K
# Fix Vc value and Column name
COMP_DB["critical"].rename(columns={"Vc_m3/kgmole": "vc"}, inplace=True) # in m3/kmol
# fix molecular weight and Column name
COMP_DB["critical"].rename(columns={"MW": "mw"}, inplace=True)

# lets ckeck the data
COMP_DB["critical"].head()


# load the binary ineteraction parameters
COMP_DB['kij'] = pd.read_excel(excelPath, sheet_name="kij")
# make all column and index as lower case
COMP_DB['kij'].columns = COMP_DB['kij'].columns.str.lower()
#set the fast column as index
COMP_DB['kij'].set_index('component', inplace=True)
# make comp name as lower case
COMP_DB['kij'].index = COMP_DB['kij'].index.str.lower()

COMP_DB['kij'].head()

# load the ideal gas heat capacity data
COMP_DB['Cp_ideal'] = pd.read_excel(excelPath, sheet_name="Cp_ideal")
# make all column and index as lower case
COMP_DB['Cp_ideal'].columns = COMP_DB['Cp_ideal'].columns.str.lower()
# set the fast column as index
COMP_DB['Cp_ideal'].set_index('component', inplace=True)
# make row names as lower case
COMP_DB['Cp_ideal'].index = COMP_DB['Cp_ideal'].index.str.lower()
# print(COMP_DB['Cp_ideal'].head())

### Load heat of formtion
# load the heat of formation data
COMP_DB['heatFormation'] = pd.read_excel(excelPath, sheet_name="heatFormation")
# make all column and index as lower case
COMP_DB['heatFormation'].columns = COMP_DB['heatFormation'].columns.str.lower()
# set the first column as index
COMP_DB['heatFormation'].set_index('component', inplace=True)
# make row names as lower case
COMP_DB['heatFormation'].index = COMP_DB['heatFormation'].index.str.lower()
# Constants for Peng-Robinson EOS
R = 8.314  # J/(mol*K)

# sanity checking
def check_sanity(comp):
    """
    Check the sanity of the stream composition
    """
    compName = tuple(comp.keys())
    # Check if all components in the stream are in the database
    if not all(item in COMP_DB["critical"].index for item in compName):
        raise ValueError("Some components are not in CRITICAL database.")
    if not all(item in COMP_DB["kij"].index for item in compName):
        raise ValueError("Some components are not in kij database.")
    if not all(item in COMP_DB["Cp_ideal"].index for item in compName):
        raise ValueError("Some components are not in Cp_ideal database.")
    
    # Check if the sum of the composition is 1
    if not np.isclose(sum(comp.values()), 1):
        raise ValueError("The sum of the composition is not equal to 1.")

    
    print("Sanity check completed successfully.")