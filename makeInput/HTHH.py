# Python script to create a set of input files for a Monte-Carlo sweep over input
# parameters so that Plumeria can be used to model the HTHH 15th January 2022
# eruption

# Import packages
import numpy as np

# Parameters to be swept over are the MER, exit velocity, mass fraction of added
# water, magma temperature, initial gas mass fraction of magma

# User-defined variables #######################################################

#### General inputs ##############################################
    
# Directory to store input files
fileDir = '/home/paulj/Documents/tonga2022/plumeria_wd/input/HTHH/'

# Number of simulations to be performed
nSims = 500

#### Weather inputs - if a weather file is provided, other inputs are not used
#### but still need to be provided. They can be meaningless

# Is a weather file to be read (y/n)?
metMarker = 'y'

# Name of weather file
metFile = '/home/paulj/Documents/tonga2022/plumeria_wd/input/soundings/HTHH.txt'

# Air temperature at vent (degC)
airTemp0 = 20.0

# Relative humidity of air
rHumAir = 0.0

# Lapse rate in troposphere (K/m) - needs to be -ive
tropoLapse = -0.0065

# Tropopause altitude (m)
tropoAlt = 11000.0

# Tropopause thickness (m)
tropoThick = 9000.0

# Lapse rate in the stratosphere (K/m) - should be +ive
stratoLapse = 0.001

# Wind speed (m/s)
windSp = 0.0

# Wind direction (bearing from N)
windDir = 90.0

#### Fixed vent properties ###########################

# Vent altitude (m)
ventAlt = 0.0

#### Fixed magma properties ##########################

# Specific heat of magma (J/kg K)
magSpHeat = 1000.0

# Magma density (DRE) (kg m-3)
magDens = 2500.0

#### Ranges of parameters to be swept over ##########
#### Define minimum and maximum values in an array

# MERs (kg /s) are sampled logarithmically - so the extremes are listed as
# powers of 10
logMERRange = np.array([5, 9])

# Exit velocity (m/s)
exitVelRange = np.array([300, 350])

# Mass fraction of added water
addH20Range = np.array([0.0, 0.4])

# Magma temperature (K)
magTempRange = np.array([1073, 1473])

# Initial magmatic gas mass fraction
gasMassFracRange = np.array([0.01, 0.05])

# Calculate random inputs and store in arrays ##################################
mer = np.log10(np.random.uniform(logMERRange[0], logMERRange[1], nSims))

exitVel = np.random.uniform(exitVelRange[0], exitVelRange[1], nSims)

addH2Ofrac = np.random.uniform(addH20Range[0], addH20Range[1], nSims)

magTemp = np.random.uniform(magTempRange[0], magTempRange[1], nSims)

gasMassFrac = np.random.uniform(gasMassFracRange[0], gasMassFracRange[1], nSims)

# In a loop, write each file

