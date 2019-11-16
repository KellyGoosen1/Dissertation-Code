import time
from preisSeed import PreisModel
import numpy as np
import pandas as pd
import multiprocessing as mp

# Construct PreisModel object with default parameter
np.random.seed(19950412)

# 1. s() = I()
# 2. Compute s(y_obs) = y_obs
L = 480     # TimeHorizon
N = 10      # number of simulated y_obs
p_0 = 100
MCSteps = 10**5
N_A = 250

# 3. Sample N thetas from pi(theta)
# Prior dbn:
#   delta, tilda, mu U(0,0.1)
#   alpha U(0.1, 0.5)
#   lambda_0 U(0,200)
#   C_lam  U(0,20)

def priorSimulation(numberSim):

    # generate numberSim simultations of each parameter
    deltaSim = np.random.uniform(0, 0.5, numberSim)
    muSim = np.random.uniform(0, 0.1, numberSim)
    alphaSim = np.random.uniform(0.05, 0.5, numberSim)
    lambda0Sim = np.random.uniform(50, 250, numberSim)
    C_lambdaSim = np.random.uniform(0, 50, numberSim)
    delta_SSim = np.random.uniform(0, 0.05, numberSim)

    # Store and return simultations
    priorDf = np.column_stack((deltaSim, muSim, alphaSim, lambda0Sim, C_lambdaSim, delta_SSim))
    return pd.DataFrame(priorDf, columns=("delta", "mu", "alpha", "lambda_0", "C_lambda", "delta_S"))

priorMatrix = priorSimulation(N)

# randomise seed
seed = int(np.random.uniform(0, 1000000))
np.random.seed(seed)

print(time.asctime())

def PreisSim(params, N_A=N_A, p_0=p_0, L=L, MCSteps=MCSteps):
    params = priorMatrix.iloc[i]
    # initialize preis model object with specified parameters
    p = PreisModel(N_A=N_A,
                   delta=params["delta"],
                   lambda_0=params["lambda_0"],
                   C_lambda=params["C_lambda"],
                   delta_S=params["delta_S"],
                   alpha=params["alpha"],
                   mu=params["mu"],
                   p_0=p_0,
                   T=L,
                   MC=MCSteps)

    # Start model
    p.simRun()
    p.initialize()

    # Simulate price path for T=L time-steps
    p.simulate()

    return pd.DataFrame(p.intradayPrice)



# print(time.asctime())
