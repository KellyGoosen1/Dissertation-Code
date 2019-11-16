import numpy as np
import pandas as pd
import time
import pathos.multiprocessing as mp

# Prepare data
# Construct PreisModel object with default parameter

def priorSimulation(numberSim):
    # generate numberSim simultations of each parameter
    deltaSim = np.random.uniform(0, 0.1, numberSim)
    muSim = np.random.uniform(0, 0.1, numberSim)
    alphaSim = np.random.uniform(0.1, 0.5, numberSim)
    lambda0Sim = np.random.uniform(0, 200, numberSim)
    C_lambdaSim = np.random.uniform(0, 20, numberSim)
    delta_SSim = np.random.uniform(0, 0.1, numberSim)

    # Store and return simultations
    priorDf = np.column_stack((deltaSim, muSim, alphaSim, lambda0Sim, C_lambdaSim, delta_SSim))
    return pd.DataFrame(priorDf, columns=("delta", "mu", "alpha", "lambda_0", "C_lambda", "delta_S"))


class ParallelModel:
    @staticmethod
    def preisSim(param1, param2, param3, param4, param5, param6,
                 N_A, p_0, L, MCSteps):
        # initialize preis model object with specified parameters
        from preisSeed import PreisModel
        import pandas as pd

        p = PreisModel(N_A=N_A,
                       delta=param1,
                       mu=param2,
                       alpha=param3,
                       lambda_0=param4,
                       C_lambda=param5,
                       delta_S=param6,
                       p_0=p_0,
                       T=L,
                       MC=MCSteps)

        # Start model
        p.simRun()
        p.initialize()

        # Simulate price path for T=L time-steps
        p.simulate()

        return pd.DataFrame(p.intradayPrice)

    def f(self, x):
        return self.preisSim(*x)


if __name__ == '__main__':
    np.random.seed(19950412)

    # 1. s() = I()
    # 2. Compute s(y_obs) = y_obs
    L = 2300  # TimeHorizon
    N = 40    # number of simulated y_obs
    p_0 = 100 * 1000
    MCSteps = 10 ** 5
    N_A = 250

    param = priorSimulation(N)
    # randomise seed
    seed = int(np.random.uniform(0, 1000000))
    np.random.seed(seed)

    print(time.asctime())
    print(mp.cpu_count())
    para_model = ParallelModel()
    with mp.Pool((mp.cpu_count()-1)) as pool:
        results_list = pool.map(
            para_model.f, [(param.iloc[i, 0], param.iloc[i, 1], param.iloc[i, 2],
                            param.iloc[i, 3], param.iloc[i, 4], param.iloc[i, 5], N_A, p_0, L, MCSteps) for i in range(N)]
        )

    results_df = pd.concat(results_list, axis=1)
    results_df = results_df.div(1000)

    param = pd.DataFrame(param)
    param.to_csv('param.csv', index=False)
    results_df.to_csv('out.csv', index=False)
    print(time.asctime())
