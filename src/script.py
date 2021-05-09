import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pymc3 as pm
import pandas as pd

import aesara

from pymc3.ode import DifferentialEquation
from scipy.integrate import odeint


N_BURN = 1000
N_SAMPLE = 3000

def SIR(y, t, p):
    ds = -p[0] * y[0] * y[1]
    di = p[0] * y[0] * y[1] - p[1] * y[1]
    return [ds, di]


if __name__ == "__main__":
    data = pd.read_csv("../data/cases_by_status_and_phu.csv", index_col="PHU_NAME")
    toronto_data = data.loc["TORONTO"]
    del data

    times = range(0, len(toronto_data))
    N = 2.3 * 1000000
    I = toronto_data["ACTIVE_CASES"].values
    R = toronto_data["RESOLVED_CASES"].values + toronto_data["DEATHS"].values
    S = N - I - R

    sir_model = DifferentialEquation(
            func=SIR,
            times=times,
            n_states=2,
            n_theta=2
            )

    yobs = np.hstack((S.reshape(-1, 1), I.reshape(-1, 1)))
    yobs = np.log(yobs)

    with pm.Model() as model:
        sigma = pm.HalfCauchy("sigma", 1, shape=2)
                
        R0 = pm.HalfNormal("R0", 3)
        lam = pm.Lognormal("lambda", pm.math.log(2), 2)
        beta = pm.Deterministic("beta", lam*R0)

        sir_curves = sir_model(y0=yobs[0], theta=[beta, lam])

        Y = pm.Normal("Y", mu=sir_curves, sigma=sigma, observed=yobs)

        trace = pm.sample(N_SAMPLE, tune=N_BURN, target_accept=0.9, cores=1)
        pm.save_trace(trace, "sir.trace")
