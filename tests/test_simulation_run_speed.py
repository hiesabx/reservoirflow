import pandas as pd
import numpy as np
from openresim import fluids, grids, models
import matplotlib.pyplot as plt
import seaborn as sns

generate_data = False
show_result = False

if generate_data:

    def create_model(n):
        grid = grids.Cartesian(
            nx=n,
            ny=1,
            nz=1,
            dx=400,
            dy=500,
            dz=50,
            phi=0.21,
            kx=273,
            comp=0,
            dtype="double",
            unify=False,
        )
        fluid = fluids.SinglePhase(
            mu=1.5, B=1, rho=50, comp=2.5 * 10**-5, dtype="double"
        )
        model = models.Model(grid, fluid, pi=3000, dt=5, dtype="double")
        id = grid.cells_id[0]
        model.set_well(id=id, q=-400, pwf=1500, s=0, r=3)
        # model.set_boundaries({0: ("rate", 0)})
        return model

    results = {
        "not.vectorized": [],
        "vectorized": [],
        "num.steps": [],
        "num.grids": [],
    }

    for s in [10, 30]:
        for n in [5, 10, 100]:
            model = create_model(n)
            model.run(nsteps=s, vectorize=False)
            results["not.vectorized"].append(model.ctime)
            model = create_model(n)
            model.run(nsteps=s, vectorize=True)
            results["vectorized"].append(model.ctime)
            results["num.grids"].append(n)
            results["num.steps"].append(s)

    df = pd.DataFrame(results)
    df.to_csv("tests/test_simulation_run_speed.csv", index=False)

if show_result:
    df = pd.read_csv("tests/test_simulation_run_speed.csv")
    df1 = df[["vectorized", "num.grids", "num.steps"]]
    df1.rename(columns={"vectorized": "ctime"}, inplace=True)
    df1["Type"] = "Vectorized"
    df2 = df[["not.vectorized", "num.grids", "num.steps"]]
    df2.rename(columns={"not.vectorized": "ctime"}, inplace=True)
    df2["Type"] = "Not Vectorized"
    df = pd.concat([df1, df2], axis=0)
    print(df)
    df["num.steps"] = df["num.steps"].astype(str)
    df["num.grids"] = df["num.grids"].astype(str)
    df["steps|grids"] = df["num.steps"].astype(str) + "|" + df["num.grids"].astype(str)
    # plt.bar(df["steps|grids"], height=df["not.vectorized"], label="Not vectorized")
    # plt.bar(df["steps|grids"], height=df["vectorized"], label="Vectorized")
    g = sns.catplot(
        col="num.steps",
        y="ctime",
        hue="Type",
        x="num.grids",
        data=df,
        kind="point",
    )
    g.set_axis_labels("Number of Grids", "Total Simulation Run Time [seconds]")
    # plt.legend()
    plt.show()
