import os
import numpy as np
from .simulate import simulate, nr_valid
from .distance import Tumor2DDistance 
from .log_transform import log_transform

__all__ = ["simulate", "nr_valid", "log_model"]


stored_data_db = os.path.join(os.path.dirname(__file__), "db1.db")
stored_data_db_2 = os.path.join(os.path.dirname(__file__), "db2.db")


def load_default():
    path = os.path.join(os.path.dirname(__file__), "default_samples.npz")
    raw_data = np.load(path)
    data_var = {key: val.var(axis=0) for key, val in raw_data.items()}
    data_mean = {key: val.mean(axis=0) for key, val in raw_data.items()}    

    for key in ["extra_cellular_matrix_profile", 'proliferation_profile']:
        n_valid_points = nr_valid(data_mean[key])
        nonzero = np.count_nonzero(data_mean[key])
        assert n_valid_points >= nonzero
        data_mean[key] = data_mean[key][:n_valid_points]
        data_var[key] = data_var[key][:n_valid_points]
        
    return raw_data, data_mean, data_var


def log_model(x):
    return log_transform(simulate)(**x)

    
    
distance = Tumor2DDistance(load_default()[2])


def animated_gif(db, output):
    from pyabc import History
    from pyabc.visualization import plot_kde_matrix
    import matplotlib.pyplot as plt
    import subprocess
    import tempfile
    tempdir = tempfile.mkdtemp()
    print("tmpdir", tempdir)
    import os
    h_loaded = History("sqlite:///" + db)

    limits = dict(log_division_rate=(-3, -1),
                  log_division_depth=(1, 3),
                  log_initial_spheroid_radius=(0, 1.2),
                  log_initial_quiescent_cell_fraction=(-5, 0),
                  log_ecm_production_rate=(-5, 0),
                  log_ecm_degradation_rate=(-5, 0),
                  log_ecm_division_threshold=(-5, 0))

    for t in range(h_loaded.n_populations):
        print("Plot population {t/h_loaded.n_populations}")
        df, w = h_loaded.get_distribution(m=0, t=t)
        plot_kde_matrix(df, w, limits=limits)
        plt.savefig(os.path.join(tempdir, f"{t:0>2}.png"))

    subprocess.run(["convert", "-delay", "15",
                    os.path.join(tempdir, "%02d.png"),
                    output + ".gif"])


