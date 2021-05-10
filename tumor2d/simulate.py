from ._tumorutil import tumor2d_simulate
import numpy as np
from multiprocessing import Process, Pipe
MAX_SEED = 2147483647

def nr_valid(arr):
    return len(arr) - len(np.nonzero((arr[::-1] == 0).cumprod())[0])


def simulate(division_rate=4.17e-2,
             initial_spheroid_radius=1.2e1,
             initial_quiescent_cell_fraction=7.5e-1,
             division_depth=100,
             ecm_production_rate=5e-3,
             ecm_degradation_rate=8e-4,
             ecm_division_threshold=1e-2,
             randseed=None):
    """

    Parameters
    ----------


    .. notes::

        Jagiella, Nick, Dennis Rickert, Fabian J. Theis, and Jan Hasenauer.
        “Parallelization and High-Performance Computing Enables Automated
        Statistical Inference of Multi-Scale Models.”
        Cell Systems 0, no. 0 (2017). doi:10.1016/j.cels.2016.12.002.

        First table on page e4

        The zeros are not missing data.
        They are actual biology.
        Do not discard them.


    Returns
    -------

    The return vectors (can) contain lots of zeros at the end.
    These can be ignored.

    For

      * growth_curve: the array entries correspond to time
       * extra_cellular_matrix and proliferation: array entries correspond
         to radial distance, measured at the time point profiletime

    """
    if randseed is None:
        randseed = np.random.randint(MAX_SEED)

    kwargs = dict(max_celldivision_rate=division_rate,
                  initial_radius=initial_spheroid_radius,
                  initial_quiescent_fraction=initial_quiescent_cell_fraction,
                  division_depth=division_depth,
                  emc_productionrate=ecm_production_rate,
                  ecm_degradationrate=ecm_degradation_rate,
                  ecm_threshold_quiescence=ecm_division_threshold,
                  randseed=randseed)

    parent, child = Pipe()

    def pipe_call(pipe, kwargs):
        res = tumor2d_simulate(**kwargs)
        pipe.send(res)

    proc = Process(target=pipe_call, args=(child, kwargs), daemon=True)
    proc.start()
    res = parent.recv()
    proc.join()


    # convert tu numpy
    growth_curve = np.array(res.growth_curve)
    extra_cellular_matrix_profile = np.array(res.extra_cellular_matrix)
    proliferation_profile = np.array(res.proliferation)

    return dict(growth_curve=growth_curve,
                extra_cellular_matrix_profile=extra_cellular_matrix_profile,
                proliferation_profile=proliferation_profile)
