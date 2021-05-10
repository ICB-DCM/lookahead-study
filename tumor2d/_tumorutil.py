import numpy as np
import logging

max_seed = 2147483647

#logging.basicConfig(level="DEBUG")
logger = logging.getLogger("TUMOR2D")
#logger.setLevel("ERROR")


class Tumor2dExperiment:
    def __init__(self, mean_gc, mean_ecm, mean_prolif, std_gc, std_ecm, std_prolif, full_data_gc=None,
                 full_data_ecm=None, full_data_prolif=None):
        self.full_data_gc = full_data_gc
        self.full_data_ecm = full_data_ecm
        self.full_data_prolif = full_data_prolif
        self.mean_gc = mean_gc
        self.mean_ecm = mean_ecm
        self.mean_prolif = mean_prolif
        self.std_gc = std_gc
        self.std_ecm = std_ecm
        self.std_prolif = std_prolif

    def compare_with_simulation(self, sim):
        score_gc = 0
        score_ecm = 0
        score_prolif = 0
        return 0


class Tumor2dSimulation:
    def __init__(self, val_gc, val_ecm, val_prolif):
        self.growth_curve = val_gc
        self.extra_cellular_matrix = val_ecm
        self.proliferation = val_prolif


def tumor2d_simulate(initial_radius=12.0, initial_quiescent_fraction=0.75,  max_celldivision_rate=0.0417,
                     division_depth=100, ecm_threshold_quiescence=0.010, emc_productionrate=0.005,
                     ecm_degradationrate=0.0008, endtime=500, outputrate=24, profiletime=408, profiledepth=1000,
                     randseed=None):
    """
    Tumor2d simulation.
    *Not* according to the published paper.

    Parameters
    ----------

    """
    # don't load at module level due to memory leaks in the original code
    import tumor2d.src.nixTumor2d as nixTumor2d
    if randseed is None:
        raise Exception("Randseed necessary.")
    profiletime /= 24
    pars = str(locals())
    logger.debug(f"START:{pars}")
    growth_curve, ecm_prof, prolif_prof = nixTumor2d.tumor2d_interface(initial_radius, initial_quiescent_fraction,
                                                                       max_celldivision_rate, division_depth,
                                                                       ecm_threshold_quiescence,
                                                                       emc_productionrate, ecm_degradationrate,
                                                                       endtime, outputrate, profiletime, profiledepth,
                                                                       randseed)
    logger.debug(f"DONE:{pars}")
    result = Tumor2dSimulation(growth_curve, ecm_prof, prolif_prof)
    return result


def tumor2d_statistic(num_reps=10, initial_radius=12.0, initial_quiescent_fraction=0.75,  max_celldivision_rate=0.0417,
                      division_depth=100, ecm_threshold_quiescence=0.010, emc_productionrate=0.005,
                      ecm_degradationrate=0.0008, endtime=500, outputrate=24, profiletime=408, profiledepth=1000,
                      randseed=np.nan):
    if np.isnan(randseed):
        randseed = np.random.randint(max_seed)
    np.random.seed(randseed)
    full_data_gc = np.empty([num_reps, np.int(np.floor(endtime / outputrate))])
    full_data_ecm = np.empty([num_reps, profiledepth])
    full_data_prolif = np.empty([num_reps, profiledepth])
    for i in range(num_reps):
        seed_simu = np.random.randint(max_seed)
        simu = tumor2d_simulate(initial_radius=initial_radius, initial_quiescent_fraction=initial_quiescent_fraction,
                                max_celldivision_rate=max_celldivision_rate,
                                division_depth=division_depth,
                                ecm_threshold_quiescence=ecm_threshold_quiescence,
                                emc_productionrate=emc_productionrate, ecm_degradationrate=ecm_degradationrate,
                                endtime=endtime, outputrate=outputrate, profiletime=profiletime,
                                profiledepth=profiledepth, randseed=seed_simu)
        full_data_gc[i] = simu.growth_curve
        full_data_ecm[i] = simu.extra_cellular_matrix
        full_data_prolif[i] = simu.proliferation
    mean_gc = np.mean(full_data_gc, axis=0)
    mean_ecm = np.mean(full_data_ecm, axis=0)
    mean_prolif = np.mean(full_data_prolif, axis=0)
    std_gc = np.max(np.std(full_data_gc, axis=0))
    std_ecm = np.max(np.std(full_data_ecm, axis=0))
    std_prolif = np.max(np.std(full_data_prolif, axis=0))
    result = Tumor2dExperiment(mean_gc, mean_ecm, mean_prolif, std_gc, std_ecm, std_prolif, full_data_gc, full_data_ecm,
                               full_data_prolif)
    return result

