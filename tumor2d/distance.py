from .simulate import nr_valid
import numpy as np

class Tumor2DDistance:
    __name__ 
    def __init__(self, variances: dict):
        self.variances = {key: val[:nr_valid(val)]
                          for key, val in variances.items()}
        self.inv_variances = {}
        for key, val in self.variances.items():
            inv = np.zeros(len(val))
            inv[val != 0] = 1 / val[val != 0]
            self.inv_variances[key] = inv
    
    def initialize(self, sample_from_prior):
        pass
    
    def __call__(self, x: dict, y: dict) -> float:
        length = {key: min([len(x[key]),
                            len(self.variances[key]),
                            len(y[key])]) for key in y}
        return sum(np.sum((x[key][:length[key]] - y[key][:length[key]])**2
                              * self.inv_variances[key][:length[key]])
                   for key in length.keys())
    
    def get_config(self):
        return {}
