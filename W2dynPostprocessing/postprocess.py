import W2dynPostprocessing.hdf5 as hdf5
import numpy as np
import W2dynPostprocessing.symmetries as sym


class Calculation():

    def __init__(self, filename):
        self.w2dfile = hdf5.W2dynfile(filename)
        self.beta = self.w2dfile.get_beta()

    def get_chi_pp(self):
        if self.w2dfile.check_for_statiter():
            print("found statiters")
        pass

    def get_chi_ph(self):
        pass
