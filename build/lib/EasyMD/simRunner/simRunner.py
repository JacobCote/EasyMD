from EasyMD.runners.gbisRunner import GBISRunner
from EasyMD.runners.simAnnealingRunner import AnnealingRunner
from EasyMD.runners.solvatedRunner import SolvatedRunner

class SimRunner():
    def __init__(self,config,modeller,system):
        self.config = config
        self.sim = self._setupSim()

    def _setupSim(self,):
        if self.config.simulated_annealing:
            return AnnealingRunner(self.config)
        elif self.config.gbis : 
            return GBISRunner(self.config)
        else :
            return SolvatedRunner()
        

