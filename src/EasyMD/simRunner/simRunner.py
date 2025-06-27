from EasyMD.runners import SolvatedRunner,AnnealingRunner,GBISRunner,RestartRunner

class SimRunner():
    def __init__(self,config,modeller,system):
        self.config = config
        self.sim = self._setupSim()

    def _setupSim(self,):
        if self.config.restart:
            return RestartRunner
        if self.config.simulated_annealing:
            return AnnealingRunner(self.config)
        elif self.config.gbis : 
            return GBISRunner(self.config)
        
        else :
            return SolvatedRunner()
        

