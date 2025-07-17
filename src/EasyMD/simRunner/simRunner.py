from EasyMD.runners import SolvatedRunner,AnnealingRunner,GBISRunner,Restarter

class SimRunner():
    def __init__(self,config,modeller,system,setup=None):
        self.config = config
        self.modeller = modeller
        self.system = system
        self.setup = setup
        self.sim = self._setupSim()
        

    def _setupSim(self,):
        if self.config.restart != None:
            return Restarter(self.setup,self.config,modeller=self.modeller,system=self.system)
        if self.config.simulated_annealing:
            return AnnealingRunner(self.config, modeller=self.modeller,system=self.system)
        elif self.config.GBIS : 
            return GBISRunner(self.config, modeller=self.modeller,system=self.system)
        
        else :
            return SolvatedRunner(self.config, modeller=self.modeller,system=self.system)
        
    def run(self,):
        self.sim.run()
        

