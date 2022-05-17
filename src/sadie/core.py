from sadie.hmmer import HMMER
from sadie.airr import Airr


class API:
    
    hmmer = HMMER
    airr = Airr
    
    def __init__(self, alignment_algo='hmmer', *args, **kwargs):
        if alignment_algo.lower().strip() == 'hmmer':
            self.alignment_algo = self.hmmer(*args, **kwargs)
        elif alignment_algo.lower().strip() == 'airr':
            self.alignment_algo = self.airr(*args, **kwargs)
            
    def run(self, *args, **kwargs):
        return self.alignment_algo.run_single(*args, **kwargs)
        
    