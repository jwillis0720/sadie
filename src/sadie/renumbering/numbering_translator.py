from pathlib import Path


class NumberingTranslator:
    def __init__(self):
        # Pathing for Numbering HMMs
        self.hmm_folder = Path(__file__).parent / "data/anarci/HMMs"
        self.hmm_paths = list(self.hmm_folder.glob("*_[a-zA-Z].hmm"))
        # Numbering specific species and chains
        self.species = set()
        self.chains = set()
        self.species_chain_to_paths = self.get_species_chain_to_paths()

    def get_species_chain_to_paths(self) -> dict:
        """
        Return a dictionary with the available species and their paths.

        Returns
        -------
        dict
            Dictionary with the available species and their paths.
        """
        species_chain_to_paths = {}

        for path in self.hmm_paths:
            species, chain = path.stem.split("_")
            species = species.strip().lower().replace(" ", "_")
            chain = chain.strip().upper()
            try:
                species_chain_to_paths[(species, chain)].append(path)
            except KeyError:
                species_chain_to_paths[(species, chain)] = [path]

            self.species.add(species)
            self.chains.add(chain)

        return species_chain_to_paths
