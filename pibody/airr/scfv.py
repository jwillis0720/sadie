# class AirrPhagemidHeavy(AirrChain):
#     def __init__(self, airr_entry):
#         super().__init__(airr_entry)
#         self.locus = self.get_locus()
#         if self.locus != "IGH":
#             raise TypeError("Airr entry not locus IGH %s", self.locus)


# class AirrPhagemidKappa(AirrChain):
#     def __init__(self, airr_entry):
#         super().__init__(airr_entry)
#         self.locus = self.get_locus()
#         if self.locus != "IGK":
#             raise TypeError("Airr entry not locus IGK %s", self.locus)


# class AirrPhagemidLambda(AirrChain):
#     def __init__(self, airr_entry):
#         super().__init__(airr_entry)
#         self.locus = self.get_locus()
#         if self.locus != "IGL":
#             raise TypeError("Airr entry not locus IGL %s", self.locus)


# class AirrPhagemidPair:
#     """Class for handling phagemid pairs that are both airr entries

#     This is perhaps the highest level abstraction you might want to deal with

#     """

#     def __init__(self, heavy_chain, light_chain):
#         """Constructure for AirrPhagemidPair

#         Args:
#             heavy_chain (AirrPhagemidHeavy): AirrphagemidHeavy Object
#             light_chain (AirrPhagemidKappa or AirrPhagemidLambda):AirrphagemidLight

#         Raises:
#             HeavyChainException: AirrphagemidHeavy is not right type
#             LightChainException: AirrphagemidLight is not right type
#         """
#         if not isinstance(heavy_chain, AirrPhagemidHeavy):
#             raise HeavyChainException("{} not heavy chain".format(heavy_chain))
#         if not isinstance(light_chain, (AirrPhagemidKappa, AirrPhagemidLambda)):
#             raise LightChainException("{} not lambda or Kappa chain".format(light_chain))
#         self.heavy_chain = heavy_chain
#         self.light_chain = light_chain

#     def get_json(self):
#         heavy_dict = json.loads(self.heavy_chain.get_json())
#         light_dict = json.loads(self.light_chain.get_json())
#         dictionary = {"HEAVY": heavy_dict, "LIGHT": light_dict}
#         return json.dumps(dictionary)

#     def to_json(self, file):
#         dictionary = json.loads(self.get_json())
#         json.dump(dictionary, open(file, "w"))

#     @staticmethod
#     def from_json(file):
#         """Read JSON from file"""
#         dictionary = json.load(open(file))
#         heavy = dictionary["HEAVY"]
#         light = dictionary["LIGHT"]
#         airr_chain_heavy = AirrChain.read_json(json.dumps(heavy))
#         airr_chain_light = AirrChain.read_json(json.dumps(light))
#         if airr_chain_light.get_locus() == "IGK":
#             light_obj = AirrPhagemidKappa
#         elif airr_chain_light.get_locus() == "IGL":
#             light_obj = AirrPhagemidLambda
#         else:
#             raise TypeError(f"{airr_chain_light.get_locus()} undefined for {file}")

#         # Use entry to reconstruc
#         return AirrPhagemidPair(AirrPhagemidHeavy(airr_chain_heavy.entry), light_obj(airr_chain_light.entry))
