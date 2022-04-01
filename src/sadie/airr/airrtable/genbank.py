from typing import Dict, List, Union
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import datetime
from Bio.SeqFeature import FeatureLocation, SeqFeature


# Qualifier Dictionary
example_qualifiers_dict = {
    "gene": "gene",
    "latin": "latin",
    "organism": "species",
    "functional": "functional",
    "scheme": "IMGT",
}

# Example Annotations
annotations = {
    "organism": "latin",
    "source": "species",
    "date": (datetime.date.today().strftime("%d-%b-%Y")),
}

feature_types = [
    "FWR1",
    "FWR2",
    "FWR3",
    "FWR4",
    "CDR1",
    "CDR2",
    "CDR3",
    "VGene",
    "JGene",
    "DGene",
    "IGK",
    "IGK_Non-Productive",
    "IGH",
    "IGL",
]


class GenBankFeature:
    def __init__(
        self,
        start: int,
        end: int,
        feature_type: str,
        id: Union[str, None] = None,
        qualifier_dict: Union[Dict[str, str], None] = None,
    ):
        self.start = start
        self.end = end
        self.id = id

        # what type of feature is this
        self.feature_type = feature_type

        # Can have other info about our feature
        self.qualifier_dict = qualifier_dict

        # Feature Location
        self.location = FeatureLocation(self.start, self.end)

        # our main feature
        self._feature = SeqFeature(self.location, type=self.feature_type, qualifiers=self.qualifier_dict)

    @property
    def feature_type(self) -> str:
        return self._feature_type

    @feature_type.setter
    def feature_type(self, t: str) -> None:
        if t not in feature_types:
            raise TypeError(f"{t} must be in {feature_types}")
        else:
            self._feature_type = t

    @property
    def feature(self) -> SeqFeature:
        return self._feature


class GenBank:
    def __init__(
        self, sequence: Union[str, Seq], id: str, name: Union[str, None] = None, description: Union[str, None] = None
    ):
        self.sequence = sequence
        self.id = id
        if name:
            self.name = name
        else:
            self.name = id[0:16]
        self.description = description

        # Our main GB record
        self._record = SeqRecord(self.sequence, id=self.id, name=self.name, description=self.description)

    @property
    def record(self) -> SeqRecord:
        return self._record

    @property
    def features(self) -> List[SeqFeature]:
        _a: List[SeqFeature] = self.record.features
        return _a

    @property
    def sequence(self) -> Seq:
        return self._sequence

    @sequence.setter
    def sequence(self, seq: Union[str, Seq]) -> None:
        if isinstance(seq, str):
            self._sequence = Seq(seq)
        elif isinstance(seq, Seq):
            self._sequence = seq
        else:
            raise TypeError(f"{type(str)} must be instance of str or Bio.Seq")

    def add_feature(self, feature: GenBankFeature) -> None:
        if not isinstance(feature, GenBankFeature):
            raise TypeError(f"{feature} must be of type {GenBankFeature}")
        else:
            self.features.append(feature.feature)
