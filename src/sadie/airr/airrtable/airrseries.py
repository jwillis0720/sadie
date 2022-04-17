import pandas as pd

from sadie.airr.models import AirrSeriesModel


class AirrSeries(pd.Series):
    _metadata = ["meta"]  # add custom namespaces here

    def __init__(self, data, *args, **kwargs):
        super(AirrSeries, self).__init__(data=data, *args, **kwargs)
        if not isinstance(data, pd.core.internals.managers.SingleBlockManager):
            if isinstance(data, pd.core.series.Series):
                self._verify()

    @property
    def _constructor(self) -> "AirrSeries":  # type: ignore:
        return AirrSeries

    def _verify(self) -> None:
        data = AirrSeriesModel(**self).dict()
        self.update(data)
