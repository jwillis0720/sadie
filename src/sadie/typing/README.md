# Custom Typing for Sadie

```python
from sadie.typing import Source, Species, Chain
from pydantic import ValidationError, BaseModel, Field
from typing import List, Union

class MyClass(BaseModel):
    source: Source
    species: Union[List[Species], Species]
    chain: Chain
```
