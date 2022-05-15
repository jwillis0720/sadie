# Custom Typing for Sadie

```python
from sadie.typing import Source, Species, Chain
from pydantic import ValidationError, validate_arguments

class MyClass:

    @validate_arguments
    def test(source: Source, species: Union[List[Species], Species], chain: Chain):
        return source, species, chain
```
