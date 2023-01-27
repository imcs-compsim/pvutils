# Macros in ParaView

Add a macro in ParaView with the following code (path to pvutils must be set):
```python
# -*- coding: utf-8 -*-
import sys
sys.path.append("<path to pvutils>")
from macros import macro_xxx
macro_xxx()
```