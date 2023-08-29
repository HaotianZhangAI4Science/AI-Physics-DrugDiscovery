This repo contains four strategies collected by Odin, and used for training Delete model,  **D**eep **Le**ad Op**t**imization **E**nveloped in Protein Pocket through Unified Deleting Strategies and a Structure-aware Network. 

I provide the source code here, hope this collection could help you!

The detailed implementation could be found at `fragmentation.ipynb`

I have updated the important function to track the anchor node in the modified molecules after any decomposition methods. 

```python
from frag import remove_structures(mol, sub_mol)
```

