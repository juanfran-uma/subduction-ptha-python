# Example external mesh

This folder is **optional**. A normal run does not use it: `generate.py` builds the mesh from the SLAB model by itself.

It holds one example of the file you would give to `--mesh-file` if you already have your own quadrilateral mesh: `alaskaaleutians_quadrilateral_coors.dat`, one line per unit source with the four corners of the cell as `lon lat depth`. The same format is written by every run in `mesh/mesh.dat`.

```
python ptha_python/generate.py alaska --zone alaskaaleutians --folder alaskaaleutians --mesh-file inputs_meshes/alaskaaleutians_quadrilateral_coors.dat
```

The format and its rules are in `ptha_python/html/docs/mesh_file.html`. `ptha_python/validation/validate_mesh_file.py` reads this file too.
