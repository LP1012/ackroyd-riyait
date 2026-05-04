import meshio

spatial_mesh = meshio.read("../gmsh/square-source-coarse.msh")
angle_mesh = meshio.read("../gmsh/solid-angle.msh")

print(spatial_mesh)
