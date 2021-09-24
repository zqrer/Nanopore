make nano
mpirun -np 32 ./nano -fpetrap none  -mesh_file ./nano_mesh/nano_l1.mesh -bias_anode 5.0 -bias_surface_charge -1.0e-7 -analytic_density 1.0e-3
