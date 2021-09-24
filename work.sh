make nano
bsub -J nano -q batch -R "span[ptile=36]" -n 144 -o test1.txt "mpirun ./nano" -fpetrap none  -mesh_file ./nano_mesh/nano_l1.mesh -bias_anode 5.0 -bias_surface_charge -1.0e-7 -analytic_density 1.0e-3
