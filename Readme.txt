This repository is used to store/backup codes and results of my master's project.

General functions are the core of the algorithms : braket of MPS, 
application of MPO to MPS, canonization, compression, etc...

/Models/Heisenberg contains functions for generating MPO evolution 
operators for the closed, open (Trotter decomposition of order 2 and 4),
and open disordered Heisenberg spin chain. 
/Models/Heisenberg also contains the run_trajectories_mps function which is a time evolution algorithm.