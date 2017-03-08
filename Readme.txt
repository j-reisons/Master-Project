This repository is used to store/backup codes and results of my master's project.
The subject of the master's project is the study of open quantum systems,
and the main numerical method used is DMRG/Matrix Product States.

General functions are the core of the MPS algorithms : braket of two MPS, 
application of MPO to MPS, canonization, compression, etc...

MPS are stored as cells of 3d arrays. The order of the indices is
left bond index, right bond index, physical index.

/Models/Heisenberg contains functions for generating MPO evolution 
operators for the closed, open (Trotter decomposition of order 2 and 4),
and open disordered Heisenberg spin chain. 
/Models/Heisenberg also contains the run_trajectories_mps function which is a time evolution algorithm.