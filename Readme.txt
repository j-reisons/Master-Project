This is my attempt to keep track of the mess that has become my MPS code.

General functions are the core of the algorithms : braket of MPS, 
application of MPO to MPS, canonization, compression, etc...

/Models/Heisenberg contains functions for generating MPO evolution 
operators for the closed, open (to order 2 and 4), and open disordered 
Heisenberg spin chain. /Models/Heisenberg also contrains the 
run_trajectories_mps function which is a time evolution algorithm.

Test folders contain tests that were performed as we were programming.