Novice.Original.Slope1 code has been associated to find optimal solution of the basis pursuit denoising problem with bounded variables as follows:
min ||Sv||+mu||v||_1
s.t. l_1<v<u_1,
 obtained by considering convex relaxation of the following compressed sensing problem:
 min_{v\in R^n} ||v||_0,
s.t. Sv=0,
l_1<=v<=u_1.
To run this program, you may run driver_cs file. To do so:    
please enter your file name, your sheet and your data range for matrix S on lines 11-13.
please enter  your data range for lower bound vector l_1 on the line 24.
please enter  your data range for lower bound vector u_1 on the line 31.

Input files must be in .xlsx (or .xls) format.
