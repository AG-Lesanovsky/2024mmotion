# 2024mmotion
The codes contained in this folder refer to the work "Impact of micromotion and field-axis misalignment on the excitation of Rydberg states of ions in a Paul trap", Wilson S. Martins, Joseph W. P. Wilkinson, Markus Hennrich, and Igor Lesanovsky.

All codes are written in Python 3.11.5, and the libraries required to run the codes are Matplotlib, NumPy and QuTip, in the following versions:

QuTiP Version: 4.7.3
Numpy Version: 1.25.2
Matplotlib Version: 3.7.2

All codes start with mm_, for micromotion. Codes in which data is generated are mm_quasienergy_.... In this code, the quasienergy operator is built up to a certain number of blocks, corresponding to a truncation of the operator, in principle, infinite. Thus, its eigenvectors and eigenvalues ​​are obtained, in addition to calculating the overlaps with the state |\downarrow>. The operating times are also calculated. Basically, the code is built so that there are five inputs: Nx = number of bosons, num_blocks = number of blocks in the matrix, z = frequency of the oscillatory field and delta = misalignment of the electric field null.

All the data in the work are already saved in the data folder but can be retrieved by simply changing the parameters in the script already mentioned.

The other scripts construct figures. The codes generating plots are mm_thermal_overlap_....

The figures folder contains .svg files, with vectorized images obtained using Inskscape, as well as figures obtained from scripts that are later processed using the same program.