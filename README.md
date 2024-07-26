# 2024mmotion
The codes contained in this folder refer to the work Micromotion Effects on Rydberg Excitation of Trapped Ions, Wilson S. Martins, Joseph W. P. Wilkinson, and Igor Lesanovsky.

All codes are written in Python 3.11.5, and the libraries required to run the codes are Matplotlib, NumPy and QuTip, in the following versions:

QuTiP Version: 4.7.3
Numpy Version: 1.25.2
Matplotlib Version: 3.7.2

All codes start with mm_, for micromotion, and the only code that generates data is mm_quasienergy_op_funits_sm.py. In this code, the quasienergy operator is built up to a certain number of blocks, corresponding to a truncation of the operator, in principle, infinite. Thus, its eigenvectors and eigenvalues ​​are obtained, in addition to calculating the overlaps with the state |nS_1/2>. The operating times are also calculated. Basically, the code is built so that there are five inputs: Nx = number of bosons, num_blocks = number of blocks in the matrix, eta = frequency of the oscillatory field and delta = misalignment of the electric field cancellations.

All the data in the work are already saved in the data folder but can be retrieved by simply changing the parameters in the script already mentioned.

The other scripts construct figures or functions obtained by perturbation theory. The former ones concern mm_quasienergy_op_vardx.py, mm_quasienergy_op_vareta.py, mm_overlap_proj_vardx.py, mm_overlap_proj_vareta.py, while the latter one is contained in the code mm_overlap_pert.py.

The figures folder contains .svg files, with vectorized images obtained using Inskscape, as well as figures obtained from scripts that are later processed using the same program.