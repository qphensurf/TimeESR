# TimeESR Examples

# 1_CW_ESR
A model simulation of a continous-wave ESR experiment, using a single spin 1/2 transport site and several DC measurements with several different driving frequencies which range between 16.9 to 17.1 GHz. The DC current vs driving frequency plot is then computed.

# 2_CW_ESR_RESOL_FREQ
Similar to the above example, covering a driving frequency range between 16.8 to 17.2 GHz, but for a range of differential frequency steps. This example shows the importance of choosing a fine enough differential frequency in order to resolve important features in the spectra.

# 3_CW_ESR_RESOL_TIME
Similar to the above example, but for a range of time points. This example shows the importance of choosing a fine enough time step.

# 4_ONE_SPIN
Driving a spin 1/2 transport site from the prepared |0> state, taken to be the ground state with an applied B field in the x direction with strength 0.5 T (i.e., a |down> state in the x basis), to a (|0> + |1>)/Sqrt[2] state, i.e., a Hadamard operation on a single spin 1/2 qubit. The calculated driving frequency is 13.9953 Ghz for this operation. A rotating frame of reference is used in the post-processing step to remove the inherent larmor frequency of the system, and the system slowly decoheres after reaching the final state.

# 5_TWO_SPINS
Driving a spin 1/2 transport site exchange coupled to another spin 1/2 with strength J = -5.0 GHz (isotropically applied) and an applied B field in the x direction with strength 0.5 T. The first four states are, in the x basis of both spins with the transport electron first, | down down >, (| up down > + | down up >)/Sqrt[2], (| up down > - | down up >)/Sqrt[2], and | up up >. The intent is to transition all populations in the | down down > state to the | up up > state. This is achieved by driving the system from the ground state to the first excited state at 19.0 GHz, then transitioning from the first excited state to the third excited state with a drive at 9.0 GHz. Again, the system slowly decoheres after reaching the final state.

# 6_FIVE_SPINS
A model simulation of a DC sweeped-voltage experiment for a spin 1/2 transport site coupled to a spin system of the following geometry: transport site - spin 1 - spin 1 - spin 1/2, with another spin 1/2 site exchange coupled to the first spin 1 site. This represents a radical porphyrin (Open Acess https://onlinelibrary.wiley.com/doi/full/10.1002/advs.202105906) except that there are two metallic centers (two Fe atoms, each with S=1). A driving frequency of 50.0 GHz is used, with bias values swept from -15.0 to 15.0 mV, and then the current vs voltage plot is computed.