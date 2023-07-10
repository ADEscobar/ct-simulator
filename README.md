# ct-simulator

The **ct-simulator** is a MATLAB script used to characterize the Packet Error Rate (PER) of radio transmissions using a Frequency-Shift Keying (FSK) modulation scheme when *two* concurrent transmitters (ct) overlap in the air sending *the same* bitstream in the presence of Additive White Gaussian Noise (AWGN). 

It provides practical insights to understand the performance of flooding protocols based on concurrent transmissions.


## Features

- Quick simulations based on a modified baseband-equivalent model tweaked to introduced the amplitude distortion of the envelope caused by the overlap of the transmissions
- Explained step-by-step process, including comparisons with theoretical results to validate the correctness of the simulations
- Obtention of PER vs Eb/N0 figures of merit
- Configurable FSK modulation parameters, such as the modulation order, sampling rate and frequency separation
- Easily extendable to different modulations and noise models
- Configurable packet lengths
- The relative energy of the concurrent transmissions can be modified, e.g. one transmitter can be set as dominant


## Run

Execute the `ct_simulator.m` script to trigger the execution of the sample simulations and the generation of the different figures that show the results. The markup of the script supports MATLAB publishing.

Tested using MATLAB R2022b.


## References

The **ct-simulator** has been used as a base for results presented in the following peer-reviewed articles:

- Escobar-Molero, A., 2019. Improving reliability and latency of wireless sensor networks using concurrent transmissions. *at-Automatisierungstechnik*, 67(1), pp.42-50.
- Baddeley, M., Boano, C.A., Escobar-Molero, A., Liu, Y., Ma, X., Raza, U., Römer, K., Schuß, M. and Stanoev, A., 2020. The impact of the physical layer on the performance of concurrent transmissions. In *2020 IEEE 28th International Conference on Network Protocols (ICNP)* (pp. 1-12). IEEE.
- Nahas, B.A., Escobar-Molero, A., Klaue, J., Duquennoy, S. and Landsiedel, O., 2021. BlueFlood: Concurrent Transmissions for Multi-hop Bluetooth 5—Modeling and Evaluation. *ACM Transactions on Internet of Things*, 2(4), pp.1-30.
- Baddeley, M., Boano, C. A., Escobar-Molero, A., Liu, Y., Ma, X., Marot, V., ... & Stanoev, A., 2023. Understanding Concurrent Transmissions. ACM Transactions on Sensor Networks.

And also in my dissertation:

- Escobar-Molero, A., 2020. Using concurrent transmissions to improve the reliability and latency of low-power wireless mesh networks (No. RWTH-2020-05541). Lehrstuhl für Integrierte Analogschaltungen und Institut für Halbleitertechnik.
