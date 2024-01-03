# BNET: Neural Network-Based Wireless Communication Framework

## Overview
"BNET: A Neural Network Approach for LLR-based Detection in the Presence of Bursty Impulsive Noise" introduces a neural network-based framework designed to enhance the reliability of wireless communications in challenging noise environments. The BNET framework approximates the log-likelihood ratio (LLR) function of the BCJR algorithm, offering a more computationally efficient approach while maintaining, or in some cases surpassing, the Bit Error Rate (BER) performance of conventional methods.

## Features
- **Neural Network-Based LLR Approximation**: Utilizes a neural network to approximate the LLR function, providing efficiency in computation.
- **Adaptability to Noise**: Performs effectively in environments with bursty impulsive noise, particularly suited for single-carrier communication systems.
- **Dynamic Performance**: Capable of adapting to changing noise model parameters over time.

## Contents
- `BER_tester.m`: Tests the BER for given LLR values using the BNET framework or the BCJR algorithm.
- `Testing BER for Conventional methods.m`: Compares BER performance of conventional methods in different scenarios.
- `generate_test_data.m` and `generate_train_data.m`: Scripts to generate test and training data for the BNET model. This data includes LLRs from MAP-BCJR as targets and received signals as inputs.
- `MAPdecoding_TSMGaussianNoise_Real_BPSK.m`: Implements MAP-BCJR decoding in the presence of Two-State Markov Gaussian Noise for BPSK modulation. (Empty file for IP/NDA agreement purposes)
- `BNET Training and Plotting BER.ipynb` and `BNET-Transformer Training and Plotting BER.ipynb`: Jupyter notebooks for training BMET models and testing their performance.
- `TSMG.m`: Generates samples of Two-State Markov Gaussian noise.

## Getting Started
1. **Prerequisites**: Ensure you have MATLAB installed, as most scripts are MATLAB-based. For Jupyter notebooks, an environment supporting Jupyter (like Anaconda) is required.
2. **Installation**: Clone the repository using `git clone [URL]`.
3. **Running the Scripts**: Each script can be run individually in MATLAB. For Jupyter notebooks, open them in your Jupyter environment and run the cells sequentially.

## Contributing
Contributions to this project are welcome. Please fork the repository and submit a pull request with your proposed changes.

## Citation
If you use this work in your research, please cite:

```
@ARTICLE{BNET,
author={Barka, Hazem and Alam, Md Sahabul and Kaddoum, Georges and Sacuto, Fabien and Agba, Basile L.},
journal={IEEE Wireless Communications Letters},
title={{BNET: A Neural Network Approach for LLR-Based Detection in the Presence of Bursty Impulsive Noise}},
year={2023},
volume={12},
number={1},
pages={80-84},
doi={10.1109/LWC.2022.3217675}}
```

## License
This project is licensed under the MIT License.
