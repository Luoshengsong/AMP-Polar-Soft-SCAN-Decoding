# AMP-Polar_Soft_Decoding

We consider an uplink wireless channel MIMO model where there are 

- $M$ antennas at BS
- $K$ single antenna terminals

Thus, we denote the effective channel as $\boldsymbol H \in \mathbb C^{M \times K}$, which is Rayleigh distributed. This system model is characterized as 
$$\boldsymbol Y = \boldsymbol H \boldsymbol X + \boldsymbol W$$
where the element of $\boldsymbol X$ is BPSK modulated, and the element of $\boldsymbol W$ is AWGN.

Assume perfect CSI is known in BS side, we construct an iterative Bayesian receiver based on Turbo structure:

[**Module-1**: MIMO Detector] <-- $\Pi$ --> [**Module-2**: Polar SCAN decoding]

The simulation result:

![image](https://github.com/Luoshengsong/AMP-Polar_Soft_SCAN_Decoding/assets/73685146/f69065b3-579e-4342-b585-a53f97fdf5da)

## Reference
[1] U. U. Fayyaz and J. R. Barry, "Polar codes for partial response channels," 2013 (ICC).
[2] Z. Pan, E. Li, et,al., "Joint Iterative Detection and Decoding Receiver for Polar Coded SCMA System," 2018 (ICC Workshops).


