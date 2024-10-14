# Spatial Mapping

This repository contains the code used to produce the figures in X. 

<div align="center">
  
![Rotating fundamental QNM](figs/mapping_animation.gif)

*Figure 1: An animation of the reconstructed fundamental QNM at a range of t0 values. The derivative of the mismatch with respect to the start time is related to the real part of the QNM frequency (the blue dashed line).*

</div>

# Requirements & usage 

To use this code you will need to install [qnmfits](https://github.com/eliotfinch/qnmfits). 

This code uses NR waveforms from the Spectral Einstein Code. CCE waveforms can be obtained from the [SXS Gravitational Waveform Database](https://data.black-holes.org/waveforms/extcce_catalog.html). All simulations were transformed into the superrest frame using the [`scri`](https://pypi.org/project/scri/) package. 

Waveforms are imported into the workbooks using the function `CCE.SXS_CCE()`. This imports a `qnmfits.Custom()` object containing the transformed CCE data. This function is not included in this repo and will need to be replaced with a suitable method for importing the waveform data.  