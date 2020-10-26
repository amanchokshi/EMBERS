---
title: "EMBERS: Experimental Measurement of BEam Responses with Satellites"
tags:
  - Python
  - Radio astronomy
  - Satellites
  - Experiment
  - Telescope
authors:
  - name: A. Chokshi
    orcid: 0000-0003-1130-6390
    affiliation: "1, 3"
  - name: J. L. B. Line
    orcid: 0000-0002-9130-5920
    affiliation: "2, 3"
  - name: B. McKinley
    orcid: 0000-0002-9006-1450
    affiliation: "2, 3"
affiliations:
 - name: School of Physics, University of Melbourne, Parkville, Victoria, 3010, Australia
   index: 1
 - name: International Centre for Radio Astronomy Research, Curtin University, Perth, WA 6845, Australia
   index: 2
 - name: ARC Centre of Excellence for All Sky Astrophysics in 3 Dimensions (ASTRO 3D)
   index: 3
date: XXX
bibliography: paper.bib
---


# Introduction

*EMBERS* is a python package which provides a modular framework for radio telescopes and interferometric arrays such as the
MWA^[https://www.mwatelescope.org], HERA^[http://reionization.org], and the upcoming SKA^[https://www.skatelescope.org] to accurately measure the all sky
beam responses of their antennas using weather and communication satellites. Such measurements can reveal environmental factors which have perturbed ideal,
simulated beam shapes in complex ways. Results from recent work such as [@JLBLine_2018] reveal the presence of gradients in ground screens,
dead dipoles, and the effect of foliage near the antennas. Telescopes such as the MWA, HERA & SKA, which may benefit from in-situ beam measurements, are involved
in large scale surveys and the search for some of the earliest signals in our Universe. Such studies push the boundaries of precision calibration where undetermined
beam errors could potentially introduce spurious contaminants and hinder detections. *EMBERS* could form the backbone of a passive parallel monitoring system for large
radio telescopes, concurrently measuring beam shape without any disruption to regular observations, providing astronomers with additional information to include in their
instrumental models.


# Statement of Need

Large scale implementation of a satellite based beam measurement system have become necessary for the accurate calibration of cutting-edge and next generation radio telescopes. 
Previous efforts at such techniques have validated this method on a significantly smaller scale. The large data volume recorded over months of observations and
the absence of modular, parallelized software packages, necessitated the creation of an automated pipeline for such analyses. *EMBERS* contains modules to pre-process and temporally align raw
RF data, download large batches of satellite ephemerides from Space-Track.org^[https://www.space-track.org] and compute the trajectories of satellites using
Skyfield [@Skyfield_2019]. *EMBERS* implements a unique cross-matching technique to automatically determine the transmission frequency of satellites based on
their trajectories and observed RF power. Satellite signals are further processed to remove modulations due to the satellites and non-linear 
amplification effects, before being projected onto all-sky beam maps.

With plans for extremely large and sensitive interferometric arrays such as the SKA being proposed, the need for a in-situ beam calibration system is critical. EMBERS
presents a simple framework for the analysis of large volumes of satellite data, enabling in-situ beam measurements with ease. 


# Theory

At the heart of *EMBERS* is a simple concept. Satellite measurements are simultaneously made using multiple Antennas Under Test (AUTs) and reference antennas
(ref) which enables us to account for any modulation in transmitted satellite power. The power received by the AUT is the product of the beam response $B_{AUT}$
and the flux transmitted by the satellite $F$. A reference antenna with a simple, well known beam response $B_{ref}$ is used to record the modulation of the
transmitted flux, and can subsequently be used to compute the beam shape of the AUT. The power received by the AUT and reference antenna are $P_{AUT} = B_{AUT}F$ and $P_{ref} = B_{ref}F$ respectively.
These expressions can be combined to obtain the beam response of the AUT:


\begin{equation}
    B_{AUT} = \frac{P_{AUT}}{P_{ref}}B_{ref}.
	\label{eq:beam_eq}
\end{equation}


With each satellite pass, we measure a cross sectional slice of the AUT beam response. With sufficient observation time, an all-sky beam response is built up.

![MWA beam maps generated using *EMBERS* with data from (Chokshi et al, in prep). The first row ((i) - (iii)) represent all sky measured beam maps, while the second row ((iv) - (vi)) represent residuals between measurements and cutting edge simulations, with the gray regions denoting the nulls surrounding the primary and secondary beams.](https://raw.githubusercontent.com/amanchokshi/EMBERS/master/docs/_static/imgs/beam_maps_joss.jpg)


# Acknowledgements

I would like to thank Rachel Webster for introducing me to this project and supporting my work. Further, I would like to thank Nichole Barry for her valuable
suggestion for improving and expanding *EMBERS*. Parts of this research were supported by the Australian Research Council Centre of Excellence for All Sky
Astrophysics in 3 Dimensions (ASTRO 3D), through project number CE170100013.


# References
