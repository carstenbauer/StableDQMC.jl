---
title: 'StableDQMC.jl: Numerical stabilization routines for determinant quantum Monte Carlo'
tags:
  - Julia
  - physics
  - condensed matter
  - quantum monte carlo
authors:
  - name: Carsten Bauer
    orcid: 0000-0001-6719-0978
    affiliation: 1
affiliations:
 - name: Institute for Theoretical Physics, University of Cologne, 50937 Cologne, Germany
   index: 1
date: 27 February 2020
bibliography: paper.bib
---

# Introduction

At the heart of condensed matter physics is the desire to gain a solid understanding of complex interacting many-particle systems, such as topological insulators, unconventional metals, and magnetic materials. Quantum Monte Carlo simulations are a powerful tool for making significant progress in this direction as they, in contrast to analytical approaches, allow for a numerically exact and unbiased study of these materials. In particular, the determinant quantum Monte Carlo (DQMC) method [@Blankenbecler1981] has proven to be effective for the analysis of itinerant electron systems like iron- or copper-based high-temperature superconductors [@Schattner2016; @Assaad2016; @Bauer2020]. However, its application requires careful handling of numerical instabilities which inevitably occur because of an exponentially spread of numerical scales due to the intrinsic energy spectra of such compounds.

Although  It is a Julia [@Bezanson2017] library.

# Acknowledgements

We thank Peter Bröcker, Frederick Freyer, Snir Gazit, Yoni Schattner, and Simon Trebst for useful discussions during the genesis of this project.

# References