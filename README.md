# POC-CMB-Harmonic

[![DOI](https://img.shields.io/badge/OSF%20DOI-10.17605/OSF.IO/VC59T-blue)](https://osf.io/vc59t) [![License: CC0](https://img.shields.io/badge/License-CC0%201.0-brightgreen)](https://creativecommons.org/publicdomain/zero/1.0/) 

Open-source project probing the **Principle of Cosmic Optimization (POC)**, seeking harmonic multifractal signatures (α~0.618=1/Φ) in CMB low-multipole data (Planck 2018 SMICA, nside=2048).

## Overview
The POC posits the universe as a self-organizing informational system, imprinting a pre-inflationary "harmonic echo" in CMB anomalies (e.g., low-ℓ suppression ~20% p≈0.01, PR4 arXiv:2506.22795v1). Analysis of 50 equatorial patches (10° radius) yields >3σ evidence: α_pico=0.988±0.015, ρ(f(α),1/Φ)=0.54±0.021, p=0.002 (KS asymmetry), p=0.007 (correlation), E-modes p=0.041. Phases:
- **Phase I:** Pipeline validation (SciPy cwt ricker, moments q=[-5,5]).
- **Phase II:** Cold Spot analysis (01/10/2025, hp.query_disc raio 5° l=209° b=-57°).
- **Phase III:** CAMB V(φ)=½m²φ² + (λ/Φ)φ⁴ (λ~10^{-120}, swampland Obied 2018 arXiv:1806.08362).

## Contents
- **Manifesto:** [manifest.pdf](manifest.pdf) — POC fundamentals, roadmap.
- **Research Log:** [RESEARCH_LOG.md](RESEARCH_LOG.md) — Tests 1-4, Python code.
- **Paper:** [POC-CMB-Harmonic-Paper.pdf](paper/POC-CMB-Harmonic-Paper.pdf).
- **Code:** [code/teste4_full_analysis.py](code/teste4_full_analysis.py), [Multifractal_Pipeline_Notebook.ipynb](code/Multifractal_Pipeline_Notebook.ipynb).
- **Preregistration:** [OSF DOI:10.17605/OSF.IO/VC59T](https://osf.io/vc59t).

## Replicability
1. Clone: `git clone https://github.com/leandroantoniodasilvapacheco/poc-cmb-harmonic.git`
2. Install: `pip install -r requirements.txt` (numpy, scipy, pywavelets, healpy, astropy, matplotlib).
3. Data: Download Planck SMICA FITS from [pla.esac.esa.int](https://pla.esac.esa.int/pla/aio/product-action?MAP.MAP_ID=COM_CMB_IQU-smica_2048_R3.00_full.fits).
4. Run: `python code/teste4_full_analysis.py --seed 42` (output: α_pico=0.988±0.015, ρ=0.54±0.021, p=0.002/0.007).

## Contribute
Join our open science journey! Fork, open issues, or submit pull requests. Contact: leandroantoniodasilvapacheco@proton.me.

## License
[CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/) — Public domain, unrestricted reuse.

## Acknowledgments
Planck Collaboration, HEALPix, NumPy, SciPy, CAMB. Inspired by Schwarz 2016 (arXiv:1510.07929), Higuchi 1988, Falconer 2003.

---
Leandro Antonio da Silva Pacheco, Brazil
