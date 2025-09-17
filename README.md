# POC-CMB-Harmonic

[![DOI](https://img.shields.io/badge/OSF%20DOI-10.17605/OSF.IO/VC59T-blue)](https://osf.io/vc59t) [![License: CC0](https://img.shields.io/badge/License-CC0%201.0-brightgreen)](https://creativecommons.org/publicdomain/zero/1.0/) 

* **OSF Project DOI**: 10.17605/OSF.IO/VC59T (osf.io/vc59t) – Dynamic repository for collaborations, code, and data (SMICA maps from Planck 2018, ID: COM_CMB_IQU-smica_2048_R3.00_full.fits). Created: September 15, 2025; Last Updated: September 16, 2025; Contributor: LEANDRO ANTONIO DA SILVA PACHECO. This DOI ensures persistence for the editable project, always directing to the current content.

* **OSF Registration DOI**: 10.17605/OSF.IO/M7HYX – Immutable preregistration snapshot (OSF Preregistration type, registered on September 16, 2025), ensuring transparency for analyses and hypotheses (e.g., H1/H0 golden ratio tests). Associated with project osf.io/vc59t; Internet Archive Link: https://archive.org/details/osf-registrations-m7hyx-v1. Note: Accessible via OSF login; make public for global indexing on doi.org and Google Scholar.
  
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

## Update
Zenodo v2.1: DOI 10.5281/zenodo.17148203 (Is new version of v1 DOI 10.5281/zenodo.17126843; OSF Projeto: 10.17605/OSF.IO/VC59T; Registro: 10.17605/OSF.IO/M7HYX)

## Contribute to POC-CMB-Harmonic
Open issues or PRs using our templates for bugs, new features, or collaboration proposals (hypotheses, code, data). See Issues for templates.

## Contact: leandroantoniodasilvapacheco@proton.me

## License
[CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/) — Public domain, unrestricted reuse.

## Acknowledgments
Planck Collaboration, HEALPix, NumPy, SciPy, CAMB. Inspired by Schwarz 2016 (arXiv:1510.07929), Higuchi 1988, Falconer 2003.

---
Leandro Antonio da Silva Pacheco, Brazil
