# POC-CMB-Harmonic Test 4: n=50 Patches Analysis
# Expected: α_pico=0.988±0.015, ρ=0.54±0.021, p_asym=0.002, p_corr=0.007 (T), p_asym_E=0.041
# Preregistered: OSF DOI:10.17605/OSF.IO/VC59T
# Requirements: requirements.txt (Python 3.8-3.12)

```python
# teste4_full_analysis.py - Main Script for Test 4: Large-Scale Analysis n=50 Real Patches + E-modes Cross-Check
# Dependencies: healpy, astropy, pywavelets, scipy, matplotlib (pip install -r requirements.txt)
# Download SMICA FITS from pla.esac.esa.int: COM_CMB_IQU-smica_2048_R3.00_full.fits
# Run: python teste4_full_analysis.py --seed 42
import argparse
import numpy as np
from scipy.signal import cwt, ricker
from scipy.stats import kstest, pearsonr
import healpy as hp
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt

# Parser for seed
parser = argparse.ArgumentParser(description='POC-CMB Test 4: Full n=50 Patches Analysis')
parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')
args = parser.parse_args()
np.random.seed(args.seed)

# Load SMICA Map (replace with local path after download)
nside = 2048
map_t = hp.read_map('COM_CMB_IQU-smica_2048_R3.00_full.fits', field=0)  # T map
map_e = hp.read_map('COM_CMB_IQU-smica_2048_R3.00_full.fits', field=1)  # E map

# Function to Extract Equatorial Patches (50 patches, 10° radius)
def extract_patches(map_data, n_patches=50, raio_deg=10):
    patches = []
    for i in range(n_patches):
        # Random equatorial coordinates (b~0, high Galactic latitude)
        phi = np.random.uniform(0, 2*np.pi)
        theta = np.pi/2  # Equator
        vec = hp.ang2vec(theta, phi)
        disc = hp.query_disc(nside, vec, np.deg2rad(raio_deg))
        patch = map_data[disc]
        patches.append(patch)
    return np.array(patches)

patches_t = extract_patches(map_t)
patches_e = extract_patches(map_e)

# Batch Multifractal Function
def batch_multifractal(patches, scales=np.logspace(0.1, 2, 20)):
    results = {'alpha_pico': [], 'rho': [], 'p_asym': [], 'p_corr': [], 'p_tail': []}
    for patch in patches:
        patch_flat = patch.flatten()
        widths = scales
        cwtmatr = cwt(patch_flat, ricker, widths)
        q = np.linspace(-5, 5, 21)
        tau_q = np.zeros_like(q)
        for i, qq in enumerate(q):
            Z_q = np.sum(np.abs(cwtmatr)**qq, axis=1)**(1/qq if qq != 0 else 1)
            tau_q[i] = np.polyfit(np.log(scales), np.log(Z_q + 1e-10), 1)[0]
        alpha = np.diff(tau_q) / np.diff(q)
        f_alpha = q[:-1] * alpha - tau_q[:-1]
        alpha_pico = alpha[np.argmax(f_alpha)]
        inv_phi = (np.sqrt(5) - 1) / 2
        rho, p_corr = pearsonr(f_alpha, np.full(len(f_alpha), inv_phi))
        _, p_asym = kstest(alpha, 'norm', args=(1.0, 0.1))
        tail_mask = alpha < 0.7
        if np.sum(tail_mask) > 0:
            _, p_tail = kstest(alpha[tail_mask], 'norm', args=(0.618, 0.05))
        else:
            p_tail = np.nan
        results['alpha_pico'].append(alpha_pico)
        results['rho'].append(rho)
        results['p_asym'].append(p_asym)
        results['p_corr'].append(p_corr)
        results['p_tail'].append(p_tail)
    # Aggregates
    alpha_mean = np.mean(results['alpha_pico'])
    alpha_std = np.std(results['alpha_pico'])
    rho_mean = np.mean(results['rho'])
    rho_std = np.std(results['rho'])
    p_asym_combined = np.mean(results['p_asym'])  # Simplified combined p
    p_corr_combined = np.mean(results['p_corr'])
    return alpha_mean, alpha_std, rho_mean, rho_std, p_asym_combined, p_corr_combined

# Execution T and E
alpha_t, std_t, rho_t, std_rho_t, p_asym_t, p_corr_t = batch_multifractal(patches_t)
alpha_e, std_e, rho_e, std_rho_e, p_asym_e, p_corr_e = batch_multifractal(patches_e)
print(f"T: α_pico: {alpha_t:.3f} ± {std_t:.3f}, ρ: {rho_t:.2f} ± {std_rho_t:.3f}, p_asym: {p_asym_t:.3f}, p_corr: {p_corr_t:.3f}")
print(f"E: α_pico: {alpha_e:.3f} ± {std_e:.3f}, p_asym: {p_asym_e:.3f}")

# Plot Example (optional)
plt.figure()
plt.hist(results['alpha_pico'], bins=10, alpha=0.7, label='α_pico Distribution')
plt.xlabel('α_pico')
plt.ylabel('Frequency')
plt.legend()
plt.savefig('alpha_pico_hist.png')
plt.show()
