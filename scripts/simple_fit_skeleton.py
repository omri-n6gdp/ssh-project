#!/usr/bin/env python3
"""
simple_fit_skeleton.py

Skeleton script for a simple 1D fit using iminuit.

"""

import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit

# ------------------------------
# 1. Made-up 1D data
# ------------------------------

# For simplicity, we'll create a 1D array of "measurements".
# In a real analysis, this might be something you have measured.
# Here we just generate some fixed numbers.

# Example: data drawn from a Gaussian with mean ~1.0 and sigma ~0.2 
rng = np.random.default_rng(42)
true_mean = 1.0
true_sigma = 0.2
n_data = 200
data = rng.normal(true_mean, true_sigma, size=n_data)

# 1D binning for the histogram
x_min, x_max = 0.0, 2.0
n_bins = 40
bin_edges = np.linspace(x_min, x_max, n_bins + 1)
bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# Build data histogram
data_hist, _ = np.histogram(data, bins=bin_edges)
print("[INFO] Data histogram counts:", data_hist)


# ------------------------------
# 2. Simple 1D model in parameter space
# ------------------------------

# We will use a simple Gaussian model with two parameters:
#   mu: mean
#   sigma: width
#
# The PDF (before normalisation) is:
#   f(x) = exp(-(x - mu)**2 / (2 * sigma**2))
#
# We will turn this into a binned expected count model for the histogram.


def unnormalised_gaussian(x, mu, sigma):
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2)


# ------------------------------
# 3. Binned Poisson NLL for the 1D histogram
# ------------------------------

def nll_binned_factory(data_hist, bin_edges):
    """
    Create a binned Poisson NLL function for the histogram 'data_hist'.

    Parameters:
      data_hist : 1D array of observed counts in each bin
      bin_edges : edges of the bins (length = n_bins+1)

    The model will be a Gaussian in x with parameters mu and sigma,
    normalised so that the total expected counts equals the total
    number of observed counts.
    """
    data_counts = data_hist.astype(float)
    n_events = data_counts.sum()
    x_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    def nll(mu, sigma):
        # Protect against non-positive sigma
        if sigma <= 0:
            return 1e12

        # Unnormalised model at bin centres
        f_vals = unnormalised_gaussian(x_centres, mu, sigma)
        f_vals_sum = np.sum(f_vals)
        if f_vals_sum <= 0 or not np.isfinite(f_vals_sum):
            return 1e12

        # Normalise so that sum of expected counts = n_events
        pdf_vals = f_vals / f_vals_sum
        E = n_events * pdf_vals  # expected counts per bin
        O = data_counts

        # Poisson NLL = sum(E - O log E)
        mask = E > 0
        E_masked = E[mask]
        O_masked = O[mask]
        return np.sum(E_masked - O_masked * np.log(E_masked + 1e-12))

    return nll

nll = nll_binned_factory(data_hist, bin_edges)

# ------------------------------
# 3b. Toy dataset generator
# ------------------------------

def generate_toy_histogram(mu, sigma, n_data, bin_edges, rng=None):
    """
    Generate one toy dataset drawn from a Gaussian(mu, sigma),
    histogram it using bin_edges, and return the counts.

    Parameters:
      mu       : Gaussian mean
      sigma    : Gaussian width
      n_data   : number of points to draw
      bin_edges: bin edges (length = n_bins+1)
      rng      : numpy Generator instance (optional; creates one if None)

    Returns:
      toy_hist : 1D array of integer counts, shape (n_bins,)
    """
    if rng is None:
        rng = np.random.default_rng()

    toy_data = rng.normal(mu, sigma, size=n_data)
    toy_hist, _ = np.histogram(toy_data, bins=bin_edges)
    return toy_hist

# ------------------------------
# 4. Set up and run a simple Minuit fit
# ------------------------------

initial_mu    = 0.5
initial_sigma = 0.3

m = Minuit(
    nll,
    mu=initial_mu,
    sigma=initial_sigma,
)
m.errordef = Minuit.LIKELIHOOD
m.limits['sigma'] = (0.17, None)   # <<< CHANGED: lower bound on sigma

print("\n[INFO] Starting MIGRAD fit on simple 1D model...")
m.migrad()
m.hesse()

print("Fit results (bounded sigma):")
print(f"  mu    = {m.values['mu']:.6f} ± {m.errors['mu']:.6f}")
print(f"  sigma = {m.values['sigma']:.6f} ± {m.errors['sigma']:.6f}")   # <<< CHANGED: watch this value
print(f"  sigma at lower bound: {m.values['sigma'] <= 0.17 + 1e-6}")     # <<< CHANGED: boundary flag
cov = m.covariance
cov_np = np.array(cov)
print("\nCovariance matrix (mu, sigma):\n", cov_np)

# ------------------------------
# 5. Generate one toy and inspect it
# ------------------------------

toy_rng = np.random.default_rng(0)   # separate RNG so toys don't interfere with data generation

toy_hist = generate_toy_histogram(
    mu=m.values['mu'],
    sigma=m.values['sigma'],
    n_data=n_data,
    bin_edges=bin_edges,
    rng=toy_rng,
)

print("\nToy histogram counts:")
print(toy_hist)
print(f"Total counts in toy: {toy_hist.sum()}  (expected {n_data})")

# ------------------------------
# 6. Toy experiment loop (bounded sigma)
# ------------------------------

n_toys = 500   # <<< CHANGED: more toys to resolve the asymmetric pull shape

toy_mu        = np.zeros(n_toys)
toy_sigma     = np.zeros(n_toys)
toy_err_mu    = np.zeros(n_toys)
toy_err_sigma = np.zeros(n_toys)

n_failed      = 0
n_at_boundary = 0   # <<< CHANGED: count toys where sigma hits the bound

for i in range(n_toys):

    t_hist = generate_toy_histogram(
        mu=m.values['mu'],
        sigma=m.values['sigma'],   # <<< CHANGED: truth is now the bounded fit result
        n_data=n_data,
        bin_edges=bin_edges,
        rng=toy_rng,
    )

    t_nll = nll_binned_factory(t_hist, bin_edges)

    tm = Minuit(
        t_nll,
        mu=m.values['mu'],
        sigma=m.values['sigma'],
    )
    tm.errordef  = Minuit.LIKELIHOOD
    tm.limits['sigma'] = (0.17, None)   # <<< CHANGED: same bound applied to every toy fit
    tm.migrad()
    tm.hesse()

    if not tm.valid or not tm.accurate:
        n_failed += 1
        toy_mu[i] = toy_sigma[i] = toy_err_mu[i] = toy_err_sigma[i] = np.nan
        continue

    # <<< CHANGED: flag toys whose sigma landed at (or within 1e-4 of) the boundary
    if tm.values['sigma'] <= 0.17 + 1e-4:
        n_at_boundary += 1

    toy_mu[i]        = tm.values['mu']
    toy_sigma[i]     = tm.values['sigma']
    toy_err_mu[i]    = tm.errors['mu']
    toy_err_sigma[i] = tm.errors['sigma']

print(f"\n[INFO] Toy loop done.  {n_toys - n_failed}/{n_toys} converged.")
print(f"[INFO] Toys with sigma at lower bound: {n_at_boundary}/{n_toys - n_failed}")  # <<< CHANGED

# ------------------------------
# 7. Pull distributions with Gaussian fits
# ------------------------------
from scipy.stats import norm                     # <<< CHANGED: new import (move to top of file ideally)
from scipy.optimize import curve_fit             # <<< CHANGED: for errors on the fit

mu_fit    = m.values['mu']
sigma_fit = m.values['sigma']

pull_mu    = (toy_mu    - mu_fit)    / toy_err_mu
pull_sigma = (toy_sigma - sigma_fit) / toy_err_sigma

valid      = np.isfinite(pull_mu) & np.isfinite(pull_sigma)
pull_mu    = pull_mu[valid]
pull_sigma = pull_sigma[valid]

x_ref     = np.linspace(-4, 4, 300)
gauss_ref = (1 / np.sqrt(2 * np.pi)) * np.exp(-0.5 * x_ref**2)

# <<< CHANGED: helper to fit a Gaussian and return errors via curve_fit
def fit_gaussian_to_pulls(pulls):
    """
    Fit a Gaussian to a pull array.

    Uses scipy.stats.norm.fit (MLE) for the central values, then
    curve_fit on the binned histogram to get approximate parameter errors.

    Returns
    -------
    loc, scale         : MLE mean and width
    err_loc, err_scale : 1-sigma errors from the covariance of curve_fit
    x_dense, y_dense   : smooth curve for plotting
    """
    # MLE fit — unbinned, exact
    loc, scale = norm.fit(pulls)

    # Binned curve_fit purely to extract a covariance matrix for error estimates.
    # We bin finely so the histogram approximates the PDF well.
    counts, edges = np.histogram(pulls, bins=20, density=True)
    centres = 0.5 * (edges[:-1] + edges[1:])

    def gaussian(x, mu, sigma):
        return norm.pdf(x, mu, sigma)

    try:
        popt, pcov = curve_fit(gaussian, centres, counts, p0=[loc, scale])
        err_loc, err_scale = np.sqrt(np.diag(pcov))
    except RuntimeError:
        # curve_fit failed to converge — fall back to NaN errors
        err_loc = err_scale = np.nan

    x_dense = np.linspace(-4, 4, 300)
    y_dense = norm.pdf(x_dense, loc, scale)

    return loc, scale, err_loc, err_scale, x_dense, y_dense


fig, axes = plt.subplots(1, 2, figsize=(11, 4))
fig.suptitle("Pull distributions — bounded sigma", fontsize=13)

for ax, pulls, label, is_sigma in zip(
    axes,
    [pull_mu, pull_sigma],
    [r"$\mu$", r"$\sigma$"],
    [False, True],
):
    n_valid = len(pulls)
    pmean   = np.mean(pulls)
    pstd    = np.std(pulls, ddof=1)

    # <<< CHANGED: Gaussian fit
    loc, scale, err_loc, err_scale, x_g, y_g = fit_gaussian_to_pulls(pulls)

    ax.hist(pulls, bins=20, density=True,
            color="steelblue", edgecolor="white", alpha=0.75,
            label="toy pulls")
    ax.plot(x_ref, gauss_ref, "k--", linewidth=1.2, label=r"$\mathcal{N}(0,1)$ ref")
    ax.plot(x_g,   y_g,       "tomato", linewidth=1.8, label="Gaussian fit")    # <<< CHANGED

    # <<< CHANGED: annotation now includes fit results
    textstr = (
        f"$n = {n_valid}$\n"
        f"sample mean $= {pmean:+.3f}$\n"
        f"sample std  $= {pstd:.3f}$\n"
        f"fit mean $= {loc:+.3f} \\pm {err_loc:.3f}$\n"
        f"fit width $= {scale:.3f} \\pm {err_scale:.3f}$"
    )
    if is_sigma:
        frac = n_at_boundary / max(n_valid, 1)
        textstr += f"\nat bound: {frac:.1%}"

    ax.text(0.97, 0.95, textstr, transform=ax.transAxes,
            fontsize=8.5, va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

    ax.set_xlabel(f"Pull of {label}", fontsize=11)
    ax.set_ylabel("Density", fontsize=11)
    ax.set_xlim(-4, 4)
    ax.legend(fontsize=9)
    ax.axvline(0, color="grey", linewidth=0.8, linestyle=":")

    if is_sigma:
        pull_at_bound = (0.17 - sigma_fit) / np.nanmean(toy_err_sigma[valid])
        ax.axvline(pull_at_bound, color="purple", linewidth=1.2,
                   linestyle="--", label=f"bound (pull≈{pull_at_bound:.1f})")
        ax.legend(fontsize=9)

    # <<< CHANGED: print fit results to terminal as well
    print(f"\nGaussian fit to pull_{label.strip('$')}:")
    print(f"  mean  = {loc:+.4f} ± {err_loc:.4f}")
    print(f"  width = {scale:.4f} ± {err_scale:.4f}")

fig.tight_layout()
plt.savefig("pull_distributions_bounded.png", dpi=150)
plt.show()