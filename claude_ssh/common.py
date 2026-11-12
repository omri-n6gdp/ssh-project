"""Shared helpers: paths, plot styling, results bookkeeping."""
import os, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT    = os.path.dirname(os.path.abspath(__file__))
FIGDIR  = os.path.join(ROOT, "figures")
DATADIR = os.path.join(ROOT, "data")
NETDIR  = os.path.join(ROOT, "netlists")
for d in (FIGDIR, DATADIR, NETDIR):
    os.makedirs(d, exist_ok=True)

RESULTS = os.path.join(DATADIR, "results.json")

plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "font.size": 11,
    "axes.titlesize": 10.5,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "axes.axisbelow": True,
    "legend.framealpha": 0.92,
    "lines.linewidth": 1.8,
})

# colour-blind friendly palette (Okabe-Ito)
# the three Monte-Carlo disorder levels use maximally distinct hues
CLR = dict(tbr="#0072B2", theory="#D55E00", bulk="#999999",
           odd="#0072B2", even="#009E73", warn="#CC79A7",
           s1="#56B4E9",     # 1%  -- sky blue
           s5="#009E73",     # 5%  -- green
           s10="#D55E00")    # 10% -- vermillion


def newfig(w=7.4, h=4.3):
    """One figure sized to sit above its bullet block on a report page."""
    return plt.subplots(figsize=(w, h))


def save(fig, name):
    path = os.path.join(FIGDIR, name)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print(f"  figure -> figures/{name}")
    return path


def load_results():
    if os.path.exists(RESULTS):
        with open(RESULTS) as fh:
            return json.load(fh)
    return {}


def update_results(section, payload):
    r = load_results()
    r[section] = payload
    with open(RESULTS, "w") as fh:
        json.dump(r, fh, indent=2)
    print(f"  results['{section}'] written")
