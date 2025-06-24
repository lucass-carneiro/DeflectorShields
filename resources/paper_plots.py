import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


# Fonts and text size
mpl.rcParams["font.family"] = "cm"
mpl.rcParams["font.size"] = 20
mpl.rcParams["text.usetex"] = True
mpl.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"

# Figure sizes
mpl.rcParams["lines.linewidth"] = 2.0


def f(sigma, R, x):
    if R < x and x < (R + sigma):
        return 1.0 / (1.0 + np.exp(sigma * (1.0 / (R - x)  + 1.0 / (R - x + sigma))))
    elif R < x:
        return 0.0
    else:
        return 1.0


def plot_f():
    R = 0.25
    sigma = 0.5
    
    x = np.arange(0.0, 1.0, 1.0 / 100.0)
    
    plt.close("all")
    plt.plot(x, np.vectorize(f)(sigma, R, x), color="black")

    plt.xlabel(r"$r$")
    plt.ylabel(r"$f(r;\, 0.5, 0.25)$")

    plt.tight_layout()
    plt.savefig("img/sharp_transition.pdf")


def main():
    plot_f()


main()
