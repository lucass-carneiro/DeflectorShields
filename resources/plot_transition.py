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

EPS = 1.0e-12


def f(x, y0, x0, dx):
    if x < x0:
        return y0
    elif x > (dx + x0):
        return 0.0
    else:
        return ((dx**3 + 4*dx**2*(x - x0) + 10*dx*(x - x0)**2 + 20*(x - x0)**3)*(dx - x + x0)**4*y0)/dx**7


def dfdx(x, y0, x0, dx):
    if x < x0 or (dx + x0) < x:
        return 0.0
    else:
        return (-140*(x - x0)**3*(dx - x + x0)**3*y0)/dx**7


def plot_f():
    y0 = 1.0
    x0 = 0.25
    dx = 0.5

    x = np.arange(0.0, 1.0, 1.0 / 100.0)

    plt.close("all")

    plt.plot(
        x,
        np.vectorize(f)(x, y0, x0, dx),
        color="black",
        label=r"$f(x;\, y_0, x_0, \Delta_x)$"
    )

    plt.plot(
        x,
        np.vectorize(dfdx)(x, y0, x0, dx),
        color="tab:red",
        label=r"$f^\prime(x;\, y_0, x_0, \Delta_x)$"
    )

    plt.axvline(x=x0, color="tab:blue",  linestyle='--')
    plt.axvline(x=x0 + dx, color="tab:blue",  linestyle='--')

    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")

    plt.legend(loc="upper right", prop={"size": 14})

    plt.tight_layout()
    plt.savefig("img/sharp_transition.pdf")


def main():
    plot_f()


main()
