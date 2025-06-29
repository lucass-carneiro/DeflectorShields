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

def E(x, y0, x0, dx):
    x1 = 1.0 / (x0 - x)
    x2 = 1.0 / (x0 - x + dx)
    return dx * (x1 + x2)


def R(x, y0, x0, dx):
    x1 = x0 - x
    num = y0 * dx * (2.0 * x1**2 + 2.0 * x1 * dx + dx**2)
    den = 4.0 * x1**2 * (x1 + dx)**2
    return num / den


def f(x, y0, x0, dx):
    if x0 < x and x < (x0 + dx):
        return y0 / (1.0 + np.exp(E(x, y0, x0, dx)))
    elif x > x0:
        return 0.0
    else:
        return 1.0


def dfdx(x, y0, x0, dx):
    if x0 < x and x < (x0 + dx):
        sech = 1.0 / np.cosh(E(x, y0, x0, dx) / 2.0)
        return R(x, y0, x0, dx) * sech**2
    else:
        return 0.0


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
        linestyle="--",
        label=r"$f^\prime(x;\, y_0, x_0, \Delta_x)$"
    )
    
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")

    plt.legend(loc="upper right", prop={"size": 14})

    plt.tight_layout()
    plt.savefig("img/sharp_transition.pdf")


def main():
    plot_f()


main()
