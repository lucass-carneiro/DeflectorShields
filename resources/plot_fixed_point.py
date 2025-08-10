import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import os


# Fonts and text size
mpl.rcParams["font.family"] = "cm"
mpl.rcParams["font.size"] = 20
mpl.rcParams["text.usetex"] = True
mpl.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"

# Figure sizes
mpl.rcParams["lines.linewidth"] = 2.0


def main():
    index = 1
    lims = 20
    point_size = 1
    legend_size = 7
    R = 4.0
    sigma = 4.0

    path_1 = os.path.join(f"fixed_point_{index}", "fixed_point_1.txt")
    path_2 = os.path.join(f"fixed_point_{index}", "fixed_point_2.txt")

    data_1 = np.loadtxt(path_1)
    data_2 = np.loadtxt(path_2)

    plt.close("all")

    plt.scatter(
        data_1[:, 0],
        data_1[:, 1],
        color="tab:blue",
        s=point_size
    )

    plt.scatter(
        data_2[:, 0],
        data_2[:, 1],
        color="tab:red",
        s=point_size
    )

    radius = plt.Circle(
        (0.0, 0.0),
        R + sigma,
        color="black",
        linestyle="--",
        fill=False,
        linewidth=2
    )

    radius_p_sigma = plt.Circle(
        (0.0, 0.0),
        R + 2.0 * sigma,
        color="black",
        linestyle="--",
        fill=False,
        linewidth=2
    )

    plt.gca().add_patch(radius)
    plt.gca().add_patch(radius_p_sigma)

    custom_lines = [
        Line2D([0], [0], color="tab:blue",),
        Line2D([0], [0], color="tab:red"),
    ]
    plt.legend(
        custom_lines,
        [
            r"$\mathrm{d}Y/\mathrm{d}t = 0$",
            r"$\mathrm{d}V^x/\mathrm{d}t = \mathrm{d}V^y/\mathrm{d}t = 0$",
        ],
        prop={"size": legend_size},
        loc="upper right"
    )

    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")

    plt.xlim(-lims, lims)
    plt.ylim(-lims, lims)

    plt.gca().set_aspect("equal")

    plt.tight_layout()
    plt.savefig(f"img/fixed_point_{index}.png", bbox_inches="tight", dpi=300)


main()
