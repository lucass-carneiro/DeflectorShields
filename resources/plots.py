"""Deflector Shield Plotter.

Usage:
  plots.py single <parameter-file> <data-file>
  plots.py (-h | --help)
  plots.py --version

Options:
  -h --help     Show this screen.
  --version     Show version.
"""
import logging
from docopt import docopt

import os
import sys
import subprocess
import shutil
import json
import concurrent.futures as ccf
import polars as pl
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use("agg")  # This prevents a GUI to lunch on multiple threads

logger = logging.getLogger(__name__)


# Fonts and text size
mpl.rcParams["font.size"] = 20
mpl.rcParams["text.usetex"] = True

# Figure sizes
mpl.rcParams['lines.linewidth'] = 1.3


def read_ipc(file):
    logger.info(f"Reading ipc file {file}")

    basename = os.path.basename(file)
    name, _ = os.path.splitext(basename)

    df = pl.read_ipc(file, memory_map=False)
    header = df.columns

    return (name, header, df)


def read_parameter_file(file_path):
    with open(file_path, "r") as file:
        return json.load(file)


def alcubierre_r(v, t, x, y, z):
    return np.sqrt((x - v*t)**2 + y**2 + z**2)


def alcubierre_f(v, sigma, radius, t, x, y, z):
    r = alcubierre_r(v, t, x, y, z)
    return (np.tanh(sigma * (r + radius)) - np.tanh(sigma * (r - radius)))/(2.0 * np.tanh(sigma * radius))


def plot_single_kernel(anim_file_name, v, t, x, y, X, Y, Z):
    plt.close("all")

    # Bubble
    plt.contourf(X, Y, Z, 1000, cmap="viridis")

    # Particle
    plt.scatter(x, y, marker="o", color="tab:red", s=2)

    # Ship
    plt.scatter(v * t, 0.0, marker="o", color="black", s=2)

    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")

    title_a = r"$t = "
    title_b = r"$"
    plt.title(f"{title_a}{t:.2f}{title_b}")

    plt.colorbar()

    plt.tight_layout()
    plt.savefig(anim_file_name, dpi=300)


def plot_single(data, v, sigma, radius):
    logger.info("Plotting single particle data")

    name, header, df = data

    anim_folder = f"{name}_anim"
    if not os.path.exists(anim_folder):
        logger.info(f"Creating temporary frames folder")

        os.mkdir(anim_folder)

    xy = np.linspace(-10.0, 10.0, 100)
    X, Y = np.meshgrid(xy, xy)

    nt = len(df)

    logger.info(f"Submitting plotting jobs to pool")

    with ccf.ProcessPoolExecutor(max_workers=8) as executor:
        for it in range(0, nt):
            anim_file_name = os.path.join(anim_folder, f"{name}_{it:08d}.png")
            t = df.item(it, header[2])
            x = df.item(it, header[3])
            y = df.item(it, header[4])

            Z = alcubierre_f(v, sigma, radius, t, X, Y, 0.0)

            executor.submit(
                plot_single_kernel,
                anim_file_name,
                v,
                t,
                x,
                y,
                X,
                Y,
                Z
            )

        logger.info(f"Waiting for plotting jobs to finish")

    logger.info(f"Animating frames with ffmpeg")
    subprocess.run([
        "ffmpeg",
        "-y",
        "-framerate",
        "15",
        "-i",
        os.path.join(anim_folder, f"{name}_%08d.png"),
        "-c:v",
        "libx264rgb",
        "-crf",
        "0",
        f"{name}_anim.mp4"
    ])

    logger.info(f"Removing temporary frames folder")
    shutil.rmtree(anim_folder)

    logger.info(f"Done")


def main(args):
    logging.basicConfig(stream=sys.stdout, level=logging.INFO)

    if args["single"]:
        parameter_file = args["<parameter-file>"]
        output_file = args["<data-file>"]

        single_particle_par_file = read_parameter_file(parameter_file)
        single_particle_data = read_ipc(output_file)

        plot_single(
            single_particle_data,
            single_particle_par_file["alcubierre_data"]["v"],
            single_particle_par_file["alcubierre_data"]["sigma"],
            single_particle_par_file["alcubierre_data"]["radius"]
        )


if __name__ == '__main__':
    args = docopt(__doc__, version="Deflector Shield Plotter 1.0")
    main(args)
