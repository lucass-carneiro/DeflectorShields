"""Deflector Shield Plotter.

Usage:
  plots.py single <parameter-file> <data-file>
  plots.py multiple <parameter-file> <data-file-folder>
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
from PIL import Image
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

    with ccf.ProcessPoolExecutor(max_workers=4) as executor:
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
        "libx264",
        "-crf",
        "0",
        f"{name}_anim.mp4"
    ])

    logger.info(f"Removing temporary frames folder")
    shutil.rmtree(anim_folder)

    logger.info(f"Done")


def plot_multiple_kernel(prefix, ipc_file_name, anim_folder, v, radius_x, radius_y, sigma_x, sigma_y, ent_img):
    ipc_file = os.path.join(prefix, ipc_file_name)
    df = pl.read_ipc(ipc_file, memory_map=False)

    it = int(df.item(0, 0))
    anim_file_name = os.path.join(anim_folder, f"{prefix}_{it:08d}.png")

    t = df.item(2, 0)

    plt.close("all")

    fig, ax = plt.subplots(1, 1, layout="tight")

    # Ship
    ent_x = v * t
    ent_y = 0.0

    ax.scatter(ent_x, ent_y, marker="o", color="black", s=2)

    ent_width, ent_height = ent_img.size
    ent_scale = 0.002

    ent_new_width = ent_width * ent_scale
    ent_new_height = ent_height * ent_scale

    ent_left = ent_x - ent_new_width / 2
    ent_right = ent_x + ent_new_width / 2
    ent_bottom = ent_y - ent_new_height / 2
    ent_top = ent_y + ent_new_height / 2

    ax.imshow(ent_img, extent=(ent_left, ent_right, ent_bottom, ent_top))

    # Particles
    for i in range(0, len(df.columns)):
        x = df.item(3, i)
        y = df.item(4, i)
        ax.scatter(x, y, marker="o", color="tab:red", s=2)

    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")

    title_a = r"$t = "
    title_b = r"$"
    ax.set_title(f"{title_a}{t:.2f}{title_b}")

    # Radii
    warp_bubble_start = plt.Circle(
        (v * t, 0.0),
        radius_x,
        fill=False,
        color="black"
    )

    warp_bubble_end = plt.Circle(
        (v * t, 0.0),
        sigma_x,
        fill=False,
        color="black",
        linestyle="--"
    )

    shield_start = plt.Circle(
        (v * t, 0.0),
        radius_y,
        fill=False,
        color="tab:blue"
    )

    shield_end = plt.Circle(
        (v * t, 0.0),
        sigma_y,
        fill=False,
        color="tab:blue",
        linestyle="--"
    )

    ax.add_patch(warp_bubble_start)
    ax.add_patch(warp_bubble_end)

    ax.add_patch(shield_start)
    ax.add_patch(shield_end)

    # Ranges
    ax.set_ylim(-2.0 * sigma_y, 2.0 * sigma_y)
    ax.set_xlim(ent_x - 2.0 * sigma_y, ent_x + 2.0 * sigma_y)

    # Save figure
    fig.savefig(anim_file_name, dpi=300)


def plot_multiple(prefix, parameters):
    logger.info("Plotting multiple particle data")

    if "AlcubierreSharp" in parameters["warp_drive_solution"]:
        v = parameters["warp_drive_solution"]["AlcubierreSharp"]["v"]

        radius_x = parameters["warp_drive_solution"]["AlcubierreSharp"]["radii_x"]["radius"]
        sigma_x = parameters["warp_drive_solution"]["AlcubierreSharp"]["radii_x"]["sigma"]

        radius_y = parameters["warp_drive_solution"]["AlcubierreSharp"]["radii_y"]["radius"]
        sigma_y = parameters["warp_drive_solution"]["AlcubierreSharp"]["radii_y"]["sigma"]
    else:
        raise RuntimeError(
            f"Parameter retrival no implemented for this warp drive"
        )

    anim_folder = f"{prefix}_anim"

    if not os.path.exists(anim_folder):
        logger.info(f"Creating temporary frames folder")
        os.mkdir(anim_folder)

    ipc_file_list = os.listdir(prefix)

    # Enterprise image
    ent_img = Image.open("resources/uss_enterprise.jpg")

    # for every file
    logger.info(f"Submitting plotting jobs to pool")

    with ccf.ProcessPoolExecutor(max_workers=8) as executor:
        for ipc_file in ipc_file_list:
            executor.submit(
                plot_multiple_kernel,
                prefix,
                ipc_file,
                anim_folder,
                v,
                radius_x,
                radius_y,
                sigma_x,
                sigma_y,
                ent_img
            )

        logger.info(f"Waiting for plotting jobs to finish")

    logger.info(f"Animating frames with ffmpeg")
    subprocess.run([
        "ffmpeg",
        "-y",
        "-framerate",
        "15",
        "-pattern_type",
        "glob",
        "-i",
        os.path.join(anim_folder, r"*.png"),
        "-c:v",
        "libx264",
        "-crf",
        "0",
        f"{prefix}_anim.mp4"
    ])

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
            single_particle_par_file["warp_drive_solution"]["AlcubierreSharp"]["v"],
            single_particle_par_file["warp_drive_solution"]["AlcubierreSharp"]["sigma"],
            single_particle_par_file["warp_drive_solution"]["AlcubierreSharp"]["radius"]
        )
    elif args["multiple"]:
        parameter_file = args["<parameter-file>"]
        output_file_prefix = args["<data-file-folder>"]

        parameters = read_parameter_file(parameter_file)

        plot_multiple(
            output_file_prefix,
            parameters
        )


if __name__ == '__main__':
    args = docopt(__doc__, version="Deflector Shield Plotter 1.0")
    main(args)
