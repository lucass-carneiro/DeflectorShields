"""Deflector Shield Plotter.

Usage:
  plots.py [options] <parameter-file> <data-file-folder>
  plots.py (-h | --help)
  plots.py --version

Options:
  -h --help           Show this screen.
  --version           Show version.
  -b --follow-bubble  Follows the bubble.
  -s --save-pdf       Saves a PDF file along a regular PNG file
"""
import logging
from docopt import docopt

import os
import sys
import subprocess
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


def plot_multiple_kernel(prefix, ipc_file_name, anim_folder, bubble_speed, radius, sigma, ent_img, follow_bubble, shutdown_t, save_pdf):
    ipc_file = os.path.join(prefix, ipc_file_name)
    df = pl.read_ipc(ipc_file, memory_map=False)

    it = int(df.item(0, 0))
    anim_file_name = os.path.join(anim_folder, f"{prefix}_{it:08d}.png")

    if save_pdf:
        anim_file_name_pdf = os.path.join(
            anim_folder,
            f"{prefix}_{it:08d}.pdf"
        )

    # Current time
    t = df.item(2, 0)

    # Ship position
    ent_x = df.item(3, -1)
    ent_y = df.item(4, -1)

    # Bubble position
    bubble_x = bubble_speed * t
    bubble_y = 0.0

    plt.close("all")

    fig, ax = plt.subplots(1, 1, layout="tight")

    # Ship
    ax.scatter(ent_x, ent_y, marker="*", color="tab:blue", s=2)

    # Enterprise
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
    for i in range(0, len(df.columns) - 1):
        x = df.item(3, i)
        y = df.item(4, i)
        ax.scatter(x, y, marker="o", color="tab:red", s=2)

    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")

    title_a = r"$t = "
    title_b = r"$"
    ax.set_title(f"{title_a}{t:.2f}{title_b}")

    # Bubble
    warp_bubble_start = plt.Circle(
        (bubble_x, bubble_y),
        radius,
        fill=False,
        color="black"
    )

    warp_bubble_end = plt.Circle(
        (bubble_x, bubble_y),
        radius + sigma,
        fill=False,
        color="black",
        linestyle="--"
    )

    # Bubble center
    if shutdown_t is None or t <= shutdown_t:
        ax.scatter(bubble_x, bubble_y, marker="o", color="black", s=2)

        ax.add_patch(warp_bubble_start)
        ax.add_patch(warp_bubble_end)

    # Ranges
    range_factor = 3.0
    ax.set_ylim(-range_factor * sigma, range_factor * sigma)

    if follow_bubble:
        ax.set_xlim(
            bubble_x - range_factor * sigma,
            bubble_x + range_factor * sigma
        )
    else:
        ax.set_xlim(
            ent_x - range_factor * sigma,
            ent_x + range_factor * sigma
        )

    # Save figure
    fig.tight_layout()

    fig.savefig(anim_file_name, dpi=300, bbox_inches="tight")

    if save_pdf:
        fig.savefig(anim_file_name_pdf, dpi=300, bbox_inches="tight")


def plot_multiple(prefix, parameters, follow_bubble, save_pdf):
    logger.info("Plotting multiple particle data")

    if "Ours" in parameters["warp_drive_solution"]:
        u = parameters["warp_drive_solution"]["Ours"]["u"]

        radius = parameters["warp_drive_solution"]["Ours"]["radius"]
        sigma = parameters["warp_drive_solution"]["Ours"]["sigma"]

        shutdown_time = parameters["warp_drive_solution"]["Ours"]["ts"]
        shutdown_duration = parameters["warp_drive_solution"]["Ours"]["ds"]
        if shutdown_duration < 0.0:
            shutdown_t = None
        else:
            shutdown_t = shutdown_time + shutdown_duration
    else:
        raise RuntimeError(
            f"Parameter retrival not implemented for this warp drive"
        )

    anim_folder = f"{prefix}_anim"

    if not os.path.exists(anim_folder):
        logger.info(f"Creating temporary frames folder")
        os.mkdir(anim_folder)

    ipc_file_list = os.listdir(prefix)

    # Enterprise image
    ent_img = Image.open("resources/ship.png")

    # for every file
    logger.info(f"Submitting plotting jobs to pool")

    with ccf.ProcessPoolExecutor(max_workers=8) as executor:
        for ipc_file in ipc_file_list:
            executor.submit(
                plot_multiple_kernel,
                prefix,
                ipc_file,
                anim_folder,
                u,
                radius,
                sigma,
                ent_img,
                follow_bubble,
                shutdown_t,
                save_pdf
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

    parameter_file = args["<parameter-file>"]
    output_file_prefix = args["<data-file-folder>"]
    follow_bubble = bool(args["--follow-bubble"])
    save_pdf = bool(args["--save-pdf"])

    parameters = read_parameter_file(parameter_file)

    plot_multiple(
        output_file_prefix,
        parameters,
        follow_bubble,
        save_pdf
    )


if __name__ == '__main__':
    args = docopt(__doc__, version="Deflector Shield Plotter 1.0")
    main(args)
