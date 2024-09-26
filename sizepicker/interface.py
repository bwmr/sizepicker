from pathlib import Path

import click

from sizepicker.picker import Picker


@click.command()
@click.option(
    "--radius",
    default=80,
    show_default=True,
    help="Particle radius [A].")
@click.option(
    "--angpix",
    default=None,
    type=float,
    show_default=True,
    help="Override pixelsize in header [A/pix].")
@click.option(
    "--contamination_binning",
    default=16,
    show_default=True,
    help="Binning used to filter out contaminations.")
@click.option(
    "--gaussian_cont",
    is_flag=True,
    default=False,
    show_default=True,
    help="Smooth background by subtracting gaussian-filtered tomogram.")
@click.option(
    "--sigma_cont",
    default=300,
    show_default=True,
    help="If smoothing background, apply gaussian filter with this sigma [A].")
@click.option(
    "--stdtimes_cont",
    default=2.5,
    show_default=True,
    help="Voxels this many SD below mean are considered contamination.")
@click.option(
    "--minsize_cont",
    default=50,
    show_default=True,
    help="Keep only contaminations of this size and above [binned px].")
@click.option(
    "--dilate_cont",
    default=200,
    show_default=True,
    help="Dilate contamination mask by this much [A].")
@click.option(
    "--radius_times",
    default=4,
    show_default=True,
    help="Picks can be this close together [particle radii].")
@click.option(
    "--inhibit",
    is_flag=True,
    default=False,
    show_default=True,
    help="Use more elaborate algorithm to get more picks.")
@click.option(
    "--detection_z",
    default=200,
    show_default=True,
    help="Only consider picks in a central slab of the volume of this extent [px].")
@click.option(
    "--stdtimes_pick",
    default=1.5,
    show_default=True,
    help="Only consider picks with a central SD this many SD above mean.")
@click.option(
    "--remove_edge",
    is_flag=True,
    default=False,
    show_default=True,
    help="Only keep particles with a higher foreground than background SD.")
@click.argument(
    "tomograms",
    nargs=-1,
    type=click.Path(exists=True))
@click.argument(
    "output_dir",
    nargs=1,
    type=click.Path(writable=True, dir_okay=True))
def cli(radius,
        angpix,
        contamination_binning,
        gaussian_cont,
        sigma_cont,
        stdtimes_cont,
        minsize_cont,
        dilate_cont,
        radius_times,
        inhibit,
        detection_z,
        stdtimes_pick,
        remove_edge,
        tomograms,
        output_dir):
    """Pick particles based on diameter.

    Takes tomograms as input. Writes all outputs to output_dir.

    """
    for tomo in tomograms:

        tomo = Path(tomo)
        output_dir = Path(output_dir)

        print(f"Working on {tomo.name}")

        p = Picker(tomo,
                   output_dir,
                   radius,
                   contamination_binning,
                   angpix)

        mask = p.getcont(gaussian_cont,
                         sigma_cont,
                         stdtimes_cont,
                         minsize_cont,
                         dilate_cont)

        print(f"Contamination mask written to {output_dir.stem}")

        boxes, particles = p.detect(mask,
                                    radius_times,
                                    inhibit,
                                    detection_z)

        print(f"Found {len(boxes)} initial picks.")

        metrics, threshold = p.prefilt(particles,
                                       stdtimes_pick)

        boxs_XYZ = p.filt(boxes,
               metrics,
               threshold,
               remove_edge)

        print(f"Wrote {len(boxs_XYZ)} picks to {output_dir.stem}")
        print("\n")
