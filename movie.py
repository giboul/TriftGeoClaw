#!/usr/bin/env python
from PIL import Image
from pathlib import Path
import argparse


def make_gif_movie(pattern: str, out_file, **kwargs) -> None:
    """
    From the path pattern, assemble images into a gif movie.
    
    Parameters
    ----------
    pattern: str
        The pattern that file names have to match
    file: str
        Path of the output movie.
    **kwargs
        Parameters passed to PIL.Image.save`
    """
    paths = sorted(Path().glob(pattern))
    images = [Image.open(p) for p in paths]
    images[0].save(out_file, append_images=images[1:], **kwargs)


parser = argparse.ArgumentParser()
parser.add_argument('--pattern', default="_plots/*.png")
parser.add_argument('--out_file', default="movie.gif")
parser.add_argument('--duration', default=200, type=int)
parser.add_argument('--loop', default=0, type=int)
parser.add_argument('--save_all', default=True, type=bool)
parser.add_argument('--optimize', default=False, type=bool)
args = parser.parse_args()


if __name__ == "__main__":
    make_gif_movie(**args.__dict__)
