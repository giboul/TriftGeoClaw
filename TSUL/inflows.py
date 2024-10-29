from argparse import ArgumentParser
from pathlib import Path
from yaml import safe_load

import src_inflows
import bc_inflows

projdir = Path().absolute().parent
with open(projdir / "config.yaml") as file:
    config = safe_load(file)["TSUL"]

parser = ArgumentParser()
parser.add_argument("avid", nargs="?", default="")
parser.add_argument("-p", "--plot", action="store_true")
parser.add_argument("-m", "--movie", action="store_true")
args = parser.parse_args()

if config["inflow"] == "bc":
    module = bc_inflows
elif config["inflow"] == "src":
    module = src_inflows
else:
    raise ValueError(f"Inflow mode '{config['inflow']}' is neither 'bc' nor 'src'.")
module.args = args

if args.plot or args.movie:
    module.plot(args.movie)
else:
    module.write()
