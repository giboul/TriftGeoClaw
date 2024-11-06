from sys import argv
from subprocess import call
from pathlib import Path
from yaml import safe_load

projdir = Path(__file__).parents[1]
with open(projdir / "config.yaml") as file:
    config = safe_load(file)["TSUL"]

module = dict(bc="bc_inflows.py", src="src_inflows.py")[config["inflow"]]

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
    try:
        module.write()
    except:
        ...
