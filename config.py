"""
Module to load the configuration parameters read by `setrun.py`.
The default values are set in `defaults.yaml`.
The required and the custom parameters should be set in `config.yaml`.
"""
from pathlib import Path
from yaml import safe_load


def load_yaml(path):
    path = Path(path)
    if not path.is_file():
        return dict()
    with open(path) as file:
        contents = safe_load(file)
    return contents


defaults = load_yaml("defaults.yaml")
config = defaults | load_yaml("config.yaml")


if __name__ == "__main__":
    for k, v in defaults.items():
        print(f"{k} : {v}")
    print("---")
    for k, v in config.items():
        print(f"{k} : {v}")
