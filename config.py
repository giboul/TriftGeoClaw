"""
Module to load the configuration parameters read by `setrun.py`.
The default values are set in `defaults.yaml`.
The required and the custom parameters should be set in `config.yaml`.
"""
from pathlib import Path
from yaml import safe_load

with open("defaults.yaml") as file:
    config = safe_load(file)

if Path("config.yaml").exists():
    with open("config.yaml") as file:
        config.update(safe_load(file))

if __name__ == "__main__":
    print(config)
