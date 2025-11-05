from pathlib import Path
from yaml import safe_load

with open("defaults.yaml") as file:
    config = safe_load(file)

if Path("config.yaml").exists():

    with open("config.yaml") as file:
        custom = safe_load(file)

    config.update(custom)

if __name__ == "__main__":
    print(config)
