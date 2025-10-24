from yaml import safe_load

with open("config.yaml") as file:
    config = safe_load(file)
