"""
Simple utilities to use geoclaw from a notebook.
"""  # TODO: check_claw()
from subprocess import call


def make(cmd, *args, **kwargs):
    return call(["make", cmd] + [f"{k}={v}" for k, v in kwargs.items()] + list(args))

def make_new(*args, **kwargs):
    return make("new", *args, **kwargs)

def make_data(*args, **kwargs):
    return make("data", *args, **kwargs)

def make_output(*args, **kwargs):
    return make("output", *args, **kwargs)
