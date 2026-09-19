"""Finite correlation extension, declared before execution."""
from fluctuations import run as run_fluctuations
def run():
    return run_fluctuations(correlation=.2)
