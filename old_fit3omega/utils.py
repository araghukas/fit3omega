from typing import Sequence
from scipy.optimize import Bounds


def positive_bounds(guesses: Sequence[float],
                    min_frac: float = 0.1,
                    max_frac: float = 1.9) -> Bounds:
    """return Bounds of 10 to 190% the value for a sequence of guesses"""
    lb, ub = [], []
    for g in guesses:
        lb.append(min_frac * g)
        ub.append(max_frac * g)
    return Bounds(lb, ub, keep_feasible=True)
