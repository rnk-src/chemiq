import numpy as np

from sympy.core import Symbol
from sympy.solvers import solve
from numpy import ndarray


class GasLaws:
    @classmethod
    def update_gas_law_arguments(cls, inputs: ndarray[float], names: ndarray[str]) -> dict[str, float] and str:
        arguments = dict()
        missing = None
        for input_value, name in zip(inputs, names):
            if input_value is not None:
                arguments.update({str(name): float(input_value)})
            else:
                missing = str(name)
        return arguments, missing

    @classmethod
    def solve_boyle(cls, **kwargs) -> float:
        missing = np.array([Symbol(letter) for letter in ['P1', 'P2', 'V1', 'V2'] if letter not in kwargs])
        if len(missing) != 1:
            raise ValueError("Expected 3 out of 4 arguments.")
        missing = missing[0]
        p1, p2, v1, v2 = [kwargs.get(letter, missing) for letter in ['P1', 'P2', 'V1', 'V2']]
        return float(solve((p1 * v1) - (p2 * v2), missing)[0])

    @classmethod
    def solve_charles(cls, **kwargs) -> float:
        missing = np.array([Symbol(letter) for letter in ['V1', 'V2', 'T1', 'T2'] if letter not in kwargs])
        if len(missing) != 1:
            raise ValueError("Expected 3 out of 4 arguments.")
        missing = missing[0]
        v1, v2, t1, t2 = [kwargs.get(letter, missing) for letter in ['V1', 'V2', 'T1', 'T2']]
        return float(solve((v1/t1)-(v2/t2), missing)[0])

    @classmethod
    def solve_gay_lussac(cls, **kwargs) -> float:
        missing = np.array([Symbol(letter) for letter in ['P1', 'P2', 'T1', 'T2'] if letter not in kwargs])
        if len(missing) != 1:
            raise ValueError(f"Expected 3 out of 4 arguments.")
        missing = missing[0]
        p1, p2, t1, t2 = [kwargs.get(letter, missing) for letter in ['P1', 'P2', 'T1', 'T2']]
        return float(solve((p1 / t1) - (p2 / t2), missing)[0])

    @classmethod
    def solve_combined_gas_law(cls, **kwargs) -> float:
        missing = np.array([Symbol(letter) for letter in ['P1', 'P2', 'T1', 'T2', 'V1', 'V2'] if letter not in kwargs])
        if len(missing) != 1:
            raise ValueError("Expected 5 out of 6 arguments.")
        missing = missing[0]
        p1, p2, t1, t2, v1, v2 = [kwargs.get(letter, missing) for letter in ['P1', 'P2', 'T1', 'T2', 'V1', 'V2']]
        return float(solve(((p1 * v1)/t1 - (p2 * v2)/t2), missing)[0])

    @classmethod
    def solve_ideal_gas_law(cls, **kwargs) -> float:
        missing = np.array([Symbol(letter) for letter in ['P', 'V', 'N', 'R', 'T'] if letter not in kwargs])
        if len(missing) != 1:
            raise ValueError("Expected 4 out of 5 arguments.")
        missing = missing[0]
        p, v, n, r, t = [kwargs.get(letter, missing) for letter in ['P', 'V', 'N', 'R', 'T']]
        return float(solve(((p * v) - (n * r * t)), missing)[0])
