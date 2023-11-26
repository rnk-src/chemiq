import numpy as np
from numpy import ndarray
from _operations.calculator import Calculator
from _operations.gas_laws import GasLaws


def boyle(p1=None, p2=None, v1=None, v2=None):
    inputs: ndarray[float] = np.array([p1, p2, v1, v2])
    names: ndarray[str] = np.array(['P1', 'P2', 'V1', 'V2'])
    arguments, missing = GasLaws.update_gas_law_arguments(inputs, names)
    solution: dict[str, float] = Calculator.solve_boyle_law(arguments=arguments)
    result = {
        'solution': solution,
        missing: solution
    }
    result.update(arguments)
    return result


def charles(v1=None, v2=None, t1=None, t2=None):
    inputs: ndarray[float] = np.array([v1, v2, t1, t2])
    names: ndarray[str] = np.array(['V1', 'V2', 'T1', 'T2'])
    arguments, missing = GasLaws.update_gas_law_arguments(inputs, names)
    solution: dict[str, float] = Calculator.solve_charles_law(arguments=arguments)
    result = {
        'solution': solution,
        missing: solution
    }
    result.update(arguments)
    return result


def gay_lussac(p1=None, p2=None, t1=None, t2=None):
    inputs: ndarray[float] = np.array([p1, p2, t1, t2])
    names: ndarray[str] = np.array(['P1', 'P2', 'T1', 'T2'])
    arguments, missing = GasLaws.update_gas_law_arguments(inputs, names)
    solution: dict[str, float] = Calculator.solve_gay_lussac_law(arguments=arguments)
    result = {
        'solution': solution,
        missing: solution
    }
    result.update(arguments)
    return result


def combined(p1=None, p2=None, v1=None, v2=None, t1=None, t2=None):
    inputs: ndarray[float] = np.array([p1, p2, v1, v2, t1, t2])
    names: ndarray[str] = np.array(['P1', 'P2', 'V1', 'V2', 'T1', 'T2'])
    arguments, missing = GasLaws.update_gas_law_arguments(inputs, names)
    solution: dict[str, float] = Calculator.solve_combined_gas_law(arguments=arguments)
    result = {
        'solution': solution,
        missing: solution
    }
    result.update(arguments)
    return result


def ideal(p=None, v=None, n=None, r=None, t=None):
    inputs: ndarray[float] = np.array([p, v, n, r, t])
    names: ndarray[str] = np.array(['P', 'V', 'N', 'R', 'T'])
    arguments, missing = GasLaws.update_gas_law_arguments(inputs, names)
    solution: dict[str, float] = Calculator.solve_ideal_gas_law(arguments=arguments)
    result = {
        'solution': solution,
        missing: solution
    }
    result.update(arguments)
