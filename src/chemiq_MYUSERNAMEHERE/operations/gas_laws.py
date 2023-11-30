"""
Offers methods to aid with Gas Law calculations.
"""

from typing import Dict

from numpy import ndarray
import numpy as np

from src.chemiq_MYUSERNAMEHERE.base.gas_laws import GasLaws


def boyle(p1=None, p2=None, v1=None, v2=None) -> float:  # pylint: disable=invalid-name
    """
    Solves the Boyle gas law given all but one variable.
    :param p1: First pressure value.
    :param p2: Second pressure value.
    :param v1: First volume value.
    :param v2: Second volume value.
    :return: The missing value that was not given as an argument.
    """
    inputs: ndarray[float] = np.array([p1, p2, v1, v2])
    names: ndarray[str] = np.array(['P1', 'P2', 'V1', 'V2'])
    arguments: Dict[str, float] = GasLaws.update_gas_law_arguments(inputs, names)
    return GasLaws.solve_boyle(arguments=arguments)


def charles(v1=None, v2=None, t1=None, t2=None) -> float:  # pylint: disable=invalid-name
    """
    Solves the Charles gas law given all but one variable.
    :param v1: First volume value.
    :param v2: Second volume value.
    :param t1: First temperature value (Kelvin).
    :param t2: Second temperature value (Kelvin).
    :return: The missing value that was not given as an argument.
    """
    inputs: ndarray[float] = np.array([v1, v2, t1, t2])
    names: ndarray[str] = np.array(['V1', 'V2', 'T1', 'T2'])
    arguments: Dict[str, float] = GasLaws.update_gas_law_arguments(inputs, names)
    return GasLaws.solve_charles(arguments=arguments)


def gay_lussac(p1=None, p2=None, t1=None, t2=None) -> float:  # pylint: disable=invalid-name
    """
    Solves the Gay-Lussac gas law given all but one variable.
    :param p1: First pressure value.
    :param p2: Second pressure value.
    :param t1: First temperature value (Kelvin).
    :param t2: Second temperature value (Kelvin).
    :return: The missing value that was not given as an argument.
    """
    inputs: ndarray[float] = np.array([p1, p2, t1, t2])
    names: ndarray[str] = np.array(['P1', 'P2', 'T1', 'T2'])
    arguments: Dict[str, float] = GasLaws.update_gas_law_arguments(inputs, names)
    return GasLaws.solve_gay_lussac(arguments=arguments)


def combined(p1=None, p2=None, v1=None, v2=None, t1=None, t2=None) -> float:
    # pylint: disable=invalid-name
    # pylint: disable=too-many-arguments
    """
    Solves the combined gas law given all but one variable.
    :param p1: First pressure value.
    :param p2: Second pressure value.
    :param v1: First volume value.
    :param v2: Second volume value.
    :param t1: First temperature value (Kelvin).
    :param t2: Second temperature value (Kelvin).
    :return: The missing value that was not given as an argument.
    """
    inputs: ndarray[float] = np.array([p1, p2, v1, v2, t1, t2])
    names: ndarray[str] = np.array(['P1', 'P2', 'V1', 'V2', 'T1', 'T2'])
    arguments: Dict[str, float] = GasLaws.update_gas_law_arguments(inputs, names)
    return GasLaws.solve_combined_gas_law(arguments=arguments)


def ideal(p=None, v=None, n=None, r=None, t=None) -> float:  # pylint: disable=invalid-name
    """
    Solves the ideal gas law given all but one variable.
    :param p: The value of the pressure.
    :param v: The value of the volume.
    :param n: The number of moles.
    :param r: The value of the ideal gas constant.
    :pram t: The value of the temperature (Kelvin).
    """
    inputs: ndarray[float] = np.array([p, v, n, r, t])
    names: ndarray[str] = np.array(['P', 'V', 'N', 'R', 'T'])
    arguments: Dict[str, float] = GasLaws.update_gas_law_arguments(inputs, names)
    return GasLaws.solve_ideal_gas_law(arguments=arguments)
