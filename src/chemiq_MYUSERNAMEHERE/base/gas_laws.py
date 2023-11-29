"""
This module contains the Gas Laws operations.
"""

import numpy as np

from sympy.core import Symbol
from sympy.solvers import solve
from numpy import ndarray


class GasLaws:
    """
    This class is responsible for handling basic Gas Laws operations.
    """
    @classmethod
    def update_gas_law_arguments(cls, inputs: ndarray[float],
                                 names: ndarray[str]) -> dict[str, float]:
        """
        Prepares arguments to be passed to the gas law functions.

        :param inputs: An array of floats.
        :param names: Names of variables corresponding to the inputs.

        :return: A dictionary containing variable names as the key and floats as the
         value of the variable.
        """
        arguments = {}
        for input_value, name in zip(inputs, names):
            if input_value is not None:
                arguments.update({str(name): float(input_value)})
        return arguments

    @classmethod
    def solve_boyle(cls, **kwargs) -> float:
        """
        Solves Boyle's gas law.

        :param kwargs: Three of the four variables involved in Boyle's gas law.
        :return: Float representing the value of the missing variable.
        """
        missing = np.array([Symbol(letter) for letter in ['P1', 'P2', 'V1', 'V2']
                            if letter not in kwargs])
        if len(missing) != 1:
            raise ValueError("Expected 3 out of 4 arguments.")
        missing = missing[0]
        first_pressure, second_pressure, first_volume, second_volume = \
            [kwargs.get(letter, missing) for letter in ['P1', 'P2', 'V1', 'V2']]
        return float(
            solve(
                (first_pressure * first_volume) - (second_pressure * second_volume),
                missing
            )[0]
        )

    @classmethod
    def solve_charles(cls, **kwargs) -> float:
        """
        Solves Charles' gas law.

        :param kwargs: Three of the four variables involved in Charles' gas law.
        :return: Float representing the value of the missing variable.
        """
        missing = np.array([Symbol(letter) for letter in ['V1', 'V2', 'T1', 'T2']
                            if letter not in kwargs])
        if len(missing) != 1:
            raise ValueError("Expected 3 out of 4 arguments.")
        missing = missing[0]
        first_volume, second_volume, first_temperature, second_temperature = \
            [kwargs.get(letter, missing) for letter in ['V1', 'V2', 'T1', 'T2']]
        return float(
            solve(
                (first_volume/first_temperature)-(second_volume/second_temperature),
                missing
            )[0]
        )

    @classmethod
    def solve_gay_lussac(cls, **kwargs) -> float:
        """
        Solves Gay-Lussac's gas law.

        :param kwargs: Three of the four variables involved in Gay-Lussac's gas law.
        :return: Float representing the value of the missing variable.
        """
        missing = np.array([Symbol(letter) for letter in ['P1', 'P2', 'T1', 'T2']
                            if letter not in kwargs])
        if len(missing) != 1:
            raise ValueError("Expected 3 out of 4 arguments.")
        missing = missing[0]
        first_pressure, second_pressure, first_temperature, second_temperature = \
            [kwargs.get(letter, missing) for letter in ['P1', 'P2', 'T1', 'T2']]
        return float(
            solve(
                (first_pressure / first_temperature) - (second_pressure / second_temperature),
                missing
            )[0]
        )

    @classmethod
    def solve_combined_gas_law(cls, **kwargs) -> float:
        """
        :param kwargs: Four of the five arguments required to solve the combined gas law.

        :return: A float representing the value of the missing parameter.
        """
        missing = np.array([Symbol(letter) for letter in ['P1', 'P2', 'T1', 'T2', 'V1', 'V2']
                            if letter not in kwargs])
        if len(missing) != 1:
            raise ValueError("Expected 5 out of 6 arguments.")
        missing = missing[0]
        (first_pressure, second_pressure, first_temperature, second_temperature,
         first_volume, second_volume) = [kwargs.get(letter, missing)
										 for letter in ['P1', 'P2', 'T1', 'T2', 'V1', 'V2']]
        return float(
            solve(
                ((first_pressure * first_volume)/first_temperature -
                 (second_pressure * second_volume)/second_temperature),
                missing
            )[0]
        )

    @classmethod
    def solve_ideal_gas_law(cls, **kwargs) -> float:
        """
        :param kwargs: Four of the five arguments required to solve the ideal gas law.

        :return: A float representing the value of the missing parameter.
        """
        missing = np.array([Symbol(letter) for letter in ['P', 'V', 'N', 'R', 'T']
                            if letter not in kwargs])
        if len(missing) != 1:
            raise ValueError("Expected 4 out of 5 arguments.")
        missing = missing[0]
        pressure, volume, moles, ideal_gas_constant, temperature = \
            [kwargs.get(letter, missing) for letter in ['P', 'V', 'N', 'R', 'T']]
        return float(
            solve(
                ((pressure * volume) - (moles * ideal_gas_constant * temperature)),
                missing
            )[0]
        )
