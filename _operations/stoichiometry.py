from typing import Any

import numpy as np

from molecule import Molecule
from numpy import ndarray
from balancer import Balancer


class NotSameSizeException(Exception):
    """Raised when the number of reactants and moles are not the same."""


class Stoichiometry:

    @classmethod
    def moles_to_grams(cls, moles: float, formula: str):
        molecule = Molecule.parse(formula)
        return moles * molecule.get_molar_mass()

    @classmethod
    def grams_to_moles(cls, grams: float, formula: str):
        molecule = Molecule.parse(formula)
        return grams / molecule.get_molar_mass()

    @classmethod
    def convert_grams_to_moles(cls, reactants: ndarray[Molecule], grams: ndarray[float]) -> ndarray[float]:
        moles: ndarray[float] = np.empty(reactants.size)
        for index, gram_amount in enumerate(grams):
            moles[index] = Stoichiometry.grams_to_moles(gram_amount, reactants[index])
        return moles

    @classmethod
    def limiting_reactant_moles(cls, reactants: ndarray[Molecule], reactant_coefficients: ndarray[float],
                                moles: ndarray[float]) -> Molecule:
        if reactants.size != moles.size:
            raise NotSameSizeException(f"Size of reactants array ({reactants.size}) does not match moles array "
                                       f"({moles.size}).")
        elif reactants.size != reactant_coefficients.size:
            raise NotSameSizeException(f"Size of reactants array ({reactants.size}) does not match coefficients array "
                                       f"({reactant_coefficients.size}).")

        division = moles / reactant_coefficients
        index_of_minimum = np.argmin(division)

        return reactants[index_of_minimum].item()

    @classmethod
    def limiting_reactant_moles_without_coefficients(cls, reactants: ndarray[Molecule], products: ndarray[Molecule],
                                                     moles: ndarray[float]) -> Molecule:
        reactant_coefficients = Balancer.balance_equation(reactants, products)[:reactants.size]
        return Stoichiometry.limiting_reactant_moles(reactants, reactant_coefficients, moles)

    @classmethod
    def limiting_reactant_grams(cls, reactants: ndarray[Molecule], reactant_coefficients: ndarray[float],
                                grams: ndarray[float]) -> Molecule:
        moles = Stoichiometry.convert_grams_to_moles(reactants, grams)
        return Stoichiometry.limiting_reactant_moles(reactants, reactant_coefficients, moles)

    @classmethod
    def limiting_reactant_grams_without_coefficients(cls, reactants: ndarray[Molecule], products: ndarray[Molecule],
                                                     grams: ndarray[float]):
        moles = Stoichiometry.convert_grams_to_moles(reactants, grams)
        return Stoichiometry.limiting_reactant_moles_without_coefficients(reactants, products, moles)
