"""
This module contains the Stoichiometry class, aiding in stoichiometry calculations.
"""
from numpy import ndarray
import numpy as np

from src.chemiq_MYUSERNAMEHERE.base.molecule import Molecule
from src.chemiq_MYUSERNAMEHERE.base.balancer import Balancer


class NotSameSizeException(Exception):
    """
    Raised when the number of reactants and moles are not the same.
    """


class Stoichiometry:
    """
    This class is responsible for performing stoichiometry calculations.
    """

    @classmethod
    def moles_to_grams(cls, moles: float, formula: str):
        """
        Converts moles of a formula to grams.

        :param moles: Number of moles of the molecule formula.
        :param formula: Formula of the molecule.
        """
        molecule = Molecule.parse(formula)
        return moles * molecule.get_molar_mass()

    @classmethod
    def grams_to_moles(cls, grams: float, formula: str):
        """
        Converts grams of a formula to moles.

        :param grams: Number of grams of the molecule formula.
        :param formula: Formula of the molecule.
        """
        molecule = Molecule.parse(formula)
        return grams / molecule.get_molar_mass()

    @classmethod
    def convert_grams_to_moles(cls, molecules: ndarray[Molecule],
                               grams: ndarray[float]) -> ndarray[float]:
        """
        Converts grams of molecules to moles.

        :param molecules: Molecules whose grams are converted to moles.
        :param grams: Number of grams of each molecule in the molecule array.
        """
        moles: ndarray[float] = np.empty(molecules.size)
        for index, gram_amount in enumerate(grams):
            moles[index] = Stoichiometry.grams_to_moles(gram_amount, molecules[index])
        return moles

    @classmethod
    def limiting_reactant_moles(cls, reactants: ndarray[Molecule],
                                reactant_coefficients: ndarray[float],
                                moles: ndarray[float]) -> Molecule:
        """
        Returns the limiting reactant in a chemical equation without coefficients with moles.

        :param reactants: Array of the reactants represented as Molecules.
        :param reactant_coefficients: Array of the reactants' coefficients.
        :param moles: Array representing the number of moles given of each molecule.
        """
        if reactants.size != moles.size:
            raise NotSameSizeException(f"Size of reactants array ({reactants.size}) "
                                       f"does not match that of the moles array "
                                       f"({moles.size}).")
        if reactants.size != reactant_coefficients.size:
            raise NotSameSizeException(f"Size of reactants array ({reactants.size}) "
                                       f"does not match that of the coefficients array "
                                       f"({reactant_coefficients.size}).")

        division = moles / reactant_coefficients
        index_of_minimum = np.argmin(division)

        return reactants[index_of_minimum].item()

    @classmethod
    def limiting_reactant_moles_without_coefficients(cls, reactants: ndarray[Molecule],
                                                     products: ndarray[Molecule],
                                                     moles: ndarray[float]) -> Molecule:
        """
        Returns the limiting reactant in a chemical equation without coefficients with moles.

        :param reactants: Array of the reactants represented as Molecules.
        :param products: Array of the products represented as Molecules.
        :param moles: Array representing the number of moles given of each molecule.
        """
        reactant_coefficients = Balancer.balance_equation(reactants, products)[:reactants.size]
        return Stoichiometry.limiting_reactant_moles(
            reactants,
            reactant_coefficients,
            moles
        )

    @classmethod
    def limiting_reactant_grams(cls, reactants: ndarray[Molecule],
                                reactant_coefficients: ndarray[float],
                                grams: ndarray[float]) -> Molecule:
        """
        Returns the limiting reactant in a chemical equation with coefficients with grams.

        :param reactants: Array of the reactants represented as Molecules.
        :param reactant_coefficients: Array of the reactants' coefficients.
        :param grams: Array representing the number of grams given of each molecule.
        """
        moles = Stoichiometry.convert_grams_to_moles(reactants, grams)
        return Stoichiometry.limiting_reactant_moles(
            reactants,
            reactant_coefficients,
            moles)

    @classmethod
    def limiting_reactant_grams_without_coefficients(cls, reactants: ndarray[Molecule],
                                                     products: ndarray[Molecule],
                                                     grams: ndarray[float]):
        """
        Returns the limiting reactant in a chemical equation without coefficients with grams.

        :param reactants: Array of the reactants represented as Molecules.
        :param products: Array of the products represented as Molecules.
        :param grams: Array representing the number of grams given of each molecule.
        """
        moles = Stoichiometry.convert_grams_to_moles(reactants, grams)
        return Stoichiometry.limiting_reactant_moles_without_coefficients(
            reactants,
            products,
            moles
        )
