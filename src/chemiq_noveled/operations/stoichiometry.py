"""
Offers methods to aid with Stoichiometry calculations.
"""

import numpy as np
from numpy import ndarray

from src.chemiq_noveled.base.stoichiometry import Stoichiometry
from src.chemiq_noveled.base.balancer import Balancer
from src.chemiq_noveled.base.molecule import Molecule


def limiting_reactant_moles_coefficients(
        reactants: ndarray[str],
        coefficients: ndarray[str],
        moles: ndarray[str]):
    """
    Determines the limiting reactant in a reaction given coefficients, reactants, and moles of the
    reactants.

    :param reactants: The reactants in the chemical equation.
    :param coefficients: The coefficients of the reactants in the chemical equation.
    :param moles: Moles of the reactants in the chemical equation.
    """
    coefficients: ndarray[int] = np.array(coefficients)
    reactants: ndarray[Molecule] = np.array([Molecule.parse(reactant) for reactant in reactants])
    moles: ndarray[float] = np.array(moles).astype(dtype=float)
    solution = Stoichiometry.limiting_reactant_moles(
        reactants=reactants,
        reactant_coefficients=coefficients,
        moles=moles
    )
    return {
        'reactants': str(reactants),
        'coefficients': str(coefficients),
        'moles': str(moles),
        'solution': str(solution)
    }


def limiting_reactant_grams_coefficients(
        reactants: ndarray[str],
        coefficients: ndarray[str],
        grams: ndarray[str]):
    """
    Determines the limiting reactant in a reaction given coefficients, reactants, and moles of the
    reactants.

    :param reactants: The reactants in the chemical equation.
    :param coefficients: The coefficients of the reactants in the chemical equation.
    :param grams: Grams of the reactants in the chemical equation.
    """
    coefficients: ndarray[int] = np.array(coefficients)
    reactants: ndarray[Molecule] = np.array([Molecule.parse(reactant) for reactant in reactants])
    grams: ndarray[float] = np.array(grams).astype(dtype=float)
    solution = Stoichiometry.limiting_reactant_grams(
        reactants=reactants,
        reactant_coefficients=coefficients,
        grams=grams
    )
    return {
        'reactants': str(reactants),
        'coefficients': str(coefficients),
        'grams': str(grams),
        'moles': str(Stoichiometry.convert_grams_to_moles(reactants, grams)),
        'solution': str(solution)
    }


def limiting_reactant_moles_no_coefficients(
        reactants: ndarray[str],
        products: ndarray[str],
        moles: ndarray[str]):
    """
    Determines the limiting reactant in a reaction given reactants, products, and moles of the
    reactants.

    :param reactants: The reactants in the chemical equation.
    :param products: The products in the chemical equation.
    :param moles: Moles of the reactants in the chemical equation.
    """
    reactants: ndarray[Molecule] = np.array([Molecule.parse(reactant) for reactant in reactants])
    products: ndarray[Molecule] = np.array([Molecule.parse(product) for product in products])
    moles: ndarray[float] = np.array(moles).astype(dtype=float)
    solution = Stoichiometry.limiting_reactant_moles_without_coefficients(
        reactants=reactants,
        products=products,
        moles=moles
    )
    return {
        'reactants': str(reactants),
        'products': str(products),
        'coefficients': str(Balancer.balance_equation(reactants, products)),
        'moles': str(moles),
        'solution': str(solution)
    }


def limiting_reactant_grams_no_coefficients(
        reactants: ndarray[str],
        products: ndarray[str],
        grams: ndarray[str]):
    """
    Determines the limiting reactant in a reaction given reactants, products, and grams of the
    reactants.

    :param reactants: The reactants in the chemical equation.
    :param products: The products in the chemical equation.
    :param grams: Grams of the reactants in the chemical equation.
    """
    reactants: ndarray[Molecule] = np.array([Molecule.parse(reactant) for reactant in reactants])
    products: ndarray[Molecule] = np.array([Molecule.parse(product) for product in products])
    grams: ndarray[float] = np.array(grams).astype(dtype=float)
    solution = Stoichiometry.limiting_reactant_grams_without_coefficients(
        reactants=reactants,
        products=products,
        grams=grams
    )
    return {
        'reactants': str(reactants),
        'product': str(products),
        'coefficients': str(Balancer.balance_equation(reactants, products)),
        'grams': str(grams),
        'moles': str(Stoichiometry.convert_grams_to_moles(reactants, grams)),
        'solution': str(solution)
    }
