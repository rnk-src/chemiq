"""
Offers balancing operations.
"""

from numpy import ndarray
import numpy as np

from src.chemiq_MYUSERNAMEHERE.base.balancer import Balancer
from src.chemiq_MYUSERNAMEHERE.base.molecule import Molecule


def balance(reactants: ndarray[str], products: ndarray[str]) -> ndarray[int]:
    """
    Balances an equation given reactants and products that can be parsed and returns the correct
     coefficients.

    :param reactants: A string list of the reactants in the chemical equation.
    :param products: A string list of the products in the chemical equation.
    :return: An int list of the correct coefficients.
    """
    reactants: ndarray[Molecule] = np.array([Molecule.parse(reactant) for reactant in reactants])
    products: ndarray[Molecule] = np.array([Molecule.parse(product) for product in products])
    return Balancer.balance_equation(
        reactants_molecules=reactants,
        products_molecules=products
    )


def balance_molecules(reactants: ndarray[Molecule], products: ndarray[Molecule]) -> ndarray[int]:
    """
    Balances an equation given reactants and products and returns the correct coefficients.

    :param reactants: A list of the reactant Molecules in the chemical equation.
    :param products: A list of the product Molecules in the chemical equation.
    :return: An int list of the correct coefficients.
    """
    return Balancer.balance_equation(
        reactants_molecules=reactants,
        products_molecules=products
    )
