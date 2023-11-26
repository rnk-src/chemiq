import numpy as np

from numpy import ndarray
from _operations.balancer import Balancer
from _operations.molecule import Molecule


def balance(reactants: ndarray[str], products: ndarray[str]) -> ndarray[int]:
    reactants: ndarray[Molecule] = np.array([Molecule.parse(reactant) for reactant in reactants])
    products: ndarray[Molecule] = np.array([Molecule.parse(product) for product in products])
    return Balancer.balance_equation(
        reactants_molecules=reactants,
        products_molecules=products
    )


def balance_molecules(reactants: ndarray[Molecule], products: ndarray[Molecule]) -> ndarray[int]:
    return Balancer.balance_equation(
        reactants_molecules=reactants,
        products_molecules=products
    )


