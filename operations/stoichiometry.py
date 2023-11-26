import numpy as np

from numpy import ndarray
from _operations.calculator import Calculator
from _operations.molecule import Molecule


def limiting_reactant_moles_coefficients(reactants: ndarray[str], coefficients: ndarray[str], moles: ndarray[str]):
    coefficients: ndarray[int] = np.array(coefficients)
    reactants: ndarray[Molecule] = np.array([Molecule.parse(reactant) for reactant in reactants])
    moles: ndarray[float] = np.array(moles).astype(dtype=float)
    solution = Calculator.limiting_reactant_moles(
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


def limiting_reactant_grams_coefficients(reactants: ndarray[str], coefficients: ndarray[str], grams: ndarray[str]):
    coefficients: ndarray[int] = np.array(coefficients)
    reactants: ndarray[Molecule] = np.array([Molecule.parse(reactant) for reactant in reactants])
    grams: ndarray[float] = np.array(grams).astype(dtype=float)
    solution = Calculator.limiting_reactant_grams(
        reactants=reactants,
        reactant_coefficients=coefficients,
        grams=grams
    )
    return {
        'reactants': str(reactants),
        'coefficients': str(coefficients),
        'grams': str(grams),
        'moles': str(Calculator.convert_grams_to_moles(reactants, grams)),
        'solution': str(solution)
    }


def limiting_reactant_moles_no_coefficients(reactants: ndarray[str], products: ndarray[str], moles: ndarray[str]):
    reactants: ndarray[Molecule] = np.array([Molecule.parse(reactant) for reactant in reactants])
    products: ndarray[Molecule] = np.array([Molecule.parse(product) for product in products])
    moles: ndarray[float] = np.array(moles).astype(dtype=float)
    solution = Calculator.limiting_reactant_moles_without_coefficients(
        reactants=reactants,
        products=products,
        moles=moles
    )
    return {
        'reactants': str(reactants),
        'products': str(products),
        'coefficients': str(Calculator.balance_equation(reactants, products)),
        'moles': str(moles),
        'solution': str(solution)
    }


def limiting_reactant_grams_no_coefficients(reactants: ndarray[str], products: ndarray[str], grams: ndarray[str]):
    reactants: ndarray[Molecule] = np.array([Molecule.parse(reactant) for reactant in reactants])
    products: ndarray[Molecule] = np.array([Molecule.parse(product) for product in products])
    grams: ndarray[float] = np.array(grams).astype(dtype=float)
    solution = Calculator.limiting_reactant_grams_without_coefficients(
        reactants=reactants,
        products=products,
        grams=grams
    )
    return {
        'reactants': str(reactants),
        'product': str(products),
        'coefficients': str(Calculator.balance_equation(reactants, products)),
        'grams': str(grams),
        'moles': str(Calculator.convert_grams_to_moles(reactants, grams)),
        'solution': str(solution)
    }
