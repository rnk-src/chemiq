"""
This module contains the method to balance chemical equations.
"""

from fractions import Fraction

import numpy as np

from numpy import ndarray
from src.chemiq_MYUSERNAMEHERE.base.element import Element
from src.chemiq_MYUSERNAMEHERE.base.molecule import Molecule


class Balancer:  # pylint: disable=too-few-public-methods
    """
    The Balancer class is used in order to balance chemical equations.
    Balancer is intended to be used in solely static contexts and does not require the creation
    of a Balancer object in order to use the balance method.
    """

    @classmethod
    def balance_equation(cls, reactants_molecules: ndarray[Molecule],
                         products_molecules: ndarray[Molecule]) -> ndarray[int]:
        """
        Balances a chemical equation by returning the correct coefficients for the molecules in the
        reaction.
        All the corresponding coefficients are in the same order as the molecules passed in.

        :param reactants_molecules: A ndarray of reactants in Molecules.
        :param products_molecules: A ndarray of products in Molecules.
        :return: The coefficients of the reactants and products in the same array,
         with the coefficient of reactants first.
        """

        element_dictionary = {}

        elements_in_reaction = Balancer._get_elements_in_reaction(
            reactants_molecules,
            products_molecules
        )

        Balancer._add_to_element_dictionary(
            reactants_molecules,
            True,
            element_dictionary,
            elements_in_reaction
        )
        Balancer._add_to_element_dictionary(
            products_molecules,
            False,
            element_dictionary,
            elements_in_reaction
        )

        lhs = np.array([value for element, value in element_dictionary.items()])

        rhs = np.array([row[0] for row in lhs])
        lhs = np.array([row[1:] for row in lhs])
        rhs = -rhs

        solution = np.linalg.lstsq(a=lhs, b=rhs, rcond=None)[0]
        solution = np.insert(solution, 0, 1)

        return Balancer._change_float_array_to_int_proportionally(solution)

    @classmethod
    def _add_to_element_dictionary(cls, molecules: ndarray[Molecule], is_reactant: bool,
                                   element_dictionary: dict[str, ndarray[int]],
                                   elements_in_reaction: set[Element]):
        """
        Adds elements to the element dictionary along with a ndarray of their count for each
        of the molecules.

        :param molecules: A ndarray of the molecule(s) to be processed.
        :param is_reactant: If the molecules passed in are reactants, set to True, else if products,
        set to False.
        :param element_dictionary: An empty dictionary.
        :param elements_in_reaction: Unique elements present in the reaction.
        :return: None
        """

        multiplier = 1 if is_reactant else -1

        for molecule in molecules:
            Balancer._get_element_counts(molecule, element_dictionary,
                                         multiplier, elements_in_reaction)

    @classmethod
    def _get_element_counts(cls, molecule: Molecule, element_dictionary: dict[str, ndarray[int]],
                            multiplier: int, elements_in_reaction: set[Element]):
        """
        Updates a given element dictionary with coefficient values corresponding to elements and
        molecules in a chemical
        equation.

        :param molecule: The Molecule of which elements are to be counted.
        :param element_dictionary: A given dictionary to add elements and their counts in.
        :param multiplier: Give 1 if the molecule given is a reactant, else -1 if it is a product.
        :param elements_in_reaction: A set of unique elements in the reaction.
        """

        for element in elements_in_reaction:
            element_symbol = element.symbol
            count = molecule.element_counts.get(element_symbol, 0) * multiplier
            if element_dictionary.get(element_symbol, None) is None:
                element_dictionary[element_symbol] = np.array([])
            element_array = element_dictionary[element_symbol]
            element_dictionary[element_symbol] = np.append(
                arr=element_array,
                values=np.array([count])
            )

    @classmethod
    def _get_elements_in_reaction(cls, reactants: ndarray[Molecule],
                                  products: ndarray[Molecule]) -> set[Element]:
        """
        :param reactants: A ndarray of the reactant Molecules in the reaction.
        :param products: A ndarray of the product Molecules in the reaction.
        :return: A set of unique Elements in the reaction.
        """
        all_elements = set()
        for molecule in reactants:
            all_elements.update(molecule.get_all_elements())
        for molecule in products:
            all_elements.update(molecule.get_all_elements())
        return all_elements

    @classmethod
    def _change_float_array_to_int_proportionally(cls, array: ndarray[float]) -> ndarray[int]:
        """
        Turns array of floats into integers proportionally.

        :param array: Array of floats to turn into integers.
        :return: Integer numpy array with the smallest proportional numbers that aren't floats.
        """

        max_denominator = 100

        fractions = [Fraction(val).limit_denominator(max_denominator)
                     for val in array]
        ratios = np.array([(f.numerator, f.denominator) for f in fractions])

        factor = np.lcm.reduce(ratios[:, 1])
        return np.array([round(v * factor) for v in array])
