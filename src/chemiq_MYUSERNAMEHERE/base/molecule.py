"""
This module contains the Molecule class to aid with molecule-related operations.
"""

import chemparse
from src.chemiq_MYUSERNAMEHERE.base.element import Element


class Molecule:
    """
    This class is responsible for handling molecule related operations.
    """

    @classmethod
    def parse(cls, formula: str):
        """
        Parses a string representing of a molecule.

        :param: String representing a molecule.
        :return: A Molecule object representing the molecule formula given.
        """
        element_counts = chemparse.parse_formula(formula)
        return Molecule(element_counts)

    def __init__(self, element_counts: dict[str, int]):
        self.element_counts = element_counts

    def has_element(self, element: Element) -> bool:
        """
        Returns if the molecule has the given element.

        :param: The element to check the presence of.
        :return: A boolean value indicating whether the molecule has the given element.
        """
        return element.symbol in self.element_counts

    def element_count(self, element: Element) -> int:
        """
        Returns the count of the given element in the molecule, and returns 0 if not present.

        :param element: The element to check the count for.
        :return: The count of the given Element.
        """
        if self.has_element(element):
            return self.__element_count(element)
        return 0

    def __element_count(self, element: Element) -> int:
        """
        Returns the count of the given element in the molecule.

        :param element: The element to check the count for.
        :return: The count of the given Element.
        """
        return self.element_counts[element.symbol]

    def get_all_elements(self) -> set[Element]:
        """
        Returns all the elements in the molecule.

        :return: Set containing the elements in the molecule.
        """
        element_symbols = self.element_counts.keys()
        result = set()
        for element_symbol in element_symbols:
            result.add(Element.get_element_for_symbol(element_symbol))
        return result

    def get_molar_mass(self) -> float:
        """
        Returns the molar mass of a molecule.

        :return: The molar mass of a molecule represented as a float.
        """
        total_mass = 0
        elements = self.get_all_elements()
        for element in elements:
            total_mass += element.atomic_mass * self.element_count(element)

        return total_mass

    def get_particles(self) -> dict[str, int]:
        """
        Returns the particles (electrons, neutrons, protons) and their
         respective counts in the molecule.

        :return: A dictionary containing the type of particle (key) and their
         respective count (value)
        """
        elements = self.get_all_elements()
        particles = {'electrons': 0, 'neutrons': 0, 'protons': 0}
        for element in elements:
            particles['electrons'] += element.data['number'] * self.element_count(element)
            particles['protons'] += element.data['number'] * self.element_count(element)
            particles['neutrons'] += ((round(element.data['atomic_mass']) - element.data['number'])
                                      * self.element_count(element))
        return particles

    def __repr__(self):
        return self.element_counts.__repr__()
