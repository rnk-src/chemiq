"""
This module contains Element class.
"""

import json
import os
from typing import Dict


class Element:
    """
    This class is responsible for handling elements and their attributes/basic operations.
    """

    elements = {}

    @classmethod
    def load_data(cls):
        """
        Loads the data of the elements from the elements.json file.
        """
        current_dir = os.path.dirname(__file__)
        elements_data_file_path = os.path.join(current_dir, '../_element_data/elements.json')
        with open(file=elements_data_file_path, mode='r', encoding='utf-8') as elements_data_file:
            data: Dict = json.load(elements_data_file)
            for element_name, element_data in data.items():
                element_data['name'] = element_name
                cls.elements[element_data['symbol']] = Element(element_data)

    @classmethod
    def has_element_for_symbol(cls, element_symbol: str) -> bool:
        """
        Returns whether the element, given a symbol, is existent.

        :param element_symbol: String representation of the element symbol.
        :return: True if the element exists, False otherwise.
        """
        return element_symbol in cls.elements

    @classmethod
    def get_element_for_symbol(cls, element_symbol: str):
        """
        Returns the element given a symbol as a string.

        :param element_symbol: String representation of the element symbol.
        :return: Element object representing the element.
        """
        if cls.has_element_for_symbol(element_symbol):
            return cls.elements[element_symbol]
        raise ValueError(f"Unknown element: {element_symbol}")

    def __init__(self, element_data: Dict):
        self.data: Dict = element_data
        self.symbol: str = element_data['symbol']
        self.atomic_mass: int = element_data['atomic_mass']

    def __hash__(self):
        return self.symbol.__hash__()


Element.load_data()
