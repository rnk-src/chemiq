import json
import os


class Element:
    elements = {}

    @classmethod
    def load_data(cls):
        current_dir = os.path.dirname(__file__)
        elements_data_file_path = os.path.join(current_dir, '../_element_data/elements.json')
        with open(elements_data_file_path, 'r') as elements_data_file:
            data: dict = json.load(elements_data_file)
            for element_name, element_data in data.items():
                element_data['name'] = element_name
                cls.elements[element_data['symbol']] = Element(element_data)

    @classmethod
    def has_element_for_symbol(cls, element_symbol: str) -> bool:
        return element_symbol in cls.elements.keys()

    @classmethod
    def get_element_for_symbol(cls, element_symbol: str):
        if cls.has_element_for_symbol(element_symbol):
            return cls.elements[element_symbol]
        else:
            raise ValueError(f"Unknown element: {element_symbol}")

    def __init__(self, element_data: dict):
        self.data = element_data
        self.symbol: str = element_data['symbol']
        self.atomic_mass: int = element_data['atomic_mass']

    def __hash__(self):
        return self.symbol.__hash__()


Element.load_data()
