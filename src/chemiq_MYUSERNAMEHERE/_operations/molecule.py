import chemparse
from element import Element


class Molecule:

    @classmethod
    def parse(cls, formula: str):
        element_counts = chemparse.parse_formula(formula)
        return Molecule(element_counts)

    def __init__(self, element_counts: dict[str, int]):
        self.element_counts = element_counts

    def has_element(self, element: Element):
        return element.symbol in self.element_counts.keys()

    def element_count(self, element: Element):
        if self.has_element(element):
            return self.__element_count(element)
        else:
            return 0

    def __element_count(self, element: Element):
        return self.element_counts[element.symbol]

    def get_all_elements(self) -> set[Element]:
        element_symbols = self.element_counts.keys()
        result = set()
        for element_symbol in element_symbols:
            result.add(Element.get_element_for_symbol(element_symbol))
        return result

    def get_molar_mass(self) -> float:
        total_mass = 0
        elements = self.get_all_elements()
        for element in elements:
            total_mass += element.atomic_mass * self.element_count(element)

        return total_mass

    def get_particles(self) -> dict:
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
