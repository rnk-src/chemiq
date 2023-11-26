import unittest
import numpy as np

from _operations.balancer import Balancer
from _operations.molecule import Molecule


class TestBalancer(unittest.TestCase):

    def test_balance_equation(self):
        """ Previously Tested Equations:
        Ba(OH)2 + H3PO4 = Ba3(PO4)2 + H2O [3, 2, 1, 6] #1
        H2 + O2 -> H2O [2, 1, 2] #2
        C4H10 + O2 -> CO2 + H2O [2, 13, 8, 10] #3
        Ga + CuBr2 -> GaBr3 + Cu [2, 3, 2, 3] #4
        I2 + F2 -> IF7 [1, 7, 2] #5
        HCl + Ca(OH)2 -> CaCl2 + H2O [2, 1, 1, 2] #6
        C4O6H4 + O2 -> CO2 + H2O [1, 2, 4, 2] #7
        """
        reactant_molecules = np.array([Molecule({'C': 4, 'O': 6, 'H': 4}), Molecule({'O': 2})])
        product_molecules = np.array([Molecule({'C': 1, 'O': 2}), Molecule({'H': 2, 'O': 1})])
        result = Balancer.balance_equation(
            reactants_molecules=reactant_molecules,
            products_molecules=product_molecules
        )
        actual = np.array([1, 2, 4, 2])
        self.assertTrue(np.array_equal(result, actual))

    def test_add_to_element_dictionary_reactant(self):
        reactant_molecules = np.array([Molecule({'C': 2, 'H': 6}), Molecule({'O': 2})])
        product_molecules = np.array([Molecule({'C': 1, 'O': 2}), Molecule({'H': 2, 'O': 1})])
        elements_in_reaction = Balancer._get_elements_in_reaction(reactant_molecules, product_molecules)
        element_dictionary = dict()
        Balancer._add_to_element_dictionary(
            molecules=reactant_molecules,
            is_reactant=True,
            element_dictionary=element_dictionary,
            elements_in_reaction=elements_in_reaction
        )
        self.assertIsNone(np.testing.assert_equal(element_dictionary, {'C': np.array([2, 0]), 'H': np.array([6, 0]),
                                                                       'O': np.array([0, 2])}))

    def test_add_to_element_dictionary_reactant_and_product(self):
        reactant_molecules = np.array([Molecule({'C': 2, 'H': 6}), Molecule({'O': 2})])
        product_molecules = np.array([Molecule({'C': 1, 'O': 2}), Molecule({'H': 2, 'O': 1})])
        elements_in_reaction = Balancer._get_elements_in_reaction(reactant_molecules, product_molecules)
        element_dictionary = dict()
        Balancer._add_to_element_dictionary(
            molecules=reactant_molecules,
            is_reactant=True,
            element_dictionary=element_dictionary,
            elements_in_reaction=elements_in_reaction
        )
        Balancer._add_to_element_dictionary(
            molecules=product_molecules,
            is_reactant=False,
            element_dictionary=element_dictionary,
            elements_in_reaction=elements_in_reaction
        )
        # self.assertEquals(element_dictionary, {'C': [2, 0, -1, 0], 'H': [6, 0, 0, -2], 'O': [0, 2, -2, -1]})
        # {'C': np.array([2, 0, -1, 0]), 'H': np.array([6, 0, 0, -2]), 'O': np.array([0, 2, -2, -1])}
        self.assertIsNone(np.testing.assert_equal(element_dictionary, {'C': np.array([2, 0, -1, 0]),
                                                                       'H': np.array([6, 0, 0, -2]),
                                                                       'O': np.array([0, 2, -2, -1])}))

    def test_get_element_counts(self):
        molecule = Molecule({'C': 2, 'H': 6})
        element_dictionary = dict()
        multiplier = 1
        reactant_molecules = np.array([Molecule({'C': 2, 'H': 6}), Molecule({'O': 2})])
        product_molecules = np.array([Molecule({'C': 1, 'O': 2}), Molecule({'H': 2, 'O': 1})])
        elements_in_reaction = Balancer._get_elements_in_reaction(reactant_molecules, product_molecules)
        Balancer._get_element_counts(molecule, element_dictionary, multiplier, elements_in_reaction)
        self.assertEquals(element_dictionary, {'C': [2], 'H': [6], 'O': [0]})

    def test_get_elements_in_reaction(self):
        reactant_molecules = np.array([Molecule({'C': 2, 'H': 6}), Molecule({'O': 2})])
        product_molecules = np.array([Molecule({'C': 1, 'O': 2}), Molecule({'H': 2, 'O': 1})])
        all_elements = set()
        all_elements.update(Molecule({'C': 2, 'H': 6}).get_all_elements())
        all_elements.update(Molecule({'O': 2}).get_all_elements())
        all_elements.update(Molecule({'C': 1, 'O': 2}).get_all_elements())
        all_elements.update(Molecule({'H': 2, 'O': 1}).get_all_elements())
        self.assertEquals(Balancer._get_elements_in_reaction(reactant_molecules, product_molecules), all_elements)

    def test_change_float_array_to_int_proportionally(self):
        first = np.array([1.5, 2.5, 2, 1])
        first = Balancer._change_float_array_to_int_proportionally(first)
        self.assertTrue(np.array_equal(first, np.array([3, 5, 4, 2])))
