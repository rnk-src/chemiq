"""
The Calculator class contains functions to perform all general calculations this package offers.
"""

from numpy import ndarray

from src.chemiq_MYUSERNAMEHERE.base.element import Element
from src.chemiq_MYUSERNAMEHERE.base.molecule import Molecule
from src.chemiq_MYUSERNAMEHERE.base.balancer import Balancer
from src.chemiq_MYUSERNAMEHERE.base.stoichiometry import Stoichiometry
from src.chemiq_MYUSERNAMEHERE.base.gas_laws import GasLaws


class NotSameSizeException(Exception):
    """
    Raised when the number of reactants and moles are not the same.
    """


class Calculator:
    """
    The Calculator class is to be used in order to perform various chemical calculations.
    This includes different topics of chemistry calculation, such as balancing chemical equations,
    performing stoichiometry, and performing gas laws calculations.
    """

    @classmethod
    def get_atomic_mass(cls, element_symbol: str) -> float:
        """
        Returns the atomic mass of an element given its symbol.
        :param element_symbol: String representing the symbol of the element.
        :return: Atomic mass as a float.
        """
        element = Element.get_element_for_symbol(element_symbol)
        return element.atomic_mass

    @classmethod
    def get_element_properties(cls, element: Element) -> dict:
        """
        Returns a dictionary containing the properties of an element.
        :param element: Element object containing the data.
        :return: Dictionary containing data regarding the element.
        """
        return element.data

    @classmethod
    def molar_mass_of(cls, molecule_formula: str) -> float:
        """
        Returns the molar mass of a given molecule.
        :param molecule_formula: String representing the molecule.
        :return: Molar mass of the molecule as a float.
        """
        molecule = Molecule.parse(molecule_formula)
        return molecule.get_molar_mass()

    @classmethod
    def balance_equation(cls, reactants: ndarray[Molecule],
                         products: ndarray[Molecule]) -> ndarray[float]:
        """
        Returns the coefficients of a balanced equation given reactants and products.
        :param reactants: NumPy array containing the reactant Molecules.
        :param products: NumPy array containing the product Molecules.
        :return: NumPy array containing the coefficients of the balanced chemical reaction
         in sequential order as floats.
        """
        return Balancer.balance_equation(
            reactants_molecules=reactants,
            products_molecules=products
        )

    @classmethod
    def element_properties(cls, element_symbol: str):
        """
        Returns a dictionary containing the properties of an element.
        :param element_symbol: String object representing the element.
        :return: Dictionary containing data regarding the element.
        """
        element = Element.get_element_for_symbol(element_symbol)
        return element.data

    @classmethod
    def get_particles(cls, molecule: Molecule):
        """
        Returns the particles in a molecule.
        :param molecule: Molecule object.
        :return: A dictionary containing the values of electrons, neutrons, and protons in an atom.
        """
        return molecule.get_particles()

    @classmethod
    def moles_to_grams(cls, moles: float, formula: str):
        """
        Converts moles to grams given a formula.
        :param moles: Number of moles.
        :param formula: Formula to parse.
        :return: Grams of the given molecule.
        """
        return Stoichiometry.moles_to_grams(
            moles=moles,
            formula=formula
        )

    @classmethod
    def grams_to_moles(cls, grams: float, formula: str):
        """
        Converts grams to moles given a formula.
        :param grams: Number of grams.
        :param formula: Formula to parse.
        :return: Moles of the given molecule.
        """
        return Stoichiometry.grams_to_moles(
            grams=grams,
            formula=formula
        )

    @classmethod
    def convert_grams_to_moles(cls, reactants: ndarray[Molecule],
                               grams: ndarray[float]) -> ndarray[float]:
        """
        Convert a NumPy array of reactants from grams to moles.
        :param reactants: NumPy array of reactant Molecules.
        :param grams: NumPy array of reactant amounts in grams.
        :return: NumPy array of reactant amounts in moles.
        """
        return Stoichiometry.convert_grams_to_moles(
            molecules=reactants,
            grams=grams
        )

    @classmethod
    def limiting_reactant_moles(cls, reactants: ndarray[Molecule],
                                reactant_coefficients: ndarray[float],
                                moles: ndarray[float]) -> Molecule:
        """
        Returns the limiting reactant in a chemical equation.
        :param reactants: Reactant Molecules in a chemical equation.
        :param reactant_coefficients: NumPy array of reactant coefficients in the same size
         as the "reactants" array.
        :param moles: NumPy array of moles in the same size as the "reactants" array.
        :return: The limiting reactant Molecule.
        """
        return Stoichiometry.limiting_reactant_moles(
            reactants=reactants,
            reactant_coefficients=reactant_coefficients,
            moles=moles
        )

    @classmethod
    def limiting_reactant_moles_without_coefficients(cls, reactants: ndarray[Molecule],
                                                     products: ndarray[Molecule],
                                                     moles: ndarray[float]) -> Molecule:
        """
        Returns the limiting reactant in a chemical equation.
        :param reactants: Reactant Molecules in a chemical equation.
        :param products: Product Molecules in a chemical equation.
        :param moles: NumPy array of moles in the same size as the "reactants" array.
        :return: The limiting reactant Molecule.
        """
        return Stoichiometry.limiting_reactant_moles_without_coefficients(
            reactants=reactants,
            products=products,
            moles=moles
        )

    @classmethod
    def limiting_reactant_grams(cls, reactants: ndarray[Molecule],
                                reactant_coefficients: ndarray[float],
                                grams: ndarray[float]) -> Molecule:
        """
        Returns the limiting reactant in a chemical equation.
        :param reactants: Reactant Molecules in a chemical equation.
        :param reactant_coefficients: NumPy array of reactant coefficients in the same size as the
         "reactants" array.
        :param grams: NumPy array of grams in the same size as the "reactants" array.
        :return:
        """
        return Stoichiometry.limiting_reactant_grams(
            reactants=reactants,
            reactant_coefficients=reactant_coefficients,
            grams=grams
        )

    @classmethod
    def limiting_reactant_grams_without_coefficients(cls, reactants: ndarray[Molecule],
                                                     products: ndarray[Molecule],
                                                     grams: ndarray[float]):
        """
        Returns the limiting reactant in a chemical equation.
        :param reactants: Reactant Molecules in a chemical equation.
        :param products: Product Molecules in a chemical equation.
        :param grams: NumPy array of grams in the same size as the "reactants" array.
        :return: The limiting reactant Molecule.
        """
        return Stoichiometry.limiting_reactant_grams_without_coefficients(
            reactants=reactants,
            products=products,
            grams=grams
        )

    @classmethod
    def solve_boyle_law(cls, arguments: dict[str, float]):
        """
        Returns the solution to satisfy Boyle's Law given all but one variable.

        Example: Solve for P1 given the variables P2 = 5, V1 = 3, and V2 = 17
        arguments = {'P2': 5, 'V1': 3, 'V2': 17}
        Calculator.solve_boyle_law(arguments)

        :param arguments: Dictionary of string and float, where string indicates the variable.
        :return: The solution of the missing variable.
        """

        return GasLaws.solve_boyle(**arguments)

    @classmethod
    def solve_charles_law(cls, arguments: dict[str, float]):
        """
        Returns the solution to satisfy Charles' Law given all but one variable.

        Example: Solve for V1 given the variables V2 = 5, T1 = 3, and T2 = 17
        arguments = {'V2': 5, 'T1': 3, 'T2': 17}
        Calculator.solve_charles_law(arguments)

        :param arguments: Dictionary of string and float, where string indicates the variable.
        :return: The solution of the missing variable.
        """

        return GasLaws.solve_charles(**arguments)

    @classmethod
    def solve_gay_lussac_law(cls, arguments: dict[str, float]):
        """
        Returns the solution to satisfy Gay Lussac's Law given all but one variable.

        Example: Solve for P1 given the variables P2 = 5, T1 = 3, and T2 = 17
        arguments = {'P2': 5, 'T1': 3, 'T2': 17}
        Calculator.solve_gay_lussac_law(arguments)

        :param arguments: Dictionary of string and float, where string indicates the variable.
        :return: The solution of the missing variable.
        """

        return GasLaws.solve_gay_lussac(**arguments)

    @classmethod
    def solve_combined_gas_law(cls, arguments: dict[str, float]):
        """
        Returns the solution to satisfy the Combined Gas Law given all but one variable.
        Example: Solve for P1 given the variables P2 = 5, V1 = 18, V2 = 4, T1 = 3, and T2 = 17
        arguments = {'P2': 5, 'V1': 18, 'V2': 4, 'T1': 3, 'T2': 17}
        Calculator.solve_combined_gas_law(arguments)

        :param arguments: Dictionary of string and float, where string indicates the variable.
        :return: The solution of the missing variable.
        """

        return GasLaws.solve_combined_gas_law(**arguments)

    @classmethod
    def solve_ideal_gas_law(cls, arguments: dict[str, float]):
        """
        Returns the solution to satisfy the Ideal Gas Law given all but one variable.

        Example: Solve for P given the variables V = 5, N = 3, R = 47, and T = 17
        arguments = {'V': 5, 'N': 3, 'R': 47, 'T': 17}
        Calculator.solve_ideal_gas_law(arguments)

        :param arguments: Dictionary of string and float, where string indicates the variable.
        :return: The solution of the missing variable.
        """

        return GasLaws.solve_ideal_gas_law(**arguments)
