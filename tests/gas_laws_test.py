import unittest
from src.chemiq_MYUSERNAMEHERE._operations.gas_laws import GasLaws


class GasLawsTest(unittest.TestCase):

    def test_boyle_law(self):
        answer: float = round(float(25/3), 3)
        solution: float = GasLaws.solve_boyle(P1=5, V1=5, P2=3)
        self.assertAlmostEqual(solution, answer, places=3)

    def test_charles_law(self):
        answer: float = 3
        solution: float = GasLaws.solve_charles(V1=20, T1=60, V2=1)
        self.assertAlmostEqual(solution, answer, places=3)

    def test_gay_lussac_law(self):
        answer: float = 13
        solution: float = GasLaws.solve_gay_lussac(P1=10, T1=260, P2=0.5)
        self.assertAlmostEqual(solution, answer, places=3)

    def test_combined_gas_law(self):
        answer: float = 15
        solution: float = GasLaws.solve_combined_gas_law(P1=1, V1=1, T1=1, P2=5, V2=3)
        self.assertAlmostEqual(solution, answer, places=3)

    def test_ideal_gas_law(self):
        answer: float = 17.5
        solution: float = GasLaws.solve_ideal_gas_law(P=5, V=7, N=1, R=2)
        self.assertAlmostEqual(solution, answer, places=3)
