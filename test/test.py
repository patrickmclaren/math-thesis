from sage.all import *
from polyrank import *

import unittest
import coverage

class TestPolyRankMMatrix(unittest.TestCase):
    def test_matrix_addition(self):
        m_1 = mmatrix.MMatrix(matrix([[1, 2], [3, 4]]))
        m_2 = mmatrix.MMatrix(matrix([[5, 6], [7, 8]]))

        result = m_1 + m_2
        expected = matrix([[6, 8], [10, 12]])

        self.assertEqual(result.matrix, expected)

if __name__ == "__main__":
    cov = coverage.coverage(source=["polyrank"])
    cov.start()

    unittest.main(exit=False)

    cov.stop()
    cov.save()