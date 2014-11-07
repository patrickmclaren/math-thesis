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

    def test_matrix_subtraction(self):
        m_1 = mmatrix.MMatrix(matrix([[5, 6], [7, 8]]))
        m_2 = mmatrix.MMatrix(matrix([[1, 2], [3, 4]]))

        result = m_1 - m_2
        expected = matrix([[4, 4], [4, 4]])

        self.assertEqual(result.matrix, expected)

    def test_matrix_multiplication(self):
        m_1 = mmatrix.MMatrix(matrix([[1, 2], [3, 4]]))
        m_2 = mmatrix.MMatrix(matrix([[5, 6], [7, 8]]))

        result = m_1 * m_2
        expected = matrix([[19, 22], [43, 50]])

        self.assertEqual(result.matrix, expected)

    def test_matrix_get_item(self):
        m_1 = mmatrix.MMatrix(matrix([[1, 2], [3, 4]]))

        row = list(m_1[0])
        expected = [1, 2]

        self.assertEqual(row, expected)

    def test_matrix_set_item(self):
        m_1 = mmatrix.MMatrix(matrix([[0, 0], [3, 4]]))
        m_1[0] = [1, 2]

        expected = matrix([[1, 2], [3, 4]])

        self.assertEqual(m_1.matrix, expected)

    def test_matrix_row_size(self):
        m_1 = mmatrix.MMatrix(matrix([[1, 2, 3], [4, 5, 6]]))

        row_size = m_1.row_size()
        expected = 3

        self.assertEqual(row_size, expected)

    def test_matrix_column_size(self):
        m_1 = mmatrix.MMatrix(matrix([[1], [2], [3]]))

        column_size = m_1.column_size()
        expected = 3

        self.assertEqual(column_size, expected)

    def test_matrix_row_reduce(self):
        m_1 = mmatrix.MMatrix(matrix([[1, 2], [3, 4]]))

        result = m_1.row_reduce()[0]
        expected = matrix([[1, 0], [0, -2]])

        self.assertEqual(result.matrix, expected)

    def test_matrix_scalar_mult(self):
        m_1 = mmatrix.MMatrix(matrix([[1, 2], [3, 4]]))
        m_1.scalar_mult(2)

        expected = matrix([[2, 4], [6, 8]])

        self.assertEqual(m_1.matrix, expected)

    def test_matrix_swap_rows(self):
        m_1 = mmatrix.MMatrix(matrix([[1, 2], [3, 4]]))
        m_1.swap_rows(0, 1)

        expected = matrix([[3, 4], [1, 2]])

        self.assertEqual(m_1.matrix, expected)

    def test_matrix_subtract_row(self):
        m_1 = mmatrix.MMatrix(matrix([[1, 2], [3, 4]]))
        m_1.subtract_row(2, 0, 1)

        expected = matrix([[1, 2], [1, 0]])

        self.assertEqual(m_1.matrix, expected)

    def test_matrix_get_pivots(self):
        m_1 = mmatrix.MMatrix(matrix([[1, 2], [3, 4]]))
        reduced = m_1.row_reduce()[0]

        pivots = reduced.get_pivots()
        expected = [1, -2]

        self.assertEqual(pivots, expected)

    def test_matrix_is_reduced(self):
        m_1 = mmatrix.MMatrix(matrix([[1, 2], [3, 4]]))
        self.assertFalse(m_1._is_reduced())

        reduced = m_1.row_reduce()[0]
        self.assertTrue(reduced._is_reduced())

    def test_matrix_max_pivots(self):
        m_1 = mmatrix.MMatrix(matrix([[1, 2], [3, 4], [5, 6]]))

        max_pivots = m_1.max_pivots()
        expected = 2

        self.assertEqual(max_pivots, expected)

    def test_matrix_apply_fn(self):
        m_1 = mmatrix.MMatrix(matrix([[1, 2], [3, 4]]))

        result = m_1.apply_fn(lambda x: x + 1)
        expected = matrix([[2, 3], [4, 5]])

        self.assertEqual(result.matrix, expected)

    def test_matrix_update_entry(self):
        m_1 = mmatrix.MMatrix(matrix([[1, 2], [3, 4]]))
        m_1.update_entry(0, 0, 0)

        expected = matrix([[0, 2], [3, 4]])

        self.assertEqual(m_1.matrix, expected)

if __name__ == "__main__":
    cov = coverage.coverage(source=["polyrank"])
    cov.start()

    unittest.main(exit=False)

    cov.stop()
    cov.save()