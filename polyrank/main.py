"""Determine the conditions under which a polynomial matrix has rank
less than, or equal to, n."""

# pylint: disable=import-error,wildcard-import
from sage.all import *
from polyrank.mmatrix import MMatrix

import polyrank.utils

##############################################################################
# main.py
##############################################################################

# pylint: disable=too-few-public-methods
class RowReduction:
    """
    Represents a leaf in the recursive row reduction of a matrix.
    """

    def __init__(self, matrix, ideal):
        self.matrix = matrix
        self.ideal = ideal

    def __repr__(self):
        return "Row reduced {}x{} matrix under {}".format(
            self.matrix.row_size(), self.matrix.column_size(), repr(self.ideal)
        )

    def rank(self):
        """
        Return the rank of a row reduction operation.

        :returns type: int
        """
        return self.matrix.rank()

# pylint: disable=dangerous-default-value
def recursive_reduce(mat, ring, _current_branch=None, _leaves=[]):
    """
    Recursively row reduce a matrix with entries in a given polynomial ring.

    Note that since we are not working with entries in a field, pivoting, i.e.
    multiplication by multiplicative inverses isn't well defined.

    There are various methods which can be used to ensure that rows are
    linearly independent as from prior rows, as the matrix rows are
    enumerated. However, such methods may not be valid for all points. An
    entry `a_{i,j}` may be equal to zero on some point `s`. To account for
    this, a new row-reduction is performed for every choice of pivot. The
    row-reduced matrices are returned with entries in normal form modulo a
    given ideal `I`.

    :param mat: The polynomial matrix to row reduce
    :type mat: :class:`.mmatrix.MMatrix`

    :param ring: The ring to which entries of ``mat`` belong
    :type ring: :class:`sage.rings.polynomial.multi_polynomial_ring_generic.\
MPolynomialRing_generic`

    :returns: [RowReduction, [[RowReduction, []], [RowReduction, []]]]
    :returns type: list
    """

    if _current_branch == None:
        _current_branch = ring.ideal()

    token = repr([_current_branch.groebner_basis(), mat.matrix])
    if token in _leaves:
        return []
    else:
        _leaves.append(token)

    reduced_mat, new_branches = mat.row_reduce(_current_branch)

    branched_results = [RowReduction(reduced_mat, _current_branch)]
    for branch in new_branches:
        next_branch = _current_branch + ring.ideal(branch)

        branched_results.append(
            recursive_reduce(mat, ring, next_branch, _leaves)
        )

    return branched_results

def find_conditions_rank_leq_n(mat, ring, max_rank):
    """
    Find the conditions under which a polynomial matrix, ``mat``, has rank
    less than or equal to ``max_rank``.

    :param mat: The polynomial matrix to evaluate
    :type mat: :class:`.mmatrix.MMatrix`

    :param ring: The ring to which entries of ``mat`` belong
    :type ring: :class:`sage.rings.polynomial.multi_polynomial_ring_generic.\
MPolynomialRing_generic`

    :param int max_rank: The maximum rank to check

    :returns: { 0 : [RowReduction, ...], 1 : [...], ..., n : [...] }
    :returns type: dict

    """
    results_tree = recursive_reduce(mat, ring)
    results = [x for x in polyrank.utils.flatten(results_tree)]

    output = {}
    for i in range(max_rank +1):
        output[i] = [result for result in results if result.rank() == i]

    return output

# pylint: disable=invalid-name,undefined-variable
if __name__ == "__main__":
    R = PolynomialRing(QQ, ['x%s'%i for i in range(1, 10)])
    R.inject_variables(verbose=False)

    mat_3x3 = MMatrix([[x1, x2, x3], [x4, x5, x6], [x7, x8, x9]])
    rank_n = 3

    conditions = find_conditions_rank_leq_n(mat_3x3, R, rank_n)
    for rank in conditions.keys():
        print("Rank {}:".format(rank))

        for condition in conditions[rank]:
            print(" "*4 + str(condition.ideal))
