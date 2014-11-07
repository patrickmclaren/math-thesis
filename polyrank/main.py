from collections import defaultdict
from sage.all import *

import itertools

from polyrank.mmatrix import MMatrix

##############################################################################
# Helper Functions
##############################################################################

def _join_dicts_of_lists(*ddicts):
    """
    Merge multiple ``dict`` of ``list`` into one ``dict``.

    :param \*dict ddicts: collection to merge
    :return type: dict

    Example::
        >>> _join_dicts_of_lists({ 1 : [3]}, { 1 : [4], 2 : [5]})
        { 1 : [3, 4], 2 : [5] }
    """

    new_dict = defaultdict(list)

    for ddict in ddicts:
        for k,v in ddict.items():
            new_dict[k].extend(v)

    return new_dict

def _clean_dict_of_lists(ddict):
    """
    Returns ``ddict`` with distinct values

    :param dict ddict: ``dict`` to reduce

    Example::
        >>> _clean_dict_of_lists({ 1 : [0, 0, 1, 1] })
        { 1 : [0, 1] }
    """
    for k,v in ddict.items():
        new_v = []
        for el in v:
            if not el in new_v:
                new_v.append(el)

        ddict[k] = new_v

    return ddict

##############################################################################
# main.py
##############################################################################

def poly_div_rem(f, g):
    """
    Divide ``f`` by ``g``, and return the remainder.

    :param f: The dividend
    :type f: :class:`sage.rings.polynomial.multi_polynomial_element.MPolynomial_element`
    :param g: The divisor
    :type g: :class:`sage.rings.polynomial.multi_polynomial_element.MPolynomial_element`

    :Returns:
        f/g

    :Returns Type:
        :class:`sage.rings.polynomial.multi_polynomial_element.MPolynomial_element`

    """
    return f.quo_rem(g)[1]

def find_conditions_rank_leq_n(mat, n):
    """
    Find the conditions under which a polynomial matrix, ``mat``, has rank
    less than or equal to ``n``.

    :param mat: The polynomial matrix to evaluate
    :type mat: :class:`.mmatrix.MMatrix`
    :param int n: The maximum rank to check

    :returns: { 0 : [f_1, f_2, ..., f_n], 1 : [...], ..., n : [...] }
    :returns type: dict

    """
    conditions = defaultdict(list)

    for i in reversed(range(n + 1)):
        reduced_mat, non_zero = mat.row_reduce()

        for pivot_subset in itertools.combinations([k for k in range(reduced_mat.max_pivots())], i):
            # only consider unique, non-zero pivots
            non_zero_pivots = list(set([reduced_mat[x][x] for x in pivot_subset if reduced_mat[x][x] != 0]))

            conditions[len(non_zero_pivots)].append(non_zero_pivots)

        for el in non_zero:
            next_mat = reduced_mat.apply_fn(poly_div_rem, el)
            conditions = _join_dicts_of_lists(conditions, find_conditions_rank_leq_n(next_mat, i))

    conditions = _clean_dict_of_lists(conditions)

    return conditions

if __name__ == "__main__":
    R, (x, y) = PolynomialRing(QQ, 2, 'xy').objgens()

    f = x
    g = y
    h = x
    k = 2 * y

    mat = MMatrix([[f, g], [h, k]])
    n = 2

    print(find_conditions_rank_leq_n(mat, n))
