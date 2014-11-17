# pylint: disable=missing-docstring

# pylint: disable=import-error,wildcard-import
#from sage.all import *

class MMatrix():
    """
    Decorator class for a Sage matrix of type :class:`sage.structure.element.M\
atrix`.

    Provides immutable infix operations by returning a new :class:`MMatrix` that
    represents the result. Specifically, the infix operations ``+``, ``-``, an\
d ``*``
    are immutable.

    Matrix rows can be retrived and set using the ``[]`` notation, note that
    setting a row is *not* immutable.

    """

    def __init__(self, mat):
        self.matrix = matrix(copy(mat))
        """
        The underlying Sage matrix - a subtype of :class:`sage.structure.eleme\
nt.Matrix`.
        """

    def __repr__(self):
        """String Representation"""
        return self.matrix.__repr__()

    def __add__(self, other):
        """Infix Addition"""
        return MMatrix(self.matrix + other.matrix)

    def __sub__(self, other):
        """Infix Subtraction"""
        return MMatrix(self.matrix - other.matrix)

    def __mul__(self, other):
        """Infix Multiplication"""
        return MMatrix(self.matrix * other.matrix)

    def __getitem__(self, i):
        """Returns a copy of row ``i``."""
        return copy(self.matrix[i])

    def __setitem__(self, i, new_row):
        """Set row ``i`` to ``new_row``."""
        self.matrix[i] = new_row

    def __deepcopy__(self):
        """Allows usage of ``copy`` on :class:`MMatrix` instance."""
        return MMatrix(self.matrix)

    def rank(self):
        """Get the rank of :attr:`.matrix`.

        :Returns Type:
            int
        """

        return self.matrix.rank()

    def row_size(self):
        """
        Get the row size (dimension of the image) of :attr:`.matrix`.

        :Returns Type:
            int

        """

        return self.matrix.dimensions()[1]

    def column_size(self):
        """
        Get the column size (dimension of the image of the transpose) of
        :attr:`.matrix`.

        :Returns Type:
            int
        """

        return self.matrix.dimensions()[0]

    def row_reduce(self, ideal=None):
        """
        Row reduce :attr:`.matrix`.

        :Returns:
            (reduced matrix, non-zero entries)

        :Returns Type:
            tuple

        """
        new_mat = self.__deepcopy__()

        if ideal:
            new_mat = self.apply_fn(ideal.reduce)

        branches = []

        for i in range(new_mat.max_pivots()):
            if new_mat[i][i] == 0:
                found_row = False
                for j in range(i+1, new_mat.row_size()):
                    if new_mat[j][i] != 0:
                        new_mat.swap_rows(i, j)
                        found_row = True
                        break
                if not found_row:
                    continue

            pivot = new_mat[i][i]

            for j in range(new_mat.row_size()):
                if i == j or new_mat[j][i] == 0:
                    continue

                target_lcm = lcm(pivot, new_mat[j][i])

                pivot_q = target_lcm / pivot
                clear_q = target_lcm / new_mat[j][i]

                new_mat.scalar_mult(clear_q, j)
                new_mat.subtract_row(pivot_q, i, j)

            branches.append(pivot)

        return (new_mat, branches)

    def scalar_mult(self, scalar, row=None):
        """
        Multiply :attr:`.matrix` by ``scalar``.

        :param float scalar: The scalar to multiply :attr:`.matrix` by

        :Returns Type:
            None

        .. warning::

           Not immutable.
        """
        if row is None:
            for i in range(self.row_size()):
                self[i] = [x * scalar for x in self[i]]
        else:
            self[row] = [x * scalar for x in self[row]]

    def swap_rows(self, i, j):
        """
        Swap row ``i`` with row ``j``.

        :param int i: The first row number
        :param int j: the second row number

        :Returns Type:
            None

        .. warning::

           Not immutable.

        """
        row_i = self[i]
        row_j = self[j]

        self[i] = row_j
        self[j] = row_i

    def subtract_row(self, mult, i, j):
        """
        Subtract ``mult`` times row ``i`` from row ``j``.

        :param float mult: The scalar to multiply row ``i``
        :param int i: The row to subtract from row ``j``
        :param int j: The row from which row ``i`` is subtracted

        :Returns Type:
            None

        .. warning::

           Not immutable.
        """
        row_i = [x * mult for x in self[i]]
        row_j = self[j]

        self[j] = [row_j[k] - row_i[k] for k in range(len(row_j))]

    def get_pivots(self):
        """
        Returns a list of the pivots of :attr:`.matrix`.

        :Returns:
            [pivot_1, pivot_2, ..., pivot_n]

        :Returns Type:
            list
        """
        if not self._is_reduced():
            raise IndexError

        return [self[i][i] for i in range(self.max_pivots())]

    def _is_reduced(self):
        """
        Returns ``True`` if :attr:`.matrix` is row reduced, and ``False`` othe\
rwise.

        :Returns Type:
            bool
        """
        for i in range(self.max_pivots()):
            for j in range(i+1, self.row_size()):
                if self[j][i] != 0:
                    return False

        return True

    def max_pivots(self):
        """
        Returns the maximum possible number of pivots of :attr:`.matrix`.

        :Returns Type:
            int
        """
        return min(self.row_size(), self.column_size())

    def apply_fn(self, func, *args, **kwargs):
        r"""
        Apply ``func`` with positional parameters ``args`` and keyword
        parameters ``kwargs`` to the entries of :attr:`.matrix`.

        :param function func: The function to apply to :attr:`.matrix`
        :param \*args args: Positional parameters to be applied to ``func``
        :param \*\*kwargs kwargs: Keyword parameters to be applied to ``func``

        :Returns Type:
            :class:`MMatrix`

        """
        new_mat = self.__deepcopy__()

        for i in range(new_mat.row_size()):
            for j in range(new_mat.column_size()):
                new_mat.update_entry(i, j, func(new_mat[i][j], *args, **kwargs))

        return new_mat

    def update_entry(self, i, j, val):
        """
        Update entry in row `i`, column `j` to `val`.

        :param int i: The row number of the entry to be updated
        :param int j: The column number of the entry to be updated
        :param Polynomial val: The new value for entry ``i, j``

        :Returns Type:
            None

        .. warning::

           Not immutable.
        """
        tmp_row = self[i]
        tmp_row[j] = val
        self[i] = tmp_row
