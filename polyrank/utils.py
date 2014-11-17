"""Miscellaneous Functional Operations"""

import collections

def _join_dicts_of_lists(*ddicts):
    r"""
    Merge multiple ``dict`` of ``list`` into one ``dict``.

    :param \*dict ddicts: collection to merge
    :return type: dict

    Example::
        >>> _join_dicts_of_lists({ 1 : [3]}, { 1 : [4], 2 : [5]})
        { 1 : [3, 4], 2 : [5] }
    """

    new_dict = collections.defaultdict(list)

    for ddict in ddicts:
        for key, val in ddict.items():
            new_dict[key].extend(val)

    return new_dict

def _clean_dict_of_lists(ddict):
    """
    Returns ``ddict`` with distinct values

    :param dict ddict: ``dict`` to reduce

    Example::
        >>> _clean_dict_of_lists({ 1 : [0, 0, 1, 1] })
        { 1 : [0, 1] }
    """
    for key, val in ddict.items():
        new_v = []
        for elm in val:
            if not elm in new_v:
                new_v.append(elm)

        ddict[key] = new_v

    return ddict

def flatten(tree):
    """
    Flattens irregularly nested lists.

    :param list tree: irregular list of lists

    :return type: generator
    """

    for elm in tree:
        if isinstance(elm, collections.Iterable) and \
            not isinstance(elm, basestring): # pylint: disable=undefined-variable
            for sub in flatten(elm):
                yield sub
        else:
            yield elm
