"""
This modules provides functionality for formatting strings for logging, standard output and standard error.
"""
from preconditions import preconditions


@preconditions(lambda string: isinstance(string, str))
def pad_spaces(string, num_chars = 70):
    """
    Adds empty spaces to the string so that the final string has at least the specified number of characters.
    If the num_chars argument is less than the length of the original string, the original string is returned.
    :param string: original string
    :param num_chars: number of characters of the final padded string
    :return: the padded string
    """
    num_spaces = num_chars - len(string)
    if num_spaces <= 0:
        return string
    return string + " "*num_spaces