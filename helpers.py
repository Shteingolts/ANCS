"""
A number of helper functions not relevant to the main logic.
"""

def add_spaces(string: str, width: int, indent: str = "right") -> str:
    """
    If the string is longer than provided width,
    returns the original string without change.
    """
    if width <= len(string):
        return string
    spaces_to_add = (width - len(string)) * " "
    if indent == "right":
        return spaces_to_add + string
    if indent == "left":
        return string + spaces_to_add

def table_row(items: list, widths: list, indent: str = "right") -> str:
    """
    Creates a string with the certain number of spaces between words
    alligned to either right or left
    """
    line = []
    for item, width in zip(items, widths):
        line.append(add_spaces(str(item), width, indent))

    return "".join(line) + "\n"
