from tkinter import Message


class InputError(Exception):
    """
    An error class to be raised when invalid input is given at any point in the pipeline,
    whether from improperly formatted files or from invalid UI selections.
    """

class InvalidPathError(Exception):
    """
    An error class for when a given path does not exist, does not conform to expected standards, etc.,
    but it is unclear whether this is because of user error or some other bug.
    """
    def __init__(self, path: str, message = None):
        self.path = path
        self.message = message

    def __str__(self):
        if self.message is None: return("Invalid path: " + self.path)
        else: return(self.message + "\n" + self.path)