class MissingMSISeqError(Exception):
    """
    An error class to be raised when MSIseq functionality is invoked, but the package is not installed.
    """


class MissingDeconstructSigsError(Exception):
    """
    An error class to be raised when deconstructSigs functionality is invoked, but the package is not installed.
    """
