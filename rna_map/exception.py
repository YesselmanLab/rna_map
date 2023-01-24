"""
defines exceptions for rna_map package
"""


class DREEMException(Exception):
    """
    Base class for DREEM exceptions
    """


class DREEMInputException(DREEMException):
    """
    Exception raised for errors in the input.
    """


class DREEMMissingRequirementsException(DREEMInputException):
    """
    Exception raised for missing requirements
    """

class DREEMExternalProgramException(DREEMException):
    """
    Exception raised for errors in external programs
    """
