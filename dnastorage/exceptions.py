
from dnastorage.util.stats import stats

class FaultTolerancePolicy(object):
    """ Describes policy for handling faults """
    def __init__(self):
        pass
    def allow(self, e):
        """ Override this method to give exception specific behavior """
        assert False and "Not implemented"
        return True
    def __str__(self):
        return "No policy specified"

class AllowAll(FaultTolerancePolicy):
    """ Describes policy for handling faults """
    def __init__(self):
        super(AllowAll,self).__init__()
        pass
    def allow(self, e):
        return True
    def __str__(self):
        return "Tolerate all faults as well as possible, attempt to handle"

class NoTolerance(FaultTolerancePolicy):
    def __init__(self):
        super(NoTolerance,self).__init__()
        pass
    def allow(self, e):
        return False
        

# All exceptions in dnastorage should inherit from this class
class DNAStorageError(Exception):
    """ Exception base class for error handling in this module. All """
    """ exceptions should derive from here. """
    def __init__(self,msg=None):
        if msg is None:
            msg = "A failure of some kind occurred in the dnastorage module."
        super(DNAStorageError,self).__init__(msg)
        stats.inc("Error")

class DNACodingError(DNAStorageError):
    """ An error occured while encoding or decoding a file """
    def __init__(self,msg=None):
        if msg is None:
            msg = "An error occurred while encoding/decoding a file."
        super(DNACodingError,self).__init__(msg)
        stats.inc("CodingError")

class DNAReedSolomonOuterCodeError(DNACodingError):
    """ An error occured while encoding or decoding a file """
    def __init__(self,msg=None):
        if msg is None:
            msg = "An error occurred while encoding/decoding a file."
        super(DNAReedSolomonOuterCodeError,self).__init__(msg)
        stats.inc("RSOuterCodeError")

class DNABlockBadIndex(DNACodingError):
    """ An error occured while decoding a block """
    def __init__(self,msg=None):
        if msg is None:
            msg = "Block had the wrong inter-block index."
        super(DNABlockBadIndex,self).__init__(msg)
        stats.inc("BlockBadIndex")

class DNABlockMissingIndex(DNACodingError):
    """ An error occured while decoding a block """
    def __init__(self,msg=None):
        if msg is None:
            msg = "Block is missing an index."
        super(DNABlockMissingIndex,self).__init__(msg)
        stats.inc("BlockMissingIndex")

class DNABlockTooLargeError(DNACodingError):
    """ An error occured while decoding a block """
    def __init__(self,msg=None):
        if msg is None:
            msg = "Block was too large."
        super(DNABlockTooLargeError,self).__init__(msg)
        stats.inc("BlockTooLargeError")


class DNABlockPayloadWrongSize(DNACodingError):
    """ An error occured while decoding a block """
    def __init__(self,msg=None):
        if msg is None:
            msg = "Strand in the block was the wrong size."
        super(DNABlockPayloadWrongSize,self).__init__(msg)
        stats.inc("BlockPayloadWrongSize")

class DNABadCodeword(DNACodingError):
    """ An error occured while decoding a codeword symbol """
    def __init__(self,msg=None):
        if msg is None:
            msg = "Strand in the block was the wrong size."
        super(DNABadCodeword,self).__init__(msg)
        stats.inc("BadCodeword")

class DNAStrandPayloadWrongSize(DNACodingError):
    """ An error occured while decoding a strand """
    def __init__(self,msg=None):
        if msg is None:
            msg = "Strand was the wrong size."
        super(DNAStrandPayloadWrongSize,self).__init__(msg)
        stats.inc("StrandPayloadWrongSize")


class DNAStrandMissingSequence(DNACodingError):
    """ Strand is missing a special sequence (like a cut) """
    def __init__(self,msg=None):
        if msg is None:
            msg = "Strand was missing a special sequence."
        super(DNAStrandMissingSequence,self).__init__(msg)
        stats.inc("StrandMissingSequence")

class DNAMissingPrimer(DNACodingError):
    """ Strand is missing a primer """
    def __init__(self,msg=None):
        if msg is None:
            msg = "Strand was missing a primer."
        super(DNAMissingPrimer,self).__init__(msg)
        stats.inc("MissingPrimer")


class DNAStrandPoorlyFormed(DNACodingError):
    """ Strand is poorly formed """
    def __init__(self,msg=None):
        if msg is None:
            msg = "Strand is poorly formed."
        super(DNAStrandPoorlyFormed,self).__init__(msg)
        stats.inc("StrandPoorlyFormed")

class DNAFileHeaderHasError(DNACodingError):
    """ Header is poorly formed """
    def __init__(self,msg=None):
        if msg is None:
            msg = "DNAFile header has an error."
        super(DNAFileHeaderHasError,self).__init__(msg)
        stats.inc("HeaderHasError")
        
