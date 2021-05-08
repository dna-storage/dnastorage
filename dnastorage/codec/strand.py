from random import randint

import dnastorage.exceptions as err
from dnastorage.codec.base_codec import BaseCodec
from dnastorage.codec.reedsolomon.rs import ReedSolomon,get_reed_solomon,ReedSolomonError
from dnastorage.util.stats import stats

import logging
logger = logging.getLogger('dna.storage.codec.strand')
logger.addHandler(logging.NullHandler())

class RandomizeCodec(BaseCodec):
    '''Add numRandBytes random values to the end of the strand. This is a naive way 
       of ensuring distinct strands during voting.
    '''
    def __init__(self,numRandBytes,CodecObj=None,Policy=None):
        super(RandomizeCodec,self).__init__(CodecObj=CodecObj,Policy=Policy)
        self._numRandBytes = numRandBytes

    def _encode(self, array):
        for i in range(self._numRandBytes):
            array.append( randint(0,255) )

    def _decode(self, array):
        return array[:self._numRandBytes]

class ReedSolomonInnerCodec(BaseCodec):
    """
    ReedSolomonInnerCodec takes a sequence of bytes as input to the _encode function and
    produces a Reed-Solomon encoded message as a byte array. 

    This is an 'Inner' Codec because it only can correct errors within a strand.
    """
    def __init__(self,numberECCBytes,c_exp=8,CodecObj=None,Policy=None):
        """
        numberECCBytes is the amount of error correction symbols to add. 
        c_exp is the characteristic.
        CodecObj is a nested codec, called after encoding or before decoding. See BaseCodec 
        for more info.
        Policy defines the fault handling support. By default, errors trigger an exception.        
        """
        super(ReedSolomonInnerCodec,self).__init__(CodecObj=CodecObj,Policy=Policy)

        self.rs = get_reed_solomon(c_exp=c_exp)
        self._numberECCBytes = numberECCBytes

    def _encode(self,array):
        """ Accepts list/array of bytes. Return the Reed-Solomon encoded byte array.
        """
        try:
            assert len(array) <= self.rs.field_charac
        except:
            #print "Failed RS check"
            return array
        # usually not needed, but makes the code a little more tolerant for use with
        # a variety of codecs that may pass in strings, bytearrays, or lists:
        message = [x for x in array]

        try:
            # encoded the message using the RS library
            mesecc = self.rs.rs_encode_msg(message, self._numberECCBytes)
        except ReedSolomonError as e:
            raise err.DNACodingError("Error while encoding Reed-Solomon Inner Codec.")
        except ZeroDivisionError as e:
            pass
        
        #print "RSInner:",len(mesecc)
        return mesecc

    def _decode(self,array):
        """
        This function expects a list of unsigned integers in the GF. For now, erasures are
        denoted with -1.
        """
        message = [x for x in array] 
        # find the -1s in the list
        erasures = [ i for i in range(len(message)) if message[i]==-1 ]
        # correct the message
        try:
            corrected_message, corrected_ecc = self.rs.rs_correct_msg(message,self._numberECCBytes, erase_pos=erasures)
            value = corrected_message
            #print "corrected message"
            stats.inc("RSInnerCodec::decode::succeeded")
        except ReedSolomonError as e:
            stats.inc("RSInnerCodec::decode::failed")
            #print "Inner: Couldn't correct message: {}".format(message)
            stats.inc("RSInnerCodec.ReedSolomonError")
            if self._Policy.allow(e):
                # leave erasures, may be helpful for outer decoder
                #value = message[0:(self._numberECCBytes)]
                value = [-1 for _ in range(len(array))]
                pass # nothing to do
            else:
                print (str(e))
                raise err.DNACodingError("RSInnerCodec failed to correct message.")
            # just proceed without further error correction
            pass
        except ZeroDivisionError as e:
            stats.inc("RSInnerCodec.ZeroDivision")
            if self._Policy.allow(e):
                # leave erasures, may be helpful for outer decoder
                #value = message[0:(self._numberECCBytes)]
                value = [-1 for _ in range(len(array))]                
                pass # nothing to do
            else:
                print (str(e))
                raise err.DNACodingError("RSInnerCodec failed to correct message.")

        return value
