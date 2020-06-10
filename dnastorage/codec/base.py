from dnastorage.codec.base_conversion import convertIntToBytes,convertBytesToInt
from dnastorage.exceptions import *

class BaseCodec(object):
    def __init__(self,CodecObj=None,Policy=None):
        self._Obj = CodecObj
        if Policy is None:
            self._Policy = AllowAll()
        else:
            self._Policy = Policy
    def _encode(self, s):
        return s
    def encode(self,s):
        if self._Obj != None:
            return self._encode(self._Obj.encode(s))
        else:
            return self._encode(s)
    def _decode(self, s):
        return s
    def decode(self,s):
        s = self._decode(s)
        if self._Obj != None:
            return self._Obj.decode(s)
        else:
            return s
