from dnastorage.codec.base_codec import BaseCodec

class EncodePacketizedFile(object):
    """
    Partition the file into key,value packets for encoding. This is an iterable
    object for convenient encoding. If a CodecObj is provided, the iterator
    returns the packets after they pass through the CodecObj.

    For simple encoders and decoders that operate on one strand at a time, you can simply
    override the _encode function to alter how encoding is done before the CodecObj is invoked.

    For customized keys, you need to override the class and alter how keys are
    generated.

    """

    def __init__(self,packetizedFile,CodecObj=None,minIndex=0):
        self._index = minIndex
        #print "self._index={}".format(self._index)
        self.minIndex = minIndex
        self._packetizedFile = packetizedFile
        self._iterating = False
        if CodecObj == None:
            self._Codec = BaseCodec()
        else:
            self._Codec = CodecObj

    @property
    def index(self):
        #print "@property _index={}".format(self._index)
        return self._index
            
    @property
    def bytes_encoded(self):
        return self._packetizedFile.bytes_read

    # Derived classes should override only this method, unless other
    # functionality also needs to be altered
    # returns key, value
    def _encode(self):
        #print "return self._index={}".format(self._index)
        return self._index,self._packetizedFile.next()

    # Ideally, this would not be overriden, but if a subclass wants to
    # fully control all steps of encoding, it must be overriden
    def encode(self):
        return self._Codec.encode(self._encode())

    # Sub-classes that alter key names, may need to override this
    # Right now, only using one file descriptor, so this will
    # cause iteration to break. Either use this interface or iteration,
    # not both!
    def __getitem__(self,key):
        if self._iterating:
            assert False and "Warning: should not use __getitem__ while iterating!"
        return self._packetizedFile[key]

    # iterator
    def __iter__(self):
        self._index = self.minIndex
        self._iterating = True
        self._packetizedFile.__iter__()
        # reset to begnning of file
        return self

    def __next__(self):
        packet = self.encode()
        if packet:
            self._index += 1
            return packet
        else:
            self._iterating = False
            raise StopIteration()

    next = __next__ # support for python 2.7
        
class DecodePacketizedFile(object):
    def __init__(self,packetizedFile,CodecObj=None,minIndex=0):
        self._packetizedFile = packetizedFile
        self._data = {}
        self.minIndex = minIndex
        if CodecObj == None:
            self._Codec = BaseCodec()
        else:
            self._Codec = CodecObj

    def has_key(self,key):
        return self._data.has_key(key)

    #def __setitem__(self,key,value):
    #    self._data[key] = value
    #def __getitem__(self,key):
    #    return self._data[key]

    def writeToFile(self,key,value):
        #print (key, len(value), value, type(value))
        # make sure value is a string
        # try:
        #     if (type(value) is list) and (type(value[0]) is int or type(value[0]) is long):
        #         value = "".join([chr(x) for x in value])
            
        #     if (type(value) is list) and (type(value[0]) is str):
        #         value = "".join([x for x in value])
        #     # assert type(value) is str
        #     self._packetizedFile[key] = value
        #except Exception as e:
        if (type(value) is list) and (type(value[0]) is int or type(value[0]) is long):
            value = bytearray(value)
        elif (type(value) is list) and (type(value[0]) is str):
            value = bytearray([ord(x) for x in value])
        else:
            value = bytearray(value)
        self._packetizedFile[key] = value
                
    # Ideally, derived classes will only override this implementation
    def _decode(self,key,value):
        self.writeToFile(key,value)

    def decode(self,strand,bypass=False,input_key=None,input_value=None):
        #if input_key < self.minIndex or input_key >= self._packetizedFile.maxKey:
        #    return
        if bypass is False:
            key,val = self._Codec.decode(strand)
        else:
            key=input_key
            val=input_value
        #print key,val
        self._decode(key,val)


    @property
    def complete(self):
        return self._packetizedFile.complete

    def write(self):
        self._packetizedFile.write()

    def dummy_write(self):
        return self._packetizedFile.dummy_write()


if __name__ == "__main__":
    import os
    import sys
    from random import randint
    from dnastorage.util.packetizedfile import pf
    packetFile = pf.ReadPacketizedFile(sys.argv[1])
    out = pf.WritePacketizedFile("output.d",packetFile.size,120)
    assert out.complete==False
    i = 0
    for p in packetFile:
        out[i] = p
        i += 1
    assert out.complete==True
    out.write()
