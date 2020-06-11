from dnastorage.codec.base_codec import BaseCodec
from dnastorage.codec.base_conversion import convertBase, convertFromBase

def binary_encode(packet):
    array = bytearray(packet)
    strand = []
    for b in array:
        strand.append( convertBase(2,b,8) )
    return "".join(strand)


def binary_decode(strand):
    i = 0
    array = []
    while i < len(strand):
        s = strand[i:i+8]
        byte = convertFromBase(2,s)
        array.append(chr(byte))
        i+=len(s)
    packet = bytearray(array)
    return packet


def binary_rotate_encode(strand,prev='A'):
    values = { 'A' : 0 ,
               'C' : 1 ,
               'G' : 2 ,
               'T' : 3   }
    complement = { 'A' : 'T' ,
               'C' : 'G' ,
               'G' : 'C' ,
               'T' : 'A'   }
    counts = { 'A':0, 'G':0,'C':0, 'T':0}
    l = [s for s in strand]
    for i in range(0,len(l)):
        swap = False
        if counts[l[i]] > counts[complement[l[i]]]:
            swap = True
        elif prev == l[i]:
            swap = True
        if swap:
            counts[complement[l[i]]]+=1
            l[i] = complement[l[i]]
            prev = l[i]
        else:
            counts[l[i]] += 1
            prev = l[i]
    #print("".join(l))
    return "".join(l)

def binary_unrotate_decode(strand):
    values = { 'A' : 0 ,
               'C' : 1 ,
               'G' : 2 ,
               'T' : 3   }
    complement = { 'A' : 'A' ,
               'C' : 'C' ,
               'G' : 'C' ,
               'T' : 'A'   }

    l = [s for s in strand]
    for i in range(0,len(l)):
        l[i] = complement[l[i]]
    return "".join(l)

class BinaryCodec(BaseCodec):
    def __init__(self,CodecObj=None,keyWidth=20):
        BaseCodec.__init__(self,CodecObj)
        self._keyWidth=keyWidth

    def _encode(self,packet):
        key = convertBase(2,packet[0],self._keyWidth)
        value = binary_encode(packet[1])
        return key+value

    def _decode(self,s):
        key = convertFromBase(2,s[0:self._keyWidth])
        value = binary_decode(s[self._keyWidth:])
        return key,value

class BinaryRotateCodec(BaseCodec):
    def __init__(self,CodecObj=None,prev='A'):
        BaseCodec.__init__(self,CodecObj)
        self._prev = prev

    def _encode(self,s):
        assert isinstance(s,str) or "Didn't get expected type."
        return binary_rotate_encode(s,self._prev)

    def _decode(self,s):
        assert isinstance(s,str) or "Didn't get expected type."
        return binary_unrotate_decode(s)




