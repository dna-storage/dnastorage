from io import BytesIO
import os
import sys
from random import randint
import unittest

from dnastorage.codec.base_conversion import *
from dnastorage.util.packetizedfile import *

class packetizedfile_py_test(unittest.TestCase):
    """ test packetizedfile support. """
    def test_packetizedfilestream(self):
        rbuff = BytesIO()
        tell = rbuff.tell()
        for i in range(10000):
            rbuff.write(bytearray(convertIntToBytes(i,2)))
        rbuff.seek(tell,0)            

        
        packetFile = ReadPacketizedFilestream(rbuff)
        packetFile.packetSize = 100

        wbuff = BytesIO()
        out = WritePacketizedFilestream(wbuff,10000*2,100)        
        assert out.complete==False
        
        i = 0
        for p in packetFile:
            out[i] = p
            i += 1

        assert out.complete==True
        out.write()

        y = 0
        while True:
            b = wbuff.read(2)
            if len(b)==0:
                break
            val = convertBytesToInt(b)
            assert val == y
            y = y+1


from dnastorage.util.stats import stats
class stats_py_test(unittest.TestCase):
    """ test stats. """
    def test_stats(self):
        ''' test arbitrary counter '''
        stats.inc("askdfjakjalk")
        stats.inc("askdfjakjalk")
        stats.inc("askdfjakjalk")
        assert stats["askdfjakjalk"]==3

    def test_unique(self):
        ''' test unique '''
        stats.inc("askdfjakjalk2")
        stats.unique("askdfjakjalk2",1)
        assert stats["askdfjakjalk2"]==1

            
if __name__ == "__main__":
    unittest.main()
