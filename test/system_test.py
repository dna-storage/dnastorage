from io import BytesIO
import os
import sys
from random import randint
import unittest


from dnastorage.system.header import *
class header_py_test(unittest.TestCase):
    """ test stats. """
    def test_header(self):
        strands = encode_file_header("",0xA,2,[1,2,3,4],"A"*19+"G","T"*19+"G")
        #for s in strands:
        #    print ("{}: strand={}".format(len(s), s))    
        header = decode_file_header(strands,"A"*19+"G","T"*19+"G")
        print (header)
        assert header['filename'] == "".join([ chr(0) ])
        assert header['size'] == 2
        assert header['formatid'] == 0xA
        assert header['other_data'] == [1,2,3,4]
        
from dnastorage.system.formats import *
class formats_py_test(unittest.TestCase):
    """ test stats. """
    def test_formats(self):
        f = file_system_formats()
        assert len(f) > 0

    def test_format_by_abbrev(self):
        assert file_system_encoder_by_abbrev("RS+CFC8") == ENC_RS_CFC8_200
        assert file_system_format_packetsize(32) == 16

from dnastorage.system.dnafile import *
class dnafile_py_test(unittest.TestCase):
    """ test stats. """
    def test_dnafile(self):
        wf = DNAFile.open("out.dna","w",primer3='TTTG',primer5='AAAG',format_name='RS+CFC8')
        for i in range(1000):
            wf.write( bytearray(convertIntToBytes(i,4)) )
        wf.close()
        rf = DNAFile.open("out.dna","r",primer3='TTTG',primer5='AAAG')
        i = 0
        while True:
            b = rf.read(4)
            if len(b)==0:
                break
            val = convertBytesToInt(b)
            assert val == i
            i = i+1
        rf.close()

from dnastorage.codec.base_conversion import convertIntToBytes, convertBytesToInt      
class segmentedfile_py_test(unittest.TestCase):
    """ test stats. """
    def test_dnafile(self):
        wf = SegmentedWriteDNAFile(primer3='T'*19+'G',primer5='A'*19+'G',format_name='RS+CFC8+RE1',output="out.dna",fsmd_abbrev='FSMD-1')    
        #for i in range(1000):
        #    wf.write( bytes([x for x in convertIntToBytes(i,4)]) )
        wf.new_segment('RS+CFC8+RE2','AT'+'A'*17+'G','TA'+'T'*17+'G')        
        for i in range(30):
            wf.write( bytes([x for x in convertIntToBytes(i,4)]) )
        wf.close()
        rf = SegmentedReadDNAFile(primer3='T'*19+'G',primer5='A'*19+'G',input="out.dna",fsmd_abbrev='FSMD-1')
        print("Should print out 0 to 30: ")
        out = []
        while True:
            s = rf.read(4)
            if len(s)==0:
                break
            n = convertBytesToInt([x for x in s])
            out.append(n)
            print(n)            
        print("Done.")
        rf.close()
        assert out == [_ for _ in range(30)]
        
if __name__ == "__main__":
    unittest.main()
