#!/usr/bin/python
import os
import sys

import logging
logger = logging.getLogger("dna.storage.util.packetizedfile")
logger.addHandler(logging.NullHandler())

class WritePacketizedFilestream:
    """
    Write a file by receiving (index,value) packets of data of a specific size at a given index. Value
    is of packetSize length.  index specifies where in the file to write the data.  Packets can be
received in any order. But, the file cannot be written until the file is complete, in other words,
    it has received all size/packetSize packets.  The indices are assumed to be densely numbered from 0
    to size/packetSize by 1.  Pythonically: range(0,size/packetSize+(size-size/packetSize*packetSize),1).
    """
    def __init__(self,fd,size,packetSize,minKey=0,zeroFillMissing=False):
        self.__fd = fd
        self.size = size
        self.__set_packetSize(packetSize)
        self.__data = {}
        self.minKey = int(minKey)                      # inclusive
        self.zeroFillMissing = zeroFillMissing
        
    def has_key(self,key):
        return key in self.__data
    
    def __setitem__(self,key,value):        
        if (key >= self.minKey) and (key < self.maxKey):
            self.__data[key] = value
        else:
            print ("not in range:",key,value,self.minKey,self.maxKey)
    def __getitem__(self,key):
        assert key >= self.minKey and key < self.maxKey
        return self.__data[key]

    @property
    def maxKey(self):
        return int(self.minKey+self.numberOfPackets) # exclusive
    
    @property
    def numberOfPackets(self):
        if (self.size % self.packetSize) > 0:
            return int(self.size / self.packetSize + 1)
        else:
            return int(self.size / self.packetSize)

    # packet property
    def __set_packetSize(self,val):
        self.__packetSize = val
    def __get_packetSize(self):
        return self.__packetSize
    packetSize = property(__get_packetSize,__set_packetSize)

    @property
    def lastPacketSize(self):
        # this function should never return 0, not the same as modulo
        # should return a number >= 1 and <= self.packetSize
        return self.size - (self.numberOfPackets-1)*self.packetSize

    @property
    def complete(self):
        if len(self.__data.keys()) < self.numberOfPackets:
            #print "too few keys {}".format(self.numberOfPackets)
            return False
        for i in range(self.minKey,self.maxKey):
            if not (i in self.__data.keys()):
                #print "missing key {}".format(i)
                return False
        return True

    def getMissingKeys(self):
        keys = self.__data.keys()
        missing = []
        for i in range(self.minKey,self.maxKey):
            if not i in keys:
                missing.append(i)
        return missing

    def hasMissingKeys(self):
        return len(self.getMissingKeys()) > 0

    ## Warning: requires buffering the whole file!
    def write(self):
        #emptyString = '\x00'*self.packetSize
        emptyString = bytearray(self.packetSize)
        for i in range(self.minKey,self.maxKey):
            #print (self.__data[i])
            if i in self.__data:
                if i == self.maxKey-1:
                    self.__fd.write(self.__data[i][0:self.lastPacketSize])
                else:
                    self.__fd.write(self.__data[i])
            else:
                if self.zeroFillMissing:
                    self.__fd.write(emptyString)
                
        # for key,value in items:
        #     if i < key:
        #         while i < key:
        #             self.__fd.write(emptyString)
        #             i+=1
        #     if i == self.numberOfPackets-1:
        #         self.__fd.write(value[0:self.lastPacketSize])
        #     else:
        #         self.__fd.write(value)
        #     i+=1
        self.__fd.flush()
            
    def close(self):
        if self.__fd != sys.stdout and self.__fd != sys.stderr:
            self.__fd.close()

    def dummy_write(self):
        items = self.__data.items()
        items.sort(cmp=lambda x,y: cmp(x[0],y[0]))
        i = 0
        emptyString = '\x00'*self.packetSize
        output_data=""
        for key,value in items:
            if value is not type(str):
                value2="".join(value)
            else:
                value2 = value
            if i < key:
                while i < key and i < self.numberOfPackets-1:
                    output_data+=emptyString
                    i+=1
            if i == self.numberOfPackets-1:
                output_data+=value2[0:self.lastPacketSize]
                i+=1
                break
            else:
                output_data+=value2
                i+=1
            if i==self.numberOfPackets:
                break
        while(i<self.numberOfPackets):
            if i == self.numberOfPackets-1:
                output_data+=emptyString[0:self.lastPacketSize]
            else:
                output_data+=emptyString
            i+=1
        #print i,self.numberOfPackets,len(output_data),self.size
        assert i == self.numberOfPackets and len(output_data)==self.size 
        return output_data



class WritePacketizedFile(WritePacketizedFilestream):
    def __init__(self,filename,size,packetSize):
        WritePacketizedFilestream.__init__(self,open(filename,"wb"),size,packetSize)
        self.__filename = filename



class ReadPacketizedFilestream:
    """
    Read a file by breaking it up into packetSize pieces. If the last piece is smaller
    than packetSize, it's padded with zeros. Clients are guaranteed all packets are uniform
    size.
    
    This is an iterable object. Packets can be read using the iterator or by requesting
    a specific index between 0 and size/packetSize.
    
    This class expects a file descriptor. There is a derived class that accepts a filename.
    """
    def __init__(self,fd):
        self.__fd = fd
        self.__set_packetSize(120)
        self.__read_size = 0
        self._RS=False

    def read(self):
        b = bytes(self.__fd.read(self.packetSize))
        self.__read_size += len(b)
        #only pad out if it is not RS, RS encodes by blocks so padding may lead to many useless strands
        if b and len(b) != self.packetSize and not self._RS:
            b = b.ljust(self.packetSize,bytes(1))

        return b

    # packet property
    def __set_packetSize(self,val):
        self.__packetSize = val
    def __get_packetSize(self):
        return self.__packetSize
    packetSize = property(__get_packetSize,__set_packetSize)

    # iterator
    def __iter__(self):
        # reset to begnning of file
        self.__fd.seek(0)
        return self

    def __next__(self):
        b = self.read()
        if b:
            return b
        else:
            raise StopIteration()

    next = __next__  # preserve python 2.7 functionality
        
    def __getitem__(self,key):
        self.__fd.seek(key*self.packetSize)
        return self.read()

    @property
    def size(self):
        return os.fstat(self.__fd.fileno()).st_size
    @property
    def bytes_read(self):
        return self.__read_size

    @property
    def numberOfPackets(self):
        if (self.size % self.packetSize)==0:
            return self.size / self.packetSize
        else:
            return self.size / self.packetSize + 1


class ReadPacketizedFile(ReadPacketizedFilestream):
    def __init__(self,filename):
        ReadPacketizedFilestream.__init__(self,open(filename,"rb"))
        self.__filename = filename
    @property
    def filename(self):
        return self.__filename


if __name__ == "__main__":
    import os
    import sys
    from random import randint
    packetFile = ReadPacketizedFile(sys.argv[1])

    out = WritePacketizedFile("output.d",packetFile.size,120)

    assert out.complete==False

    i = 0
    for p in packetFile:
        out[i] = p
        i += 1

    assert out.complete==True
    out.write()
