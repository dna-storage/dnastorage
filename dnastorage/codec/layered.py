from random import randint
from math import ceil,log

import dnastorage.exceptions as err
from dnastorage.codec.codecfile import EncodePacketizedFile, DecodePacketizedFile
from dnastorage.codec.phys import CombineCodewords
from dnastorage.codec.block import doMajorityVote, partitionStrandsIntoBlocks
from dnastorage.codec.block import reportBlockStatus

from dnastorage.util.stats import stats

import logging
logger = logging.getLogger('dna.storage.codec.LayeredCodec')

class LayeredEncoder(EncodePacketizedFile):

    def __init__(self,packetizedFile,minIndex=0,\
                 strandSizeInBytes=10,blockSizeInBytes=200*10,\
                 blockCodec=None,\
                 blockToStrandCodec=None,\
                 strandCodec=None,\
                 strandToCodewordCodec=None,\
                 codewordToPhysCodec=CombineCodewords(),\
                 physCodec=None,\
                 Policy=None):
        
        EncodePacketizedFile.__init__(self,packetizedFile,minIndex=minIndex)
        # set packetSize (can't hurt to do it again here)
        self._packetizedFile.packetSize = blockSizeInBytes

        # make sure sizes are consistent with each other
        # blocks must be whole multiple of strand size
        assert blockSizeInBytes % strandSizeInBytes == 0

        self.blockCodec = blockCodec
        self.blockToStrandCodec = blockToStrandCodec
        self.strandCodec = strandCodec
        self.strandToCodewordCodec = strandToCodewordCodec
        self.codewordToPhysCodec = codewordToPhysCodec
        self.physCodec = physCodec
                
        self.strandSizeInBytes = strandSizeInBytes
        self.blockSizeInBytes = blockSizeInBytes
        return

    def _layered_encode(self, block):
        # perform block encoding
        enc_block = self.blockCodec.encode(block)

        # convert block into strands
        strands = self.blockToStrandCodec.encode(enc_block)


        # protected logical strands
        final_strands = []
        tmp_strands = []
        for s in strands:
            ecc_s = self.strandCodec.encode(s)
            tmp_strands.append(ecc_s)
            cw_s = self.strandToCodewordCodec.encode(ecc_s)
            # phys codecs expect a DNA sequnce as a string
            phys_s = self.codewordToPhysCodec.encode(cw_s)            
            final = self.physCodec.encode(phys_s)
            final_strands.append(final)

        return final_strands

    
    def dummy_encode(self):
        """ Dummy encode let's us run the sequence of codecs and produce a single
            strand. We can look at that strand to determine its length or possibly
            other properties. Helpful when designing codecs.
        """
        block = [ randint(0,255) for _ in range(self.blockSizeInBytes) ]
        strands = self._layered_encode( (randint(0,255**2),block) )
        return strands
                
    def encode(self):        
        block = self._encode()
        try:
            block = (block[0],[ ord(_) for _ in block[1] ])
        except Exception as e:
            # This is pretty horrible, kids don't try this at home!
            block = (block[0],[ _ for _ in block[1] ])
        return self._layered_encode(block) # get entire block
        

class LayeredDecoder(DecodePacketizedFile):
    
    def __init__(self,packetizedFile,minIndex=0,\
                 strandSizeInBytes=10,blockSizeInBytes=200*10,\
                 blockIndexSize=2,
                 physCodec=None,\
                 physToStrandCodec=None,\
                 strandCodec=None,\
                 strandToBlockCodec=None,\
                 blockCodec=None,\
                 Policy=None,
                 intraBlockIndexSize=1):\
                         
        DecodePacketizedFile.__init__(self,packetizedFile,minIndex=minIndex)
        # set packetSize (can't hurt to do it again here)
        self._packetizedFile.packetSize = blockSizeInBytes

        # make sure sizes are consistent with each other
        # blocks must be whole multiple of strand size
        assert blockSizeInBytes % strandSizeInBytes == 0

        self.physCodec = physCodec
        self.physToStrandCodec = physToStrandCodec
        self.strandCodec = strandCodec
        self.strandToBlockCodec = strandToBlockCodec        
        self.blockCodec = blockCodec
                
        self.strandSizeInBytes = strandSizeInBytes
        self.blockSizeInBytes = blockSizeInBytes
        self.all_strands = []
        self.blockIndexSize = blockIndexSize
        self.intraBlockIndexSize = intraBlockIndexSize
        self._Policy = Policy

        self.strand_errors = 0
        self.block_errors = 0
        
        logger.info("strandSizeInBytes = {}".format(strandSizeInBytes))
        logger.info("blockSizeInBytes = {}".format(blockSizeInBytes))
        logger.info("blockIndexSize = {}".format(blockIndexSize))
        logger.info("intraBlockIndexSize = {}".format(intraBlockIndexSize))
        return


    def _layered_decode_phys_to_strand(self, phys_strand):
        # perform block encoding
        try:
            phys_s = self.physCodec.decode(phys_strand)
            cw_s = self.physToStrandCodec.decode(phys_s)
            s = self.strandCodec.decode(cw_s)
            stats.inc("LayeredDecoder::phys_to_strand::succeeded")
        except err.DNAStorageError as p:
            stats.inc("LayeredDecoder::phys_to_strand::failed")
            s = [-1] + [ 0 for _ in range(self.strandSizeInBytes-1) ]
            self.strand_errors += 1
                
        return s

    def decode_from_phys_to_strand(self, s):
        return self._layered_decode_phys_to_strand(s)
    
    def decode(self, phys_strand, bypass=False, input_key=None, input_value=None):
        if bypass==False:
            try:            
                s = self._layered_decode_phys_to_strand(phys_strand)
                self.all_strands.append(s)
            except err.DNAMissingPrimer as p:
                # just ignore strands that don't have a required primer
                pass
            except err.DNAStorageError as e:
                if self._Policy.allow(e):
                    pass
                else:
                    raise e
        else:
            self.all_strands.append(input_value)

    def _attempt_final_decoding(self):
        # do voting here!!        
        self.voted_strands = doMajorityVote(self.all_strands,\
                                            self.blockIndexSize+self.intraBlockIndexSize)
        if self.blockIndexSize>0:
            blocks = partitionStrandsIntoBlocks(self.all_strands,self.blockIndexSize)
        else:
            blocks = [ (0,self.all_strands) ]
        #print "number of blocks=",len(blocks)

        reportBlockStatus(blocks,self.minIndex,
                          self.blockIndexSize,self.intraBlockIndexSize)
        
        #blocks.sort()
        for b in blocks:
            idx = b[0]
            if idx < self.minIndex or idx >= self._packetizedFile.maxKey:
                # this happens due to errors in strands, and we should just
                # discard these erroneous blocks
                # If we want to reclaim some of these strands, it needs
                # to happen as part of the error processing per strand
                stats.inc("LayeredCodec::_attempt_final_decoding::indexOutOfRange")
                continue

            stats.inc("LayeredCodec::_attempt_final_decoding::numberOfBlocks")
            try:
                #print "attempt",idx,len(b[1])
                b_contig = self.strandToBlockCodec.decode(b)
                #print "attempt",b_contig[0],len(b_contig[1])
                b_noecc = self.blockCodec.decode(b_contig)
                #print "attempt",b_noecc[0],len(b_noecc[1])
                self.writeToFile(b_noecc[0],b_noecc[1])

            except err.DNAStorageError as e:
                self.block_errors += 1
                if self._Policy.allow(e):
                    continue                
                else:
                    raise e
            except Exception as e:
                print("LayeredCodec._attempt_final_decoding caught error: "+str(e))
                self.block_errors += 1
                print (b)
                
            
    def write(self):
        self._attempt_final_decoding()
        super(LayeredDecoder,self).write()

    def only_write(self):
        super(LayeredDecoder,self).write()

    def dummy_write(self):
        self._attempt_final_decoding()
        return self._packetizedFile.dummy_write()
        
if __name__ == "__main__":
    from dnastorage.codec.rscodec import ReedSolomonOuterCodec
    from dnastorage.codec.strand import ReedSolomonInnerCodec
    from dnastorage.exceptions import NoTolerance
    from dnastorage.codec.commafreecodec import *
    from dnastorage.codec.phys import *
    
    pol = NoTolerance()

    payload=9
    
    outerECC = 70 # strands
    blockSize = payload*185 # 185 normal strands
    blockCodec = ReedSolomonOuterCodec(packetSize=blockSize,errorSymbols=70,payloadSize=payload,Policy=pol)

    blockToStrand = BlockToStrand(9,(185+70)*payload,Policy=pol)
    
    # take index into account
    # better to just say number of error symbols
    strandCodec = ReedSolomonInnerCodec(2,Policy=pol)

    codewords = CommaFreeCodewords(14)

    cut = InsertMidSequence('AGGTACCA')
    pre = PrependSequence('CAGGTACGCAGTTAGCACTC',CodecObj=cut, isPrimer=True)
    app = AppendSequence( 'CGTGGCAATATGACTACGGA',CodecObj=pre, isPrimer=True)

    flank = PrependSequence('CAGGAGAATGCCTTCCTAGG',CodecObj=app, isPrimer=True)
    flank = AppendSequence('CCTCGGTTCTTCTTGACCAG',CodecObj=flank, isPrimer=True)
    
    physCodec = flank
    pf = ReadPacketizedFile(sys.argv[1])
    
    enc = LayeredEncoder(pf,blockSizeInBytes=blockSize,strandSizeInBytes=payload,\
                         blockCodec=blockCodec,\
                         blockToStrandCodec=blockToStrand,\
                         strandCodec=strandCodec,\
                         strandToCodewordCodec=codewords,\
                         codewordToPhysCodec=CombineCodewords(),\
                         physCodec=physCodec)

    s = enc.dummy_encode()
    print (s[10])
    print (len(s),len(s[0]))

    all_strands = [ _ for _ in enc ]
    
    wpf = WritePacketizedFile(sys.argv[2],pf.size,blockSize)

    dec = LayeredDecoder(wpf,blockSizeInBytes=blockSize,strandSizeInBytes=payload,
                         blockCodec=blockCodec,\
                         strandCodec=strandCodec,\
                         physCodec=physCodec,\
                         physToStrandCodec=codewords,\
                         strandToBlockCodec=blockToStrand,\
                         Policy=pol)

    
    for ss in all_strands:
        for s in ss:
            dec.decode(s)

    dec._attempt_final_decoding()            
    dec.write()
    assert dec.complete

    
