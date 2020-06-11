from dnastorage.codec.phys import CombineCodewords, PrependSequence, AppendSequence, InsertMidSequence
from dnastorage.codec.commafreecodec import CommaFreeCodewords
from dnastorage.codec.strand import ReedSolomonInnerCodec
from dnastorage.codec.block import ReedSolomonOuterCodec, BlockToStrand
from dnastorage.codec.layered import LayeredEncoder, LayeredDecoder

def reverse_complement(seq):
    if len(seq)==0:
        return seq
    complement = {'T':'A', 'G':'C', 'C':'G', 'A':'T'}
    r = [complement[x] for x in seq]
    r.reverse()
    return "".join(r)

def customize_RS_CFC8(is_enc,pf,primer5,primer3,intraBlockIndex=1,\
                      interBlockIndex=2,innerECC=2,strandSizeInBytes=15,\
                      blockSizeInBytes=15*185,Policy=None,\
                      withCut=None,outerECCStrands=None,minIndex=0):
    assert blockSizeInBytes % strandSizeInBytes == 0
    payload=strandSizeInBytes
    blockStrands = blockSizeInBytes / strandSizeInBytes
    if outerECCStrands is None:
        outerECCStrands = 255-blockStrands # strands
    else:
        assert outerECCStrands + blockStrands <= 255
    assert blockStrands + outerECCStrands < 256

    index = intraBlockIndex + interBlockIndex

    blockCodec = ReedSolomonOuterCodec(packetSize=blockSizeInBytes,\
                                       errorSymbols=outerECCStrands,payloadSize=payload,Policy=Policy)

    blockToStrand = BlockToStrand(payload,(blockStrands+outerECCStrands)*payload,Policy=Policy,\
                                  intraIndexSize=intraBlockIndex,\
                                  interIndexSize=interBlockIndex,filterZeroes=True)
    
    # take index into account
    # better to just say number of error symbols
    strandCodec = ReedSolomonInnerCodec(innerECC,Policy=Policy)

    codewords = CommaFreeCodewords(payload+innerECC+index,Policy=Policy)

    if not (withCut is None):
        cut = InsertMidSequence(withCut)
    else:
        cut = None
    
    pre = PrependSequence(primer5,isPrimer=True,Policy=Policy,CodecObj=cut)
    app = AppendSequence(reverse_complement(primer3), CodecObj=pre, isPrimer=True, \
                         Policy=Policy)    
    physCodec = app

    if is_enc==True:
        enc = LayeredEncoder(pf,blockSizeInBytes=blockSizeInBytes,\
                             strandSizeInBytes=strandSizeInBytes,\
                         blockCodec=blockCodec,\
                         blockToStrandCodec=blockToStrand,\
                         strandCodec=strandCodec,\
                         strandToCodewordCodec=codewords,\
                         codewordToPhysCodec=CombineCodewords(),\
                             physCodec=physCodec,Policy=Policy,minIndex=minIndex)
        return enc
    else:
        dec = LayeredDecoder(pf,blockSizeInBytes=blockSizeInBytes,\
                             strandSizeInBytes=strandSizeInBytes,\
                             blockIndexSize=interBlockIndex,\
                             blockCodec=blockCodec,\
                             strandCodec=strandCodec,\
                             physCodec=physCodec,\
                             physToStrandCodec=codewords,\
                             strandToBlockCodec=blockToStrand,Policy=Policy,minIndex=minIndex)
        return dec

