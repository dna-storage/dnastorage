import random
import editdistance as ed

import dnastorage.exceptions as err
from dnastorage.codec.base_codec import BaseCodec

class CombineCodewords(BaseCodec):
    def __init__(self,CodecObj=None,Policy=None):
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)

    def _encode(self, codeword_list):
        return "".join(codeword_list)

    def _decode(self, s):
        assert ("not used for decoding")


class NormalizeStrandLength(BaseCodec):
    def __init__(self,length,CodecObj=None,Policy=None):
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)
        self.length = length
        
    def _encode(self, phys_s):
        if len(phys_s) > self.length:
            e = err.DNAStrandPayloadWrongSize("NormalizeStrandLength: Strand is too long ({})".format(len(phys_s)))
            if self._Policy.allow(e):
                pass
            else:
                raise e
        elif len(phys_s) < self.length:
            add = self.length - len(phys_s)
            phys_s = phys_s + "".join([ random.choice('AGCT') for _ in range(add) ])
            
        return phys_s

    def _decode(self, s):
        assert ("not used for decoding")

        
class InsertMidSequence(BaseCodec):
    def __init__(self,seq,CodecObj=None,Policy=None):
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)
        self._seq = seq

    def _encode(self,strand):
        if strand.find(self._seq) != -1:
            e = err.DNAStrandPoorlyFormed("Found sequence already present while prepending {}"\
                                        .format(self._seq))
            if self._Policy.allow(e):
                pass
            else:
                raise e

        middle = int(len(strand)/2)
        return strand[0:middle]+self._seq+strand[middle:]

    def _decode(self,strand):
        index = strand.find(self._seq)
        if index != -1:
            return strand[0:index]+strand[index+len(self._seq):]
        else:
            er = err.DNAStrandMissingSequence("{} should have had {} inside it.".format(strand,self._seq))
            # there could be errors in the cut preventing us from seeing it
            # with an exact match, so now we look for an inexact match
            if self._Policy.allow(er):
                middle = int(len(strand)/2)
                slen = len(self._seq)
                res = []
                for m in range(middle-slen,middle):
                    sli = strand[m:m+slen]
                    res.append( ed.eval(sli,self._seq) )
                mn = min(res)
                if mn < slen/3:
                    place = middle - slen + res.index(mn)
                    return strand[0:place]+strand[place+slen:]
                else:
                    # just leave the strand along, and hopefully codewords can
                    # still be extracted properly
                    return strand
            else:
                raise er


class PrependSequence(BaseCodec):
    def __init__(self,seq,CodecObj=None,Policy=None,isPrimer=False):
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)
        self._seq = seq[:]
        self.is_primer = isPrimer

    def _encode(self,strand):
        if strand.find(self._seq) != -1:
            er = err.DNAStrandPoorlyFormed("Found sequence already present while prepending {}"\
                                        .format(self._seq))
            if self._Policy.allow(er):
                pass
            else:
                raise er
        return self._seq + strand

    def _decode(self,strand):
        index = strand.find(self._seq)
        if index != -1: # expected at beginning
            return strand[index+len(self._seq):]
        else:
            er = err.DNAStrandMissingSequence("{} should have had {} at beginning.".format(strand,self._seq))
            # there could be errors in the cut preventing us from seeing it
            # with an exact match, so now we look for an inexact match
            if self._Policy.allow(er):
                slen = len(self._seq)
                res = []
                # FIXME: how far in should we look?
                for m in range(0,50):
                    sli = strand[m:m+slen]
                    res.append( ed.eval(sli,self._seq) )
                mn = min(res)
                
                idx = res.index(mn)
                #print (res,mn)
                if mn < 5:
                    return strand[idx+slen:]
                else:
                    if self.is_primer:
                        raise err.DNAMissingPrimer("Missing primer {}".format(self._seq))
                    # just leave the strand along, and hopefully codewords can
                    # still be extracted properly
                    return strand
            else:
                raise er


class AppendSequence(BaseCodec):
    def __init__(self,seq,CodecObj=None,Policy=None,isPrimer=False):
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)
        self._seq = seq
        self.is_primer = isPrimer

    def _encode(self,strand):
        if strand.find(self._seq) != -1:
            er = err.DNAStrandPoorlyFormed("Found sequence already present while appending {}"\
                                        .format(self._seq))
            if self._Policy.allow(er):
                pass
            else:
                raise er
        
        return strand + self._seq

    def _decode(self,strand):
        index = strand.find(self._seq)
        slen = len(self._seq)
        if index != -1: # expected at end
            return strand[:index]
        else:
            er = err.DNAStrandMissingSequence("{} should have had {} at end.".format(strand,self._seq))
            # there could be errors in the cut preventing us from seeing it
            # with an exact match, so now we look for an inexact match
            if self._Policy.allow(er):
                slen = len(self._seq)
                res = []
                for m in range(len(strand)-2*slen,len(strand)):
                    sli = strand[m:m+slen]
                    res.append( ed.eval(sli,self._seq) )
                mn = min(res)
                idx = res.index(mn)
                if mn < 5:
                    return strand[:idx+len(strand)-2*slen]
                else:
                    if self.is_primer:
                        raise err.DNAMissingPrimer("Missing primer".format(self._seq))
                    # just leave the strand along, and hopefully codewords can
                    # still be extracted properly
                    return strand
            else:
                raise er


if __name__ == "__main__":
    import random
    
    cut = InsertMidSequence('AGATATAGGG')
    pre = PrependSequence('TAAAGGAAAAAG',CodecObj=cut)
    app = AppendSequence( 'CAAAATATAAAA',CodecObj=pre)
    app = app
    
    match = 0
    for _ in range(10000):
        strand = [ random.choice('AGCT') for _ in range(100) ]
        strand = "".join(strand)

        original = strand

        strand = app.encode(strand)

        copy = [ random.choice('AGCT') for _ in range(10) ] + [ _ for _ in strand ] + [ random.choice('AGCT') for _ in range(10) ]
        copy[ random.randint(0,9) ] = random.choice('AGCT')
        copy[ random.randint(140,150) ] = random.choice('AGCT')
        copy[ random.randint(75,78) ] = random.choice('AGCT')
        copy = "".join(copy)

        
        new =  app.decode(copy)
        if new == original:
            match += 1
        else:
            print ("--",len(original)-len(new))
            print (original)
            print (new)
            
    print (match / 10000.0 * 100)
