from random import randint
import editdistance as ed

from dnastorage.codec.base_codec import BaseCodec
from dnastorage.codec import base_conversion
from dnastorage.util.stats import stats
import dnastorage.exceptions as err

#In this order, no codewords have a hamming distance
#smaller than 2, and only two consecutive codewords have a hamming
#distance of 2. This helps make sure that consecutive indices vary by
#at least two positions

cfc_all = ['CGTGAGCA', 'ACACGATG', 'ATAGCTGC', 'ATAGAGAT', 'CAGTCATC',
           'TCATGCAT', 'TCGACGCA', 'CTATGTCG', 'ACGAGTCA', 'CACATGCG',
           'ATGTGTGC', 'TGTGACAG', 'TCATAGCA', 'TCAGATAG', 'ACGATCTC',
           'CGTATATA', 'CAGACTGA', 'CTCAGCTA', 'CTACTACA', 'ATGAGAGT',
           'ACTGTGCT', 'CATATGAC', 'GTCAGAGT', 'GTACGCAC', 'CTAGACTC',
           'TACGTGTC', 'CATAGAGT', 'TCGAGTCG', 'ATCATGCG', 'CGAGTGTG',
           'GCTCTCTA', 'TCGTATAC', 'TGTCTGCA', 'GTGTATCG', 'TCACATAC',
           'ATCAGACA', 'CACAGCGA', 'CTACTGAC', 'ATCATGTC', 'GCGTGATA',
           'ACGAGTGC', 'TCGAGATG', 'CGTAGTCA', 'ATACATCT', 'GTGACGCG',
           'CAGTCACG', 'CTCGTGCG', 'GCAGAGTG', 'CTCGCGCA', 'TCACTGTC',
           'GCTCTGAG', 'ATCAGTCG', 'ATCGAGCA', 'TATCATCT', 'CTCACTGA',
           'CTGTATAT', 'ATCGACAG', 'ATCAGCGA', 'ACTGTATG', 'GTAGTGTG',
           'CGTGTGCT', 'CGCGACAG', 'GCAGATAT', 'CTCGTACA', 'GCTCGCAT',
           'CTGTGCAT', 'AGTAGCGC', 'CTACTGCG', 'CTCGCACT', 'TCTCGCTA',
           'CTATGCAC', 'GTCAGCTG', 'CTATCTGA', 'GCAGATCG', 'GAGACTCT',
           'ACTGCACT', 'GCGATGCT', 'TGAGTAGA', 'CTAGATAT', 'GTACGCTA',
           'TCATACAC', 'GTAGTATA', 'TACGTACT', 'CGTAGACG', 'GCAGACAG',
           'TGCGAGAT', 'ACATAGAT', 'GTCTACAG', 'GTCTCTGA', 'GTGAGTCA',
           'CATGACAT', 'CGACTCTG', 'CGATCACG', 'CAGCGTAC', 'GTCTCTCT',
           'ATGTATGA', 'ACGTCGTC', 'TCGAGTAC', 'TCGATCGC', 'GTCTAGCT',
           'GCATCTCT', 'ATAGCGCA', 'CATATGCT', 'CACAGTAT', 'AGTAGCTG',
           'TACGACAT', 'ACTAGAGC', 'GACACGCT', 'CTATGTGC', 'AGACTGCG',
           'TCGTGCAC', 'GTGAGTGC', 'CACACGTG', 'GAGAGATG', 'ATCGTGCT',
           'CACAGTCG', 'ACACACTC', 'GTGATCTC', 'GTCAGTAT', 'TATCACGT',
           'ACTGTAGC', 'CTCGACTA', 'CAGACGCT', 'CAGAGATC', 'TACTGATA',
           'CATGCTGC', 'ATGTGATA', 'CGACGATG', 'GCTCTGCT', 'TGCGAGCA',
           'ATCACAGA', 'GCTCGCTG', 'GCTCTGTC', 'GCATCTGA', 'AGACTACA',
           'CGATCATC', 'TAGTAGAG', 'GTCATGAT', 'CACATGTC', 'ATGTGACT',
           'CTGTAGTC', 'GTGTATGC', 'CATGCACG', 'GCAGCTCA', 'GCAGCACG',
           'GTCGATCT', 'CTGTGTAC', 'GAGACATC', 'CGTGATGC', 'CAGCTCTA',
           'ATGAGCTG', 'TAGTGACA', 'GCATGTGC', 'GCGTAGCT', 'TACGAGCT',
           'ACGACTGA', 'GACGACGC', 'CTACTCTC', 'CACACAGA', 'CACATGAT',
           'ACGTAGCA', 'ATAGCGTC', 'ACTGCTCA', 'TATCTGCT', 'CTGTACAC',
           'CGTAGAGA', 'CAGTGATA', 'CTGTACTA', 'ACATGACT', 'TACGTGCA',
           'TGTCTGAT', 'ACACGAGC', 'GCAGATGC', 'TAGTAGCA', 'TCACATCT',
           'ACGTACAG', 'TGTATAGA', 'TCAGTACA', 'GTCGCGTC', 'CGACGTCT',
           'CACTACGC', 'TAGTGTAC', 'GACGACAG', 'TATAGCAC', 'CTATGCTA',
           'CTGTGCGC', 'ACGTATGC', 'GACACTGA', 'CGATAGAT', 'CGACACTC',
           'CTAGATGC', 'CAGCGCAT', 'TCGTGTCA', 'CACACATC', 'TCATATCG',
           'ATCAGATG', 'CACATATA', 'TCTCGATG', 'GTCTGCTA', 'AGTGAGCG',
           'CGTGAGAT', 'CGCGATCT', 'CGAGTATA', 'GTAGTGCA', 'CTACGAGA',
           'CAGAGCTA', 'TCAGCGTC', 'GTGTAGCA', 'TCGACGTG', 'CATGTGTC',
           'TCTAGAGT', 'CAGATCTG', 'CGTATCTC', 'TGTATCTA', 'TCAGTGCG',
           'ACGATATA', 'TGAGAGCG', 'GCATGCTA', 'ATCGATAT', 'GAGAGCAC',
           'CACTACAG', 'AGACTAGC', 'CAGTGTCG', 'AGTAGTCT', 'GCAGACTC',
           'TACTCTGA', 'TATGTAGA', 'CTGTAGAG', 'CGATACGC', 'GCGAGAGT',
           'TCGCGCTG', 'TCACACGC', 'TAGTGATG', 'CACTAGCG', 'CATATCTA',
           'TCATAGTC', 'CTATCTCT', 'GTGAGACG', 'TGACGCTG', 'TCTCGAGC',
           'TCAGCACT', 'ATAGTACA', 'ACACAGCA', 'CTACTAGC', 'CTAGACAG',
           'CAGTGCAC', 'CTAGTGTC', 'GCGTGACT', 'TGTGAGCT', 'GCAGCGCT',
           'CAGTGACT', 'GCAGAGCA', 'GCGCGCTA', 'TCGACATC', 'TATAGCTA',
           'GTAGTACT']

cfc_extra = ['TACGATGA', 'CACTGCTG', 'CGTAGATC', 'TCACTCTA']
cfc_special_header = 'ATCGATGC'

cfc = cfc_all[0:256]

cfc_inv = {}
def create_cfc_inv():
    global cfc_inv
    for i,c in enumerate(cfc):
        cfc_inv[c] = i
    cfc_inv[None] = -1000
    
class CommaFreeCodewords(BaseCodec):
    def __init__(self,numberSymbols,CodecObj=None,Policy=None):
        BaseCodec.__init__(self,CodecObj,Policy=Policy)
        self._numberBytes = numberSymbols
        global create_cfc_inv
        create_cfc_inv()

    def _encode(self,strand):
        enc_strand = [ cfc[s] for s in strand ]
        return enc_strand

    def _decode_cfc(self, s):
        global cfc_inv
        if s in cfc_inv:
            return cfc_inv[s]
        return None

    def _decode(self,s):
        dec = self._decode_helper(s)
        return dec

    def exact_vote(self, s):
        exact = []        
        for i in range(len(s)):
            exact.append(self._decode_cfc(s[i:i+8]))
        return exact

    def inexact_vote(self, s):
        #print s
        global cfc
        D = {}
        for x in cfc:
            d = ed.eval(s,x)
            D[d] = D.get(d,[]) + [x]
        if 0 in D:
            assert len(D[0])==1
            return [100,D[0][0],0]        
        elif 1 in D:
            if len(D[1])==1:
                return [100,D[1][0],8]
            elif len(D[1])>1:
                p = D[1][ randint(0,len(D[1])-1) ]
                return [1.0/len(D[1])*100,p,8]
        elif 2 in D:
            if len(D[2])==1:
                return [100,D[2][0],8]
            elif len(D[2])>1:
                p = D[2][ randint(0,len(D[2])-1) ]
                return [1.0/len(D[2])*100,p,8]
        else:
                return [0.0,None,8]


    def get_next_exact(self, i, exact):
        while i < len(exact):
            if exact[i] != None:
                return i
            i += 1
        return i

    def _decode_helper(self, s):
        numSyms = self._numberBytes
        exact = self.exact_vote(s)

        #print len(s)
        
        new_strand = []
        i = 0
        while i < len(s) and len(new_strand) < numSyms:
            #print i
            if not(exact[i] is None):
                if i%8 != 0:
                    stats.inc("CFC8::misAlignedCodeword")
                new_strand.append(exact[i])
                #print i,"->",exact[i]
                i += 8
            else:
                e = err.DNABadCodeword("Missing expected CFC8 symbol")
                if self._Policy.allow(e):
                    j = self.get_next_exact(i,exact)
                    skipped = int(round((j-i)/8.0))
                    i_tmp = i
                    for k in range(skipped):
                        #new_strand.append(-1)
                        v = self.inexact_vote(s[i_tmp:i_tmp+8])
                        if v[0]==100:
                            new_strand.append(cfc_inv[v[1]])
                            i_tmp+=v[2]
                        else:
                            i_tmp+=8
                            new_strand.append(-1)
                    i = j
                else:
                    raise e
        
        if len(new_strand) != numSyms:
            e = err.DNAStrandPayloadWrongSize("Payload wrong size in CFC8")
            stats.inc("CFC8:payloadWrongSize")
            if self._Policy.allow(e):            
                if len(new_strand) < numSyms:
                    while len(new_strand) < numSyms:
                        new_strand.append(-1)
                else:
                    new_strand = new_strand[0:numSyms]
            else:
                raise e


        cnt=0
        found=False
        for _ in new_strand:
            stats.inc("CFC8::TotalCodewords")
            if _ == -1:
                found=True
                stats.inc("CFC8::UnidentifiedCodeword")
            else:
                stats.inc("CFC8::IdentifiedCodeword")
        if found:
            #print "found error"
            stats.inc("CFC8::strandsWithErrors")
            stats.inc("CFC8::strandsTotal")
        else:
            stats.inc("CFC8::strandsTotal")
            
        return new_strand
                
            # if j<3:
            #     i += (j-i)
            #     continue
            # v = inexact_vote(s[i:i+10])
            # if v[0] == 100:
            #     new_strand.append(v[1])
            #     if j < len(s):
            #         i += (j-i)
            #     else:
            #         i += v[0]
            # else:
            #     i += 1
    
