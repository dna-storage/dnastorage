#!/usr/bin/python

from dnastorage.codec.base_codec import BaseCodec
from dnastorage.codec import base_conversion

class HuffmanTableBase:
    
    class Node:
        def __init__(self, nbase, symbol, weight=None, childlist=[]):
            self._nbase = nbase
            self.symbol = symbol
            self.weight = weight
            self.enc = None
            self._childlist = childlist
            if weight==None and len(childlist)==0:
                self.weight = -1.0
            elif weight==None:
                self.weight = sum([ x.weight for x in childlist]) 
               
        def __lt__(self, other):
            return self.weight < other.weight
        def __le__(self, other):
            return self.weight <= other.weight
        def __gt__(self, other):
            return self.weight > other.weight
        def __ge__(self, other):
            return self.weight >= other.weight
        def __eq__(self,other):
            return self.weight == other.weight
        def __ne__(self,other):
            return self.weight != other.weight
                        
        def dec(self,val,symLut):
            #print val[0],symLut,self.enc,self.symbol,symLut,self._childlist
            if self.symbol != None:
                return len(self.enc),self.symbol
            if len(val) == 0:
                return 0,None
            if len(self._childlist) <= symLut[val[0]]:
                return 0,None            
            # if len(val)==1:
            #     newval = ""
            # else:
            #     newval = val[1:]
            #print val
            return self._childlist[ symLut[val[0]] ].dec( val[1:], symLut )


        def _insert(self, enc, val, nbase, symLut):
            # figure if leaf or not
            #print enc,val
            if len(enc)==1:
                n = HuffmanTableBase.Node(nbase,val,1.0,[])
                self._childlist.append( n )
                return n
            else:
                if symLut[ enc[0] ] >= len(self._childlist):
                    self._childlist.append( HuffmanTableBase.Node(nbase,None,None,[]) )
                return self._childlist[ symLut[enc[0]] ]._insert(enc[1:],val,nbase,symLut)
                    
                    
        def symbol(self):
            return self.symbol

        def __str__(self):
            s = ""
            if len(self._childlist):
                for c in self._childlist:
                    s += "->"+str(c)

            #return "{}:({},{})".format(self.enc,self.symbol,self.weight)
            return "{}".format(self.symbol) + s

        def contains(self,sym):
            if self.symbol == sym:
                return True
            else:
                for c in self._childlist:
                    if c.contains(sym):
                        return True
                return False

        def assign_enc(self, bsyms, enc):
            self.enc = enc
            #print self._childlist
            if len(self._childlist) > 0:
                self._enc = enc
                for s,n in zip(bsyms,self._childlist):
                    n.assign_enc(bsyms,enc+s)
            #else:
            #    print str(self)

        def prevent_all_ones(self, nsyms):
            if len(self._childlist) == nsyms:
                if len(self._childlist[nsyms-1]._childlist)==0:
                    new = HuffmanTableBase.Node(self._nbase,None,None,[self._childlist[nsyms-1]])
                    self._childlist[nsyms-1] = new
                else:
                    self._childlist[nsyms-1].prevent_all_ones(nsyms)
                
            
            
    def __init__(self, nbase, base_syms, symbols, weights=None, prevent_ones=False):
        """ Build a huffman encoder/decoder table under the following assumptions:    """
        """     nbase: number of symbols in codeword (4 for DNA, 2 for binary)        """
        """     base_syms: symbols used in codeword (['A','C','G','T'] or ['0', '1']) """
        """     symbols: the symbols we'll replace with huffman codewords             """
        """     weights: frequency of the symbols (same order as symbols)             """
        self._nbase = nbase
        self._base_syms = base_syms
        self._base_syms_lookup = { b : i for i,b in enumerate(base_syms) }
        self._symbols = symbols
        self._weights = weights
        self._prevent_ones = prevent_ones
        self._nodes = []
        self._enc_table = None
        assert nbase == len(base_syms)
        self.memoize = {}
        if symbols != None:
            # we may be initialized with no symbols, in that case, don't do this
            # work to setup the huffman table, just wait
            if weights == None:
                self._weights = [ 1.0 / len(symbols) for _ in range(len(symbols)) ]
            else:
                W = sum(weights)
                #print W
                self._weights = [ float(x) / W for x in weights ]
            #print self._weights
            assert len(self._symbols) == len(self._weights)
            for s,w in zip(self._symbols,self._weights):
                n = HuffmanTableBase.Node(nbase,s,weight=w)
                self._nodes.append( n )
            self._nodes.sort(key=lambda n: n.weight)
        else:
            self.root = HuffmanTableBase.Node(self._nbase,None,None,[])


    def _insert(self, enc, val):
        n = self.root._insert(enc,val,self._nbase,self._base_syms_lookup)
        self._nodes.append(n)
        
    @classmethod
    def from_raw_table_hack(cls, table, nbase, base_syms, prevent_ones=False):
        vals = [ t[1] for t in table ]
        weights = [ 1.0 / nbase**t[0] for t in table ]
        return cls(nbase, base_syms, vals, weights, prevent_ones)

    @classmethod
    def from_raw_table(cls, table, nbase, base_syms, prevent_ones=False):
        ht = cls(nbase,base_syms,None,prevent_ones=prevent_ones)
        i = 0
        prev_l = table[0][0]
        for l,val in table:
            if l!=prev_l:
                i = (i) << (l - prev_l)
                prev_l = l
            enc = "{:0{w}b}".format(i,w=l)
            #print i,enc,val
            ht._insert(enc,val)
            i += 1
        #print ht.root
        ht.root.assign_enc(base_syms,"")
        #print ht.root
        return ht


    def get_raw_table(self,length_only=False):
        table = []
        for n in self._nodes:
            if n.symbol != None:
                if length_only:
                    table.append( [len(n.enc), n.symbol] )
                else:
                    table.append( [n.enc, n.symbol] )
                                        
        table.sort(key=lambda x: x[0])
        return table

    def get_tables(self):
        table = []
        for n in self._nodes:
            if n.symbol != None:
                table.append( [n.enc, n.symbol] )

        table.sort(key=lambda x: x[1] )

        enc = { x[1] : x[0] for x in table }
        dec = { x[0] : x[1] for x in table } 

        return (enc, dec)

    def average_length(self):
        L = 0.0
        for n in self._nodes:
            if n.symbol != None:
                L += n.weight * len(n.enc)
        return L

    def histogram(self):
        h = {}
        #H = {}
        for n in self._nodes:
            if n.symbol != None:
                l = len(n.enc)
                h[l] = h.get(l,0) + 1
                #H[l] = H.get(l,[]) + [ n.enc ]
                #H[l].sort()
        return h

    def _build_tree(self):
        _queue = []
        for n in self._nodes:
            _queue.append( n )
        _queue.sort(key=lambda n: n.weight)

        if len(_queue)==1:
            nodes = []
            nodes.append( _queue.pop(0) )
            new = HuffmanTableBase.Node(self._nbase,None,None,nodes)
            _queue.append(new)
        else:
            while len(_queue) > 1:
                nodes = []
                i = 0
                while len(_queue) > 0 and i < self._nbase:
                    nodes.append( _queue.pop(0) )
                    i += 1
                new = HuffmanTableBase.Node(self._nbase,None,None,nodes)
                self._nodes.append(new)
                _queue.append(new)
                # Note, this is ineffecient.  We should replace this dumb sort
                # with something more efficient like a min-heap
                _queue.sort(key=lambda n: n.weight)

        assert len(_queue) == 1
        r = _queue[0]

        if self._prevent_ones:
            r.prevent_all_ones(self._nbase)
        
        r.assign_enc(self._base_syms,"")
        self.root = _queue[0]

        
    def encode(self, val):
        if self._enc_table == None:
            table = []
            for n in self._nodes:
                if n.symbol != None:
                    table.append( [n.enc, n.symbol] )

            table.sort(key=lambda x: x[1])
            self._enc_table = { x[1] : x[0] for x in table }
        return self._enc_table[val]

    def decode(self, val):
        if self.memoize.has_key(val):
            return self.memoize[val]        
        ans = self.root.dec(val,self._base_syms_lookup)
        self.memoize[val] = ans
        return ans

class HuffmanTable(HuffmanTableBase):                
    def __init__(self, nbase, base_syms, symbols, weights=None,prevent_ones=False):
        HuffmanTableBase.__init__(self,nbase,base_syms,symbols,weights,prevent_ones)
        if symbols != None:
            self._build_tree()


class ErrorWhileDecodingTable(HuffmanTableBase):
    def __init__(self, nbase, base_syms, symbols=None, weights=None, prevent_ones=False):
        HuffmanTableBase.__init__(self,nbase,base_syms,symbols,weights,prevent_ones)

    def decode(self, val):
        return 0,0

    def encode(self, val):
        assert False and "Never use this for encoding"

    def average_length(self):
        return -1

    def histogram(self):
        return {}
            
class LengthLimitedHuffmanTable(HuffmanTableBase):
    # based on package-merge algorithm
    def __init__(self, L, nbase, base_syms, symbols, weights=None,prevent_ones=False):
        assert nbase==2
        assert nbase**L >= len(symbols)
        HuffmanTableBase.__init__(self,nbase,base_syms,symbols,weights=weights,prevent_ones=prevent_ones)
        self._original_nodes = self._nodes[:]
        merge = []
        new_nodes = []
        for _ in range(L,0,-1):            
            merge = new_nodes + self._nodes[:]
            #print "L = {}".format(_)
            merge.sort(key = lambda n: n.weight)
            even = len(merge)//2*2
            # drop last packet if length of merge is odd
            package = [ [merge[i],merge[i+1]] for i in range(0,even,2) ] 
            new_nodes = []
            for p in package:
                new = HuffmanTable.Node(self._nbase,None,weight=None,childlist=p)
                new_nodes.append(new)        
        merge.sort(key = lambda n: n.weight)
        code_length = []
        # adjust the weights based on code lengths
        for n in self._nodes:
            n.weight = 1/(nbase**sum([ m.contains(n.symbol) for m in merge ]))
            #print n
        # build Huffman code
        self._build_tree()
        

if __name__ == "__main__":
    import random
    R = random.Random()
    l = [ HuffmanTable.Node(2,R.randint(0,256),R.random()) for _ in range(256) ] 
    l.sort()
    
    syms = [ x for x in range(256) ]
    w = [ 1.0/256 for x in range(256) ]
    for i in range(ord('A'),ord('z'),1):
        w[i] = w[i]*2.0

    ht5 = LengthLimitedHuffmanTable(9, 2, ['0','1'], syms, w, True)
    print (ht5.average_length())
    print (ht5.histogram())


    rt = ht5.get_raw_table(True)
    print (rt)
    ht6 = HuffmanTable.from_raw_table(rt,2,['0','1'])

    print( ht6.get_raw_table())

    for x,y in zip(ht5.get_raw_table(),ht6.get_raw_table()):
        if x!=y:
            print (x,y)
    
    if ht5.get_raw_table(True) == ht6.get_raw_table(True):
        print ("They match ")

    
