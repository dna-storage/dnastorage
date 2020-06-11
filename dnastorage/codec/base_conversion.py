bases = ['A', 'C', 'G', 'T']

values = { 'A' : 0 ,
           'C' : 1 ,
           'G' : 2 ,
           'T' : 3   }

def randomTernary(length):
    return "".join([ bases[x%3] for x in range(length) ])

def convertTernaryHelper(dec,s):
    m = int(dec % 3)
    q = int(dec / 3)
    s = s + bases[m]
    if q > 0:
        return convertTernaryHelper(q,s)
    else:
        return s

def convertTernary(dec,length):
    s = convertTernaryHelper(dec,'')
    s = s.ljust(length,bases[0])
    return s

def convertQuarnaryHelper(dec,s):
    m = int(dec % 4)
    q = int(dec / 4)
    s = s + bases[m]
    if q > 0:
        return convertQuarnaryHelper(q,s)
    else:
        return s

def convertQuarnary(dec,length):
    s = convertQuarnaryHelper(dec,'')
    s = s.ljust(length,bases[0])
    return s

def convertBaseHelper(base : int, dec : int, s : str, symbols = bases):
    m = int(dec % base)
    q = int(dec / base)
    #print(s)
    s = s + symbols[m]
    #print(s)
    if q > 0:
        return convertBaseHelper(base,q,s,symbols)
    else:
        return s
    
def convertBytetoBinary(x,length):
    s=convertBytetoBinaryHelper(x,'')
    s=s.ljust(length,"0")
    #print s
    return s

def convertBytetoBinaryHelper(x,s): #convert byte x to a binary string
    binary=['0','1']
    m=int(x%2)
    q=int(x/2)
    s=s+binary[m]
    if q >0:
        return convertBytetoBinaryHelper(q,s)
    else:
        return s

def convertBase(base, dec, length, symbols=bases):
    assert base <= len(symbols)
    s = convertBaseHelper(base,dec,'', bases)
    s = s.ljust(length,bases[0])
    return s

def convertToAnyBase(base : int, dec : int, length: int, symbols=bases):
    assert base == len(symbols)
    s = convertBaseHelper(base,dec,'', symbols)
    while len(s) < length:
        s += symbols[0]
    assert length == len(s)
    return s

def convertFromBase(base,s):
    val = 0
    power = 1
    for i in range(len(s)):
        val += values[s[i]]*power
        power *= base
    return val

def convertIntToBytes(val,num_bytes):
    val = int(val)
    if val == None:
        return [-1 for _ in range(num_bytes)]
    else:
        l = [(val & (0xff << pos*8)) >> pos*8 for pos in range(num_bytes)]
        return l

def convertBytesToInt(l):
    _sum = 0
    for i,val in enumerate(l):
        _sum += val * (256**i)
    return _sum

