from dnastorage.exceptions import NoTolerance, AllowAll
from dnastorage.codec import commafreecodec
from dnastorage.codec.builder import customize_RS_CFC8

def ENC_FSMD_200(pf, primer5, primer3, bIndex=0, policy=NoTolerance(),withCut=None):
    enc = customize_RS_CFC8(True,pf,primer5,primer3,1,0,2,15,90,policy,withCut=withCut,outerECCStrands=20,minIndex=bIndex)
    return enc    

def DEC_FSMD_200(pf, primer5, primer3, bIndex=0,policy=AllowAll(),withCut=None):
    dec = customize_RS_CFC8(False,pf,primer5,primer3,1,0,2,15,90,policy,withCut=withCut,outerECCStrands=20,minIndex=bIndex)
    return dec

def ENC_RS_CFC8_200(pf, primer5, primer3, bIndex=0, policy=NoTolerance()):
    enc = customize_RS_CFC8(True,pf,primer5,primer3,1,2,2,15,185*15,policy,minIndex=bIndex)
    return enc    

def DEC_RS_CFC8_200(pf, primer5, primer3, bIndex=0, policy=AllowAll()):
    dec = customize_RS_CFC8(False,pf,primer5,primer3,1,2,2,15,185*15,policy,minIndex=bIndex)
    return dec    

### Support for preview - FIXME - move to dnapreview - create registry

def ENC_FSMD_WCUT_160(pf, primer5, primer3, bIndex=0, policy=NoTolerance(),withCut=None):
    withCut="AGGTACCA"
    enc = customize_RS_CFC8(True,pf,primer5,primer3,1,0,2,10,120,policy,withCut=withCut,outerECCStrands=20,minIndex=bIndex)
    return enc    

def DEC_FSMD_WCUT_160(pf, primer5, primer3, bIndex=0,policy=AllowAll(),withCut=None):
    withCut="AGGTACCA"
    dec = customize_RS_CFC8(False,pf,primer5,primer3,1,0,2,10,120,policy,withCut=withCut,outerECCStrands=20,minIndex=bIndex)
    return dec

def ENC_RS_CFC8_RE1_160(pf, primer5, primer3, bIndex=0, policy=NoTolerance()):
    withCut="AGGTACCA"
    enc = customize_RS_CFC8(True,pf,primer5,primer3,1,1,3,9,185*9,policy,minIndex=bIndex,
                            withCut=withCut)
    return enc    

def DEC_RS_CFC8_RE1_160(pf, primer5, primer3, bIndex=0, policy=AllowAll()):
    withCut="AGGTACCA"
    dec = customize_RS_CFC8(False,pf,primer5,primer3,1,1,3,9,185*9,policy,minIndex=bIndex,\
                            withCut=withCut)
    return dec    

def ENC_RS_CFC8_RE2_160(pf, primer5, primer3, bIndex=0, policy=NoTolerance()):
    withCut="CCTGCAGG"
    enc = customize_RS_CFC8(True,pf,primer5,primer3,1,1,3,9,185*9,policy,minIndex=bIndex,
                            withCut=withCut)
    return enc    

def DEC_RS_CFC8_RE2_160(pf, primer5, primer3, bIndex=0, policy=AllowAll()):
    withCut="CCTGCAGG"
    dec = customize_RS_CFC8(False,pf,primer5,primer3,1,1,3,9,185*9,policy,minIndex=bIndex,\
                            withCut=withCut)
    return dec    

def ENC_RS_CFC8_RE3_160(pf, primer5, primer3, bIndex=0, policy=NoTolerance()):
    withCut="GCGGCCGC"
    enc = customize_RS_CFC8(True,pf,primer5,primer3,1,1,3,9,185*9,policy,minIndex=bIndex,
                            withCut=withCut)
    return enc    

def DEC_RS_CFC8_RE3_160(pf, primer5, primer3, bIndex=0, policy=AllowAll()):
    withCut="GCGGCCGC"
    dec = customize_RS_CFC8(False,pf,primer5,primer3,1,1,3,9,185*9,policy,minIndex=bIndex,\
                            withCut=withCut)
    return dec    

def ENC_RS_CFC8_RE4_160(pf, primer5, primer3, bIndex=0, policy=NoTolerance()):
    withCut="GTTTAAAC"
    enc = customize_RS_CFC8(True,pf,primer5,primer3,1,1,3,9,185*9,policy,minIndex=bIndex,
                            withCut=withCut)
    return enc    

def DEC_RS_CFC8_RE4_160(pf, primer5, primer3, bIndex=0, policy=AllowAll()):
    withCut="GTTTAAAC"
    dec = customize_RS_CFC8(False,pf,primer5,primer3,1,1,3,9,185*9,policy,minIndex=bIndex,\
                            withCut=withCut)
    return dec    

### End support for preview - FIXME - move to dnapreview - create registry

FileSystemFormats = {
    # KEY      KEY    LEN  PacketSize, Abbrev.   Description                   ENCODER       DECODER
    0x0010 : [0x0010, 200, 90, "FSMD", "File system meta-data format", ENC_FSMD_200, DEC_FSMD_200 ],
    0x0020 : [0x0020, 200, 16, "RS+CFC8", "Reed-Solomon coded with ad hoc Comma-free codewords",
              ENC_RS_CFC8_200, DEC_RS_CFC8_200 ],

    # Support preview -- support through registry
    0x0011 : [0x0011, 160, 120, "FSMD-1", "File system meta-data format with cut", ENC_FSMD_WCUT_160, DEC_FSMD_WCUT_160 ],
    0x1000 : [0x1000, 200, 20, "Segmented", "Segmented file format to support Preview", None, None],
    0x2021 : [0x2021, 160, 9, "RS+CFC8+RE1", "Reed-Solomon coded with Comma-free codewords",
              ENC_RS_CFC8_RE1_160, DEC_RS_CFC8_RE1_160 ],
    0x2022 : [0x2022, 160, 9, "RS+CFC8+RE2", "Reed-Solomon coded with Comma-free codewords",
              ENC_RS_CFC8_RE2_160, DEC_RS_CFC8_RE2_160 ],
    0x2023 : [0x2023, 160, 9, "RS+CFC8+RE3", "Reed-Solomon coded with Comma-free codewords",
              ENC_RS_CFC8_RE3_160, DEC_RS_CFC8_RE3_160 ],
    0x2024 : [0x2024, 160, 9, "RS+CFC8+RE4", "Reed-Solomon coded with Comma-free codewords",
              ENC_RS_CFC8_RE4_160, DEC_RS_CFC8_RE4_160 ],
}

def file_system_formats():
    return [ v[3] for k,v in FileSystemFormats.items() ]

_abbrevFileSystemDict = { v[3] : v for k,v in FileSystemFormats.items() }

def file_system_format_description(formatid):
    return FileSystemFormats[formatid][4]

def file_system_format_packetsize(formatid):
    return FileSystemFormats[formatid][2]

def file_system_encoder(formatid):
    return FileSystemFormats[formatid][5]

def file_system_decoder(formatid):
    return FileSystemFormats[formatid][6]

def file_system_encoder_by_abbrev(ab):
    return _abbrevFileSystemDict[ab][5]

def file_system_decoder_by_abbrev(ab):
    return _abbrevFileSystemDict[ab][6]

def file_system_formatid_by_abbrev(ab):
    return _abbrevFileSystemDict[ab][0]



