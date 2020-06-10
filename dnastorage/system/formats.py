from dnastorage.codec import commafreecodec
from dnastorage.codec.builder import customize_RS_CFC8
from dnastorage.exceptions import *

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

FileSystemFormats = {
    # KEY      KEY    LEN  PacketSize, Abbrev.   Description                   ENCODER       DECODER
    0x0010 : [0x0010, 200, 90, "FSMD", "File system meta-data format", ENC_FSMD_200, DEC_FSMD_200 ],
    0x0020 : [0x0020, 200, 16, "RS+CFC8", "Reed-Solomon coded with ad hoc Comma-free codewords",
              ENC_RS_CFC8_200, DEC_RS_CFC8_200 ],
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



