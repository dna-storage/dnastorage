from io import BytesIO
import sys
import logging

from dnastorage.util.packetizedfile import ReadPacketizedFilestream, WritePacketizedFilestream
import dnastorage.system.formats as formats
import dnastorage.system.header as header

logger = logging.getLogger('dna.storage.system.dnafile')
logger.addHandler(logging.NullHandler())

class DNAFile:
    def __init__(self):
        return

    @classmethod
    def open(self, filename, op, primer5, primer3, format_name="", write_incomplete_file=False, fsmd_abbrev='FSMD',flanking_primer5="",flanking_primer3="",use_flanking_primer_for_decoding=False):
        '''        
        1. filename is the input file of strands
        2. op is the file operation: 'r' or 'w'
        3. format_name is optional for reading and required for writing. 
        4. primer5 is the coding primer, and primer3 is the non-coding primer.
        '''     
        # check if we are reading or writing
        if op=="r":
            
            with open(filename,"r") as fd:
                logger.debug("open {} for reading.".format(filename))
                strands = get_strands(fd)

                if use_flanking_primer_for_decoding==True:
                    h = header.decode_file_header(strands,flanking_primer5+primer5,\
                                       flanking_primer3+primer3,fsmd_abbrev=fsmd_abbrev)
                else:
                    h = header.decode_file_header(strands,primer5,primer3,fsmd_abbrev=fsmd_abbrev)
            
                logger.debug("decoded header: {}".format(h)) 
                assert h['version'][0] <= header.system_version['major']
                assert h['version'][1] <= header.system_version['minor']
            
                return ReadDNAFile(input=filename, primer5=primer5, primer3=primer3,\
                                   fsmd_abbrev=fsmd_abbrev,\
                                   flanking_primer5=flanking_primer5,\
                                   flanking_primer3=flanking_primer3,\
                                   use_flanking_primer_for_decoding=use_flanking_primer_for_decoding)

        elif op=="w":
            return WriteDNAFile(output=filename,primer5=primer5,primer3=primer3, \
                                format_name=format_name,fsmd_abbrev=fsmd_abbrev,\
                                flanking_primer5=flanking_primer5,\
                                flanking_primer3=flanking_primer3)
        else:
            return None

    def flush(self):
        return

    def close(self):
        return

    def readable(self):
        assert False
    def writable(self):
        assert False

def get_strands(in_fd):
    strands = []
    while True:
        s = in_fd.readline()        
        if len(s) == 0:
            break
        s = s.strip()
        if s.startswith('%'):
            continue
        strands.append(s)
    return strands

class ReadDNAFile(DNAFile):
    # ReadDNAFile reads a set of strands from a file.  It finds the header,
    # determines compatibility and encoding type, and then decodes the file.
    #
    def __init__(self,**kwargs):     
        DNAFile.__init__(self)

        if 'input' in kwargs:
            self.input_filename = kwargs['input']
            self.in_fd = open(self.input_filename,"r")
            need_to_close = True
        elif 'in_fd' in kwargs:
            self.in_fd = kwargs['in_fd']
            self.input_filename = ""
            need_to_close = False
            
        if not ('fsmd_abbrev' in kwargs):
            self.fsmd_abbrev = 'FSMD'
        else:
            self.fsmd_abbrev = kwargs['fsmd_abbrev']
            
        assert 'primer5' in kwargs and 'primer3' in kwargs
        self.primer5 = kwargs['primer5']
        self.primer3 = kwargs['primer3']

        if 'flanking_primer5' in kwargs:
            self.flanking_primer5 = kwargs['flanking_primer5']
        else:
            self.flanking_primer5 = ''

        if 'flanking_primer3' in kwargs:
            self.flanking_primer3 = kwargs['flanking_primer3']
        else:
            self.flanking_primer3 = ''

        if 'use_flanking_primer_for_decoding' in kwargs:
            self.use_flanking_primers = kwargs['use_flanking_primer_for_decoding']
        else:
            self.use_flanking_primers = False
            
        # get all the strands ( this will be bad for large files )
        strands = get_strands(self.in_fd)

        h = header.decode_file_header(strands,self.primer5,self.primer3,fsmd_abbrev=self.fsmd_abbrev)

        self.strands = header.pick_nonheader_strands(strands,self.primer5)
        
        assert h['version'][0] <= header.system_version['major']
        assert h['version'][1] <= header.system_version['minor']

        self.formatid = h['formatid']
        self.header = h 
        self.size = h['size']

        # set up mem_buffer 
        self.mem_buffer = BytesIO()
        
        if self.formatid == 0x1000:
            # let sub-classes handle initialization
            return

        dec_func = formats.file_system_decoder(self.formatid)
            
            
        self.mem_buffer = BytesIO()
        self.pf = WritePacketizedFilestream(self.mem_buffer,self.size,formats.file_system_format_packetsize(self.formatid))

        if self.use_flanking_primers:
            self.dec = dec_func(self.pf,self.flanking_primer5+self.primer5,\
                                self.flanking_primer3+self.primer3)
        else:
            self.dec = dec_func(self.pf,kwargs['primer5'],kwargs['primer3'])

        for s in self.strands:
            if s.startswith(self.primer5):
                self.dec.decode(s)

        self.dec.write()
        assert self.dec.complete

        self.mem_buffer.seek(0,0) # set read point at beginning of buffer
        
        if need_to_close:
            self.in_fd.close()
        return

    def read(self, n=1):        
        return self.mem_buffer.read(n)
    def readline(self, n=-1):
        return self.mem_buffer.readline(n)

    def readable(self):
        return True
    def writable(self):
        return False
        
class WriteDNAFile(DNAFile):
    # WriteDNAFile writes a set of strands.
    def __init__(self,**kwargs):     
        DNAFile.__init__(self)
        if 'formatid' in kwargs:            
            enc_func = formats.file_system_encoder(kwargs['formatid'])
            self.formatid = kwargs['formatid']
        elif 'format_name' in kwargs:
            enc_func = formats.file_system_encoder_by_abbrev(kwargs['format_name'])
            self.formatid = formats.file_system_formatid_by_abbrev(kwargs['format_name'])

        if not ('fsmd_abbrev' in kwargs):
            self.fsmd_abbrev = 'FSMD'
        else:
            self.fsmd_abbrev = kwargs['fsmd_abbrev']

        self.mem_buffer = BytesIO()
        self.pf = ReadPacketizedFilestream(self.mem_buffer)

        assert 'primer5' in kwargs and 'primer3' in kwargs
        self.primer5 = kwargs['primer5']
        self.primer3 = kwargs['primer3']

        if 'flanking_primer3' in kwargs:
            self.flanking_primer3 = kwargs['flanking_primer3']
        else:
            self.flanking_primer3 = ''

        if 'flanking_primer5' in kwargs:
            self.flanking_primer5 = kwargs['flanking_primer5']
        else:
            self.flanking_primer5 = ''

        self.enc = enc_func(self.pf,self.flanking_primer5+self.primer5,\
                            self.flanking_primer3+self.primer3)

            
        if 'output' in kwargs:
            self.output_filename = kwargs['output']
            self.out_fd = open(self.output_filename,"w")
        elif 'out_fd' in kwargs:
            self.out_fd = kwargs['out_fd']
            self.output_filename = ""
            
        self.size = 0
        self.strands = []
        return

    def _encode_buffer(self):
        for block in self.enc:
            if type(block) == list:
                for s in block:
                    self.strands.append(s)
            else:
                self.strands.append(block)
        #print "_encode_buffer: index={}".format(self.enc.index)

    def read(self, size):
        assert False and "Do not read at the same time as writing"
        
    def writable(self):
        return True
    def readable(self):
        return False
    
    def write(self, buff):
        self.size += len(buff)
        tell = self.mem_buffer.tell()
        self.mem_buffer.seek(0,2)

        #print ( len(buff), bytes([x for x in buff]), buff )

        #buff = bytes([x for x in buff])        
        #buff = bytes([ord(x) for x in buff])
        try:
            # convert string
            self.mem_buffer.write(buff)
            self.mem_buffer.seek(tell,0)
        except Exception as e:
            buff = bytearray([x for x in buff])
            self.mem_buffer.write(buff)
            self.mem_buffer.seek(tell,0)
            
        return
    
    def flush(self):
        self._encode_buffer()
        return

    def close(self):
        self.flush()
        hdr = header.encode_file_header(self.output_filename,self.formatid,self.size,\
                                 "",self.flanking_primer5+self.primer5,\
                                 self.flanking_primer3+self.primer3,\
                                 fsmd_abbrev=self.fsmd_abbrev)

        for i,h in enumerate(hdr):
            self.strands.insert(i,h)

        comment = header.encode_file_header_comments(self.output_filename,self.formatid,self.size,"",\
                                              self.primer5,self.primer3)
        self.out_fd.write(comment)
        for s in self.strands:
            self.out_fd.write("{}\n".format(s))

        if self.out_fd != sys.stdout and self.out_fd != sys.stderr:
            self.out_fd.close()
        return


    
