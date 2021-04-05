from io import BytesIO
import sys
import logging

from dnastorage.util.packetizedfile import ReadPacketizedFilestream, WritePacketizedFilestream
import dnastorage.system.formats as formats
import dnastorage.system.header as header
from dnastorage.codec.base_conversion import convertIntToBytes, convertBytesToInt

logger = logging.getLogger('dna.storage.system.dnafile')
logger.addHandler(logging.NullHandler())

def reverse_complement(seq):
    complement = {'T':'A', 'G':'C', 'C':'G', 'A':'T'}
    r = [complement[x] for x in seq]
    r.reverse()
    return "".join(r)

class DNAFile:
    def __init__(self):
        return

    @classmethod
    def open(self, filename, op, primer5, primer3, format_name="", write_incomplete_file=False, fsmd_abbrev='FSMD',flanking_primer5="",flanking_primer3="",use_flanking_primer_for_decoding=False,use_single_primer=False,preview_mode=False,reverse_primer3_from_seq=False):
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

                if h['formatid'] == formats.file_system_formatid_by_abbrev("Segmented"):
                    logger.debug("SegmentedReadDNAFile({},{},{})".format(filename,primer5,primer3))
                    return SegmentedReadDNAFile(input=filename,\
                                                primer5=primer5,primer3=primer3,\
                                                write_incomplete_file=write_incomplete_file,\
                                                fsmd_abbrev=fsmd_abbrev,\
                                                flanking_primer5=flanking_primer5,\
                                                flanking_primer3=flanking_primer3,\
                                                use_flanking_primer_for_decoding=use_flanking_primer_for_decoding,\
                                                use_single_primer=use_single_primer,\
                                                header=h,\
                                                preview_mode=preview_mode,\
                                                reverse_primer3_from_seq=reverse_primer3_from_seq)
                else:
                    return ReadDNAFile(input=filename, primer5=primer5, primer3=primer3,\
                                       fsmd_abbrev=fsmd_abbrev,\
                                       flanking_primer5=flanking_primer5,\
                                       flanking_primer3=flanking_primer3,\
                                       use_flanking_primer_for_decoding=use_flanking_primer_for_decoding,\
                                       header=h,preview_mode=preview_mode)
                
                # return ReadDNAFile(input=filename, primer5=primer5, primer3=primer3,\
                #                    fsmd_abbrev=fsmd_abbrev,\
                #                    flanking_primer5=flanking_primer5,\
                #                    flanking_primer3=flanking_primer3,\
                #                    use_flanking_primer_for_decoding=use_flanking_primer_for_decoding)
        elif "w" in op and "s" in op:
            return SegmentedWriteDNAFile(output=filename,primer5=primer5,primer3=primer3, \
                                         format_name=format_name,fsmd_abbrev=fsmd_abbrev,\
                                         flanking_primer5=flanking_primer5,\
                                         flanking_primer3=flanking_primer3)       
        elif op=="w":
            return WriteDNAFile(output=filename,primer5=primer5,primer3=primer3, \
                                format_name=format_name,fsmd_abbrev=fsmd_abbrev,\
                                flanking_primer5=flanking_primer5,\
                                flanking_primer3=flanking_primer3)
                
        # elif op=="w":
        #     return WriteDNAFile(output=filename,primer5=primer5,primer3=primer3, \
        #                         format_name=format_name,fsmd_abbrev=fsmd_abbrev,\
        #                         flanking_primer5=flanking_primer5,\
        #                         flanking_primer3=flanking_primer3)
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
            self.need_to_close = True
        elif 'in_fd' in kwargs:
            self.in_fd = kwargs['in_fd']
            self.input_filename = ""
            self.need_to_close = False
            
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
        if 'header' in kwargs and kwargs['header'] != None:
            h = kwargs['header']
        else:
            h = header.decode_file_header(strands,self.primer5,self.primer3,fsmd_abbrev=self.fsmd_abbrev)    

        self.strands = header.pick_nonheader_strands(strands,self.primer5)
        
        assert h['version'][0] <= header.system_version['major']
        assert h['version'][1] <= header.system_version['minor']

        self.header = h 
        self.formatid = h['formatid']
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
        
        if self.need_to_close:
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
            self.out_fd = open(self.output_filename,"wt")
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
        try:# fix me
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


class SegmentedWriteDNAFile(WriteDNAFile):
    # SegmentedWriteDNAFile writes a set of strands.
    def __init__(self,**kwargs):     
        WriteDNAFile.__init__(self,**kwargs)
        self.segments = []
        self.beginIndex = 0
        return

    def _record_segment(self):
        self.segments += [[ self.formatid, self.size, self.primer5, self.primer3, self.beginIndex, self.flanking_primer5, self.flanking_primer3 ]]

    def new_segment(self, format_name, primer5, primer3, flanking_primer5="", flanking_primer3=""):
        self._encode_buffer()  # write everything in the buffer to the file
        self._record_segment() # remember new segment
        self.primer5 = primer5
        self.primer3 = primer3
        self.flanking_primer5 = flanking_primer5
        self.flanking_primer3 = flanking_primer3
        # get current index + 1
        self.beginIndex = self.enc.index
        #print "beginIndex={}".format(self.enc.index)
        enc_func = formats.file_system_encoder_by_abbrev(format_name)
        self.formatid = formats.file_system_formatid_by_abbrev(format_name)
        # we consumed the prior buffer, so just make a new one to avoid
        # cursor positioning problems (JMT: not sure if this is the best way)
        self.mem_buffer = BytesIO()
        self.pf = ReadPacketizedFilestream(self.mem_buffer)
        self.enc = enc_func(self.pf,flanking_primer5+primer5,flanking_primer3+primer3,self.beginIndex)
        self.size=0

    def encode_segments_header(self,segments):
        oprimer5 = segments[0][2]
        oprimer3 = segments[0][3]
        assert len(segments) <= 256
        #if len(segments) == 0:
        #    return []
        hdr = [ len(segments) ]
        logger.debug("encode segments header : {}".format(hdr))
        for s in segments:
            hdr += convertIntToBytes(s[0],2)
            hdr += header.encode_size_and_value( s[1] )
            hdr += header.encode_size_and_value( s[4] )
            hdr += header.encode_primer_diff(oprimer5,s[2])
            hdr += header.encode_primer_diff(oprimer3,s[3])

        return hdr

    def encode_segments_header_comments(self,segments):
        comment = "% segment descriptions\n"        
        tab = "%    "
        for i,s in enumerate(segments):
            comment += tab + "{}. ".format(i) + formats.file_system_format_description(s[0]) + "\n"
            comment += tab + "    size = {}".format(s[1]) + "\n"
            comment += tab + "    5' = {}".format(s[2]) + "\n"
            comment += tab + "    3' = {}".format(s[3]) + "\n"
            comment += tab + "    beginIndex = {}".format(s[4]) + "\n"            
            comment += "%\n"
        return comment

    def close(self):
        logger.debug("WriteSegmentedDNAFile.close")
        
        self.flush()
        self._record_segment() # record last segment
        
        hdr_other = self.encode_segments_header(self.segments)

        formatid = formats.file_system_formatid_by_abbrev("Segmented")

        #print "formatid=",formatid
        
        size = sum([x[1] for x in self.segments])
        primer5 = self.segments[0][2]
        primer3 = self.segments[0][3]
        flanking5 = self.segments[0][-2]
        flanking3 = self.segments[0][-1]
        
        hdr = header.encode_file_header(self.output_filename,formatid,size,hdr_other,flanking5+primer5,flanking3+primer3,fsmd_abbrev=self.fsmd_abbrev)

        #print "Number of strands in header: ", len(hdr)

        for i,h in enumerate(hdr):
            self.strands.insert(i,h)

        comment = header.encode_file_header_comments(self.output_filename,formatid,\
                                              size,hdr_other,primer5,primer3)
        self.out_fd.write(comment)
        comment = self.encode_segments_header_comments(self.segments)
        self.out_fd.write(comment)
        for ss in self.strands:
            if type(ss) is list:
                for s in ss:
                    self.out_fd.write("{}\n".format(s))
            else:
                self.out_fd.write("{}\n".format(ss))
                    
        if self.out_fd != sys.stdout and self.out_fd != sys.stderr:
            self.out_fd.close()
        return


class SegmentedReadDNAFile(ReadDNAFile):

    def decode_segments_header(self,other_data):
        val = other_data
        numSeg = val[0]
        if numSeg==0:
            return
        pos = 1
        allSegs = []
        for i in range(numSeg):
            # get format of this segment
            seg = [convertBytesToInt(val[pos:pos+2])]
            pos+=2

            # get size in bytes
            v,p = header.decode_size_and_value(val,pos)
            pos += p
            seg += [v]

            # get begin index
            v,p = header.decode_size_and_value(val,pos)
            pos += p
            seg += [v]

            primer5,p = header.decode_primer_diff(val[pos:], self.primer5)

            pos += p
            primer3,p = header.decode_primer_diff(val[pos:], self.primer3)

            if i > 0 and self.reverse_primer3_from_seq:
                primer3 = reverse_complement(primer3)
                
            pos += p
            seg.append(primer5)
            seg.append(primer3)
            allSegs.append(seg)

        return allSegs



    # ReadDNAFile reads a set of strands from a file.  It finds the header,
    # determines compatibility and encoding type, and then decodes the file.
    #    
    def __init__(self,**kwargs):     
        ReadDNAFile.__init__(self,**kwargs)

        if 'reverse_primer3_from_seq' in kwargs:
            self.reverse_primer3_from_seq = kwargs['reverse_primer3_from_seq']
        else:
            self.reverse_primer3_from_seq = False
        
        logger.debug("sizeof other_data = {}".format(len(self.header['other_data'])))
        
        if len(self.header['other_data'])==0:
            return

        # restore cursor to end of buffer for writing
        self.mem_buffer.seek(0,2)
        
        segs = self.decode_segments_header(self.header['other_data'])
        self.segments = segs

        if 'preview_mode' in kwargs:
            self.preview_mode = kwargs['preview_mode']
        else:
            self.preview_mode = False
            
        #print "segments=",segs

        for s in segs:
            logger.debug("formatid={} size={} bindex={} primer5={} primer3={}".format(s[0],s[1],s[2],s[3],s[4]))                                             
            formatid = s[0]
            size = s[1]
            bindex = s[2]
            primer5 = s[3]
            primer3 = s[4]

            logger.info("{}".format(s))
            
            if 'use_single_primer' in kwargs and kwargs['use_single_primer']==True:
                # strands from sequencing should all have the same primer
                primer5 = self.primer5
                primer3 = self.primer3

            dec_func = formats.file_system_decoder(formatid)
            #self.mem_buffer = BytesIO()
            self.pf = WritePacketizedFilestream(self.mem_buffer,size,\
                                                formats.file_system_format_packetsize(formatid),\
                                                minKey=bindex)

            #print primer5, primer3, bindex
            self.dec = dec_func(self.pf,primer5,primer3,bindex)

            for s in self.strands:
                #if s.find(primer5)!=-1:
                #print "dnafile.py",self.dec.decode_from_phys_to_strand(s)
                self.dec.decode(s)

            self.dec._attempt_final_decoding()
            #print(self.dec.strand_errors,self.dec.block_errors," missing=",self.dec._packetizedFile.hasMissingKeys())
                            
            if self.preview_mode and (self.dec.block_errors > 0 or not self.dec.complete):
                print("found errors!",self.dec.block_errors)
                print("mising=",self.dec._packetizedFile.hasMissingKeys())
                break
            else:
                self.dec.only_write()
                write_anyway = 'write_incomplete_file' in kwargs and \
                    kwargs['write_incomplete_file']==True
                #if self.dec.complete or write_anyway:
                #    self.dec.write()
                if not write_anyway:
                    assert self.dec.complete

        #print [x for x in self.mem_buffer.getvalue()]
        self.mem_buffer.seek(0,0) # set read point at beginning of buffer

        if self.need_to_close:
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
