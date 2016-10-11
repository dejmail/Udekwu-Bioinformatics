import os

import logging
import logging.config
import yaml

from string import upper
from re import compile
import argparse

from cogent.parse.fastq import MinimalFastqParser
from skbio.sequence import DNA
 
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
#from Bio.SeqRecord import SeqRecord    
from StringIO import StringIO
        
from qiime.check_id_map import process_id_map
import subprocess

class Demultiplex(object):
    
    def __init__(self, logger=None):
        """
        """
        self.setup_logging()
        self.logger = logger or logging.getLogger(__name__)
        self.f_count = 0
        self.r_count = 0
        self.f_only_count = 0
        self.r_only_count = 0
        self.both_primers_count = 0
        self.no_seq_left = 0
        self.quality_errors = 0
        self.unmapped_count = 0
        self.total_seqs = 0
        print("Instantiating class")        
        
        #self.name = name
        #self.balance = balance

    def setup_logging(self,
                      default_path='logging.yaml',
                      default_level=logging.INFO,
                      env_key='LOG_CFG'):
        
        """Setup logging configuration"""
        
        path = default_path
        value = os.getenv(env_key, None)
        if value:
            path = value
        if os.path.exists(path):
            with open(path, 'rt') as f:
                config = yaml.safe_load(f.read())
            logging.config.dictConfig(config)
        else:
            logging.basicConfig(level=default_level)


    def check_if_files_exist(self, file_list):
        
        for files in file_list:
            try: 
                self.logger.info(type(files))
            except (OSError, IOError) as e:
                self.logger.error("The filepath {0} does not exist, cannot be found, or cannot be read".format(e))
            
                
                
    def get_primers(self, header, mapping_data):
        """ Returns lists of forward/reverse primer regular expression
    
            header:  list of strings of header data.
            mapping_data:  list of lists of mapping data
    
            Will raise error if either the LinkerPrimerSequence or ReversePrimer fields
            are not present
        """
     
        if "LinkerPrimerSequence" in header:
            primer_ix = header.index("LinkerPrimerSequence")
        else:
            raise IndexError(
                ("Mapping file is missing LinkerPrimerSequence field."))
        if "ReversePrimer" in header:
            rev_primer_ix = header.index("ReversePrimer")
        else:
            raise IndexError(("Mapping file is missing ReversePrimer field."))
     
        iupac = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': '[AG]', 'Y': '[CT]',
                 'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]',
                 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'}
     
        raw_forward_primers = set([])
     
        raw_reverse_primers = set([])
    
        for line in mapping_data:
            # Split on commas to handle pool of primers
            raw_forward_primers.update([upper(primer).strip() for
                                        primer in line[primer_ix].split(',')])
            raw_reverse_primers.update([upper(str(DNA(primer).rc())) for
                                        primer in line[rev_primer_ix].split(',')])
    
        if not raw_forward_primers:
            self.logger.critical("No forward primers detected in mapping file.")
            raise ValueError(("No forward primers detected in mapping file."))
            
        if not raw_reverse_primers:
            self.logger.critical("No reverse primers detected in mapping file.")
            raise ValueError(("No reverse primers detected in mapping file."))

     
        forward_primers = []
        reverse_primers = []
        for curr_primer in raw_forward_primers:
            forward_primers.append(compile(''.join([iupac[symbol] for
                                                    symbol in curr_primer])))
        for curr_primer in raw_reverse_primers:
            reverse_primers.append(compile(''.join([iupac[symbol] for
                                                    symbol in curr_primer])))
         
        return forward_primers, reverse_primers
     
            
    def run_demultiplex_and_trim(self, opts):
        sample_primer_dict = {}
        file_list = []
    
        self.logger.info("Starting counts \n\
                         f_count - {0} \n\
                         r_count {1} \n\
                         f_only_count {2} \n\
                         r_only_count {3} \n\
                         both_primers_count {4} \n\
                         no_seq_left {5} \n\
                         quality_errors {6} \n\
                         unmapped_count {7} \n\
                         total_seqs {8} \n".format(self.f_count, self.r_count,
                                                self.f_only_count, 
                                                self.r_only_count,
                                                self.both_primers_count,
                                                self.no_seq_left, 
                                                self.quality_errors,
                                                self.unmapped_count, 
                                                self.total_seqs))

        metafile = opts.m
        sequence_file = opts.f
       # log_out = opts.o
        
        #metafile = opts['m']
        #sequence_file = opts['f']
        #log_out = open(opts['o'], 'w')
    
        #extract the relevant data from the metadata file
        header, mapping_data, run_description, errors, warnings = process_id_map(metafile)
        # get the primer search patterns
        forward_primers, reverse_primers = self.get_primers(header, mapping_data)
    
        # replace colons with underscores in the sample_id names
        for sample_id in mapping_data:
            sample_primer_dict[sample_id[0].replace(":","_")] = (sample_id[2], sample_id[4])
            file_list.append(sample_id[0].replace(":","_"))
        self.logger.info(mapping_data)
            
        for label,seq,qual in MinimalFastqParser(sequence_file, strict=False):
            self.total_seqs += 1
            start_slice = 0
            end_slice = -1
            f_primer_found = False
            r_primer_found = False
            for curr_primer in forward_primers:
                self.logger.info("current F primer pattern...{0}".format(curr_primer.pattern))        
                # regex search using curr_primer pattern
                if curr_primer.search(seq):
                    f_primer_seq = curr_primer
                    start_slice = int(curr_primer.search(seq).span()[1])
                    self.f_count += 1
                    f_primer_found = True
            for curr_primer in reverse_primers:
                self.logger.info("current R primer pattern...{0}".format(curr_primer.pattern))
                if curr_primer.search(seq):
                    r_primer_seq = curr_primer
                    # span() gives start and end coordinates
                    end_slice = int(curr_primer.search(seq).span()[0])
                    self.r_count += 1
                    r_primer_found = True
            curr_seq = seq[start_slice:end_slice]
            # need to clip the quality sequence as well
            curr_qual = qual[start_slice:end_slice]
            if len(curr_seq) < 1:
                self.no_seq_left += 1
                continue
    
            if f_primer_found == True and r_primer_found == False:
                
                self.f_only_count += 1
    
            if f_primer_found == False and r_primer_found == True:
                self.r_only_count += 1
    
            if f_primer_found and r_primer_found:
                self.both_primers_count += 1
            
                
                # write each sequence that has a forward and reverse primer to file
                
                filename = self.get_sample_id_from_primer_sequence(sample_primer_dict, f_primer_seq.pattern, r_primer_seq.pattern)
                filename = str(filename) + ".fastq"
                #construct a SeqIO object, which checks for proper quality scores and formatting
                record = self.return_fastq_seqio_object(seq,label,qual,filename)
                # each file is appended to, not overwritten
                out_seqs = open(filename, 'a')
                if filename == "None.fastq":
                    self.unmapped_count += 1
                    out_seqs.write("@{0}\n{1}\n+{0}\n{2}\n".format(label+":sample_id_"+filename[:-6], curr_seq, curr_qual))
                else:            
                    try:
                        #logger.info(filename)
                        SeqIO.write(record, out_seqs, "fastq")
                    except TypeError as self.e:
                        self.quality_errors += 1
                        #logger.info(record, e.args, e.message)
    
                out_seqs.close()
    
        self.logger.info("Forward primer hits: {0}\n".format(self.f_count))
        self.logger.info("Only forward primer found: {0}\n".format(self.f_only_count))
        self.logger.info("Reverse primer hits: {0}\n".format(self.r_count))
        self.logger.info("Only reverse primer found: {0}\n".format(self.r_only_count))
        self.logger.info("No seq left after truncation: {0}\n".format(self.no_seq_left))
        self.logger.info("Sequences with errors in quality scores: {0}\n".format(self.quality_errors))
        self.logger.info("Sequences not mapped: {0}\n".format(self.unmapped_count))
        self.logger.info("Total sequences checked: {0}\n".format(self.total_seqs))
        unaccounted_sequences = (self.total_seqs -
                                self.f_count -
                                self.r_count -
                                self.r_only_count -
                                self.f_only_count -
                                self.no_seq_left -
                                self.quality_errors -
                                self.unmapped_count)
        self.logger.info("Sequences unaccounted for: {0}".format(unaccounted_sequences))
    
        self.logger.info("Run finished")
    

    
    
    # In[47]:
    
    def check_metadata_file_formatting(self, metafile):
    
        command = ("validate_mapping_file.py -m {0}".format(metafile))
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell='TRUE')
        output, err = p.communicate()
        
        return (output,err)
    
    def process_metadata_file(self, metafile):
    
        try:
            output = self.check_metadata_file_formatting(str(metafile))
            self.logger.info("processed ok - {0}".format(output[0]))
        except IOError as e:
            self.logger.error("The metadata file does not exist or the path is wrong - error {0}".format(e))
        
    
    
    # To be able to match the sample_id and the associated primer names and sequences with the correct reads, I need to invert the regex pattern. This is easily accomplished by making a new iupac dictionary and simply swapping the keys and values of the original dictionary. 
    
    # In[48]:
    
    iupac = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': '[AG]', 'Y': '[CT]',
             'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]',
             'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'}
    
    def invert_regex_pattern(self, sample_primer_dict):
        ''' Inverts the dictionary of iupac DNA ambiguity letters so that we can
        check the regex pattern with something else'''
        inverted_dict = {}
        
        for keys, values in sample_primer_dict.items():
            # only keep the ambiguous combinations
            if len(values) > 1:
                inverted_dict.update({values: keys})
        return inverted_dict
    
    #logger.info(inverted_iupac)
    
    
    def replace_ambiguous_pattern_with_iupac_base(self, sequence):
        
        self.logger.info("reversing pattern {0}".format(sequence))
        #inverted_iupac = self.invert_regex_pattern(sequence)
        #self.logger.info("reversed {0} pattern".format(inverted_iupac))
        #self.sequence = SeqIO.read(StringIO(sequence), "fasta")
        for key,value in inverted_iupac.items():
            if key in sequence:
                sequence = sequence.replace(key, value, 1)
            else:
                # loop does not continue if this is not here, strange?
                continue
        return str(sequence)
    
    
    # Simply reverse complement the primer sequence to reflect the fact that the Forward and Reverse reads have been combined, and thus what you're searching for needs to change. The sequence comes in as a string, I make a Seq object out of it, which then has certain methods available to it as highlighted below.
    
    
    def reverse_complement_reverse_primer(self, sequence):
        # Sequences are reverse complemented when F and R reads were merged. So same has to happen to primers.

        
        rc_sequence = Seq(sequence, IUPAC.ambiguous_dna)
        
        return rc_sequence.reverse_complement()
        
    
    # This function takes the sample_id with primers and the regex primer patterns and reverses the 
    
    def get_sample_id_from_primer_sequence(self, sample_id_primer_dict, f_primer_pattern, r_primer_pattern):
        
        
        rectified_f_primer = self.replace_ambiguous_pattern_with_iupac_base(f_primer_pattern)
        
        rectified_r_primer = self.replace_ambiguous_pattern_with_iupac_base(r_primer_pattern)
        
        # reverse complement after replacing the bases, otherwise the brackets are changed as well
        rectified_r_primer = self.reverse_complement_reverse_primer(rectified_r_primer)
        
        #values_set = set()
        count=0
        
        for key, values in sample_id_primer_dict.items():
            count+=1
            if count > len(sample_id_primer_dict):
                return "No_match"
            elif str(rectified_f_primer) in str(values[0]) and str(rectified_r_primer) in str(values[1]):
                return key
            else:
                continue

    
    def return_fastq_seqio_object(self, seq,label,qual,filename):
        
        '''Construct Biopython SeqIO object for each of the fastq reads'''
        

    
        fastq_string = ("@{0}\n{1}\n+{0}\n{2}\n".format(label+":sample_id_" + filename[:-6], seq, qual))
        try:
            record = SeqIO.read(StringIO(fastq_string), "fastq")
            return record
        except ValueError as self.e:
            #logger.info(e.args, label)
            return None
    
def test_function():

    opts = {'m' : './SI_LI_subsamples/Small_intestine_primers_and_pairs.tsv',
            'f' : './SI_LI_subsamples/SI_sub.fq',
            'o' : './output.log'}
    Demultiplex(opts)
    
    
if __name__ == '__main__':
    #test_function()
    # I don't actually know how the utilities below from QIIME work
# from qiime.util import parse_command_line_parameters, make_option, gzip_open
# Would probably be useful to add 

# 1. process the incoming arguments from the command line
# 2. check presence of the files that have been input in the 
# 3. DONE check validity of the metadata file by calling the QIIME library
#     
    parser = argparse.ArgumentParser(description='A primer demultiplexer for Illumina NGS samples when the  \
    primer sequences are not present in the headers and the files are not demultiplexed on the machine \
    it was run on.')

    parser.add_argument('-m', metavar='metadata file', required=True, type=argparse.FileType('r'))
    parser.add_argument('-f', metavar='sequence file - fastq only, not gzipped', required=True, type=argparse.FileType('r'))

    try:
        results = parser.parse_args()

        if not (results.m or results.f or results.o):
            parser.error("You have to specify -m and -f files!")
        else:
            #meta_present = check_if_files_exist(results.m)
            #ngs_present = check_if_files_exist(results.f)
            #output_present = check_if_files_exist(results.o)
            #if (meta_present, ngs_present, output_present) == True:
            start_run = Demultiplex()
            start_run.run_demultiplex_and_trim(results)
            #return results
    except IOError as e:
        parser.error(str(e))

    #opts = parse_command_line()
    #check_metadata_file_formatting(opts.m)
    #run_demultiplex_and_trim(opts)

#    if argv[2].endswith('.gz') and argv[3].endswith('.gz'):
#        f = gzip_open(argv[2])
#        r = gzip_open(argv[3])
#    else:
#        f = open(argv[2],'U')
#        r = open(argv[3],'U')
#
#    outf = open(argv[2], "w")

    
