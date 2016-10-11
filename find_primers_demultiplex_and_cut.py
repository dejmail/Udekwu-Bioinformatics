import os
import sys

import logging
import logging.config
import yaml

from string import upper
from re import compile
import argparse
import gzip
import random, string

from cogent.parse.fastq import MinimalFastqParser
from skbio.sequence import DNA
 
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
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
        self.iupac = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', \
                      'R': '[AG]', 'Y': '[CT]','S': '[GC]',\
                      'W': '[AT]', 'K': '[GT]', 'M': '[AC]',\
                      'B': '[CGT]','D': '[AGT]', 'H': '[ACT]',\
                      'V': '[ACG]', 'N': '[ACGT]'}
        self.inverted_iupac = self.invert_regex_pattern(self.iupac)
        self.open_file_list = set()
        
    def setup_logging(self,
                      default_path='logging.yaml',
                      default_level=logging.DEBUG,
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
            forward_primers.append(compile(''.join([self.iupac[symbol] for
                                                    symbol in curr_primer])))
        for curr_primer in raw_reverse_primers:
            reverse_primers.append(compile(''.join([self.iupac[symbol] for
                                                    symbol in curr_primer])))
         
        return forward_primers, reverse_primers
        
    def write_fastq_sequence(self, filename, fastq_entry):
        
        """ Takes in a filename and fastq entry as text and constructs
            a SeqIO object that is written to file. This helps ensure 
            that the output is correctly formatted when written to file.        
        """
        
        self.out_seqs = open(filename, 'a')                    
        try:
            #construct a SeqIO object, which checks for proper quality scores and formatting
            SeqIO.write(fastq_entry, self.out_seqs, "fastq")
            self.out_seqs.close()
        except TypeError as self.e:
            self.quality_errors += 1
            self.logger.debug(fastq_entry, e.args, e.message)
            
    def run_demultiplex_and_trim(self, opts, **kwargs):
        
        """
            The main part of the script that pulls all the various 
            manipulations together. It takes arguments from the command
            line as well as **kwargs (currently only specifying gzip or not)
        
        """
        sample_primer_dict = {}
        file_list = []

        
        if not opts:
            sys.exit("sequence file not getting to main method")
        self.logger.info("incoming cmdline opts - {0}".format(opts))
    
        self.logger.info("Forward primers found: {0}".format(self.f_count))
        self.logger.info("Reverse primers found: {0}".format(self.r_count))        
        self.logger.info("Samples successfully mapped F+R found): {0}".format(self.both_primers_count))        
        self.logger.info("Only forward primer found: {0}".format(self.f_only_count))
        self.logger.info("Only reverse primer found: {0}".format(self.r_only_count))
        self.logger.info("No seq left after truncation: {0}".format(self.no_seq_left))
        self.logger.info("Sequences with errors in quality scores: {0}".format(self.quality_errors))
        self.logger.info("Sequences not mapped: {0}".format(self.unmapped_count))
        self.logger.info("Total sequences checked: {0}".format(self.total_seqs))

        metafile = opts.m
        output_directory = opts.o
       
        # extract .gz to temp file location
        if 'gzipFilename' in kwargs:
            self.logger.info("Incoming kwargs detected...gzip file?")
            sequence_file = kwargs.get('gzipFilename')
        else:
            self.logger.info("No kwargs, normal Fastq file")
            sequence_file = opts.f
        
        self.logger.info("opt.f.name - {0}".format(sequence_file))
        self.logger.info("sequences being read from {0}".format(sequence_file.name))
    
        #extract the relevant data from the metadata file
        header, mapping_data, run_description, errors, warnings = process_id_map(metafile)
        self.logger.debug(mapping_data)
        
        # get the primer search patterns
        forward_primers, reverse_primers = self.get_primers(header, mapping_data)
        self.logger.debug(forward_primers)
        self.logger.debug(reverse_primers)
        
        # replace colons with underscores in the sample_id names
        for sample_id in mapping_data:
            sample_primer_dict[sample_id[0].replace(":","_")] = (sample_id[2], sample_id[4])
            file_list.append(sample_id[0].replace(":","_"))
            
        for label,seq,qual in MinimalFastqParser(sequence_file, strict=False):
            self.logger.debug(dir(seq))
            self.total_seqs += 1
            start_slice = 0
            end_slice = -1
            f_primer_found = False
            r_primer_found = False
            
            for curr_primer in forward_primers:
                # regex search using curr_primer pattern
                if curr_primer.search(seq):
                    self.logger.debug("found F primer pattern...{0}".format(curr_primer.pattern))        
                    f_primer_seq = curr_primer
                    start_slice = int(curr_primer.search(seq).span()[1])
                    self.f_count += 1
                    f_primer_found = True
            for curr_primer in reverse_primers:
                if curr_primer.search(seq):
                    self.logger.debug("found R primer pattern...{0}".format(curr_primer.pattern))
                    r_primer_seq = curr_primer
                    # span() gives start and end coordinates
                    end_slice = int(curr_primer.search(seq).span()[0])
                    self.r_count += 1
                    r_primer_found = True
            
            #  created a clipped sequence and quality score
            curr_seq = seq[start_slice:end_slice]
            curr_qual = qual[start_slice:end_slice]
            
            # can change this value so it can be specified from cmdline
            if len(curr_seq) < 1:
                self.no_seq_left += 1
                continue
    
            if (f_primer_found == True) and (r_primer_found == False):
                self.logger.debug("Only forward primer found")
                self.f_only_count += 1
                fastq_entry = self.return_fastq_seqio_object(seq,"OnlyF", qual,"Only_F")
                self.write_fastq_sequence(os.path.join(output_directory, "Only_F.fq"), fastq_entry )
                
            if (f_primer_found == False) and (r_primer_found == True):
                self.logger.debug("Only reverse primer found")
                self.r_only_count += 1
                fastq_entry = self.return_fastq_seqio_object(seq,"OnlyR", qual,"Only_R")
                self.write_fastq_sequence(os.path.join(output_directory, "Only_R.fq"), fastq_entry )
            
            if (f_primer_found == False and r_primer_found == False):
                self.logger.debug("Neither F nor R primer found")
                self.unmapped_count += 1
                fastq_entry = self.return_fastq_seqio_object(seq,"None", qual,"None")
                self.write_fastq_sequence(os.path.join(output_directory, "None.fq"), fastq_entry )                
            
            if (f_primer_found == True and r_primer_found == True):
                self.logger.debug("Both F and R primers found")
                self.logger.debug(curr_seq, "\n", f_primer_seq.pattern,
                                  start_slice, end_slice, r_primer_seq.pattern)
                self.both_primers_count += 1
                # get filename from reversing the regex pattern
                filename = self.get_sample_id_from_primer_sequence(sample_primer_dict, f_primer_seq.pattern, r_primer_seq.pattern)
                filename = str(filename) + ".fq"                
                filename = os.path.join(output_directory, filename)
                record = self.return_fastq_seqio_object(curr_seq,label,curr_qual,filename)                
                self.write_fastq_sequence(filename, record)
    
        self.logger.info("__________________________")
        self.logger.info("Forward primers found: {0}".format(self.f_count))
        self.logger.info("Reverse primers found: {0}".format(self.r_count))        
        self.logger.info("Samples successfully mapped F+R found): {0}".format(self.both_primers_count))        
        self.logger.info("Only forward primer found: {0}".format(self.f_only_count))
        self.logger.info("Only reverse primer found: {0}".format(self.r_only_count))
        self.logger.info("No seq left after truncation: {0}".format(self.no_seq_left))
        self.logger.info("Sequences with errors in quality scores: {0}".format(self.quality_errors))
        self.logger.info("Sequences not mapped: {0}".format(self.unmapped_count))
        self.logger.info("Total sequences checked: {0}".format(self.total_seqs))
    
        self.logger.info("Run finished")
     
    
    def check_metadata_file_formatting(self, metafile):
    
        command = ("validate_mapping_file.py -m {0}".format(metafile))
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell='TRUE')
        output, err = p.communicate()
        
        return (output,err)
    
    def process_metadata_file(self, metafile):
    
        try:
            output = self.check_metadata_file_formatting(str(metafile))
            self.logger.debug("processed ok - {0}".format(output[0]))
        except IOError as e:
            self.logger.critical("The metadata file does not exist or the path is wrong - error {0}".format(e)) 

    
    def invert_regex_pattern(self, sample_primer_dict):
        
        
        '''
            Invert the regex pattern by looping through the keys and values
            in the dictionary and then reversing them such that 
            key=value, and value=key. 
        '''
        
        inverted_dict = {}
        
        for keys, values in sample_primer_dict.items():
            # only keep the ambiguous combinations
            if len(values) > 1:
                inverted_dict.update({values: keys})
        return inverted_dict    
    
    def replace_ambiguous_pattern_with_iupac_base(self, sequence):
        
        '''    
            To be able to match the sample_id and the associated primer names and 
            sequences with the correct reads, one needs to invert the regex pattern
            that was used to find the primers used as they are not specified in
            the sequence header.
            
            This is accomplished by making a new IUPAC dictionary and simply 
            swapping the keys and values of the original iupac dictionary. 
        '''
        
        
        self.logger.debug("reversing pattern {0}".format(sequence))
        
        #self.logger.info("reversed {0} pattern".format(inverted_iupac))
        #self.sequence = SeqIO.read(StringIO(sequence), "fasta")
        for key,value in self.inverted_iupac.items():
            if key in sequence:
                sequence = sequence.replace(key, value, 1)
            else:
                # loop does not continue if this is not here, strange?
                continue
        return str(sequence)
    
    
    # Simply reverse complement the primer sequence to reflect the fact that the Forward and Reverse reads have been combined, and thus what you're searching for needs to change. The sequence comes in as a string, I make a Seq object out of it, which then has certain methods available to it as highlighted below.
    
    
    def reverse_complement_reverse_primer(self, sequence):
        '''
            Simply returns a reverse complement of the sequence that comes
            into the function.
            
            Used to find the correct primer sequence for the reverse primers
            as the R reads are reverse complemented when the F and R reads
            are merged.
        '''
        
        rc_sequence = Seq(sequence, IUPAC.ambiguous_dna)
        
        return rc_sequence.reverse_complement()
        
    
    # This function takes the sample_id with primers and the regex primer patterns and reverses the 
    
    def get_sample_id_from_primer_sequence(self, sample_id_primer_dict, f_primer_pattern, r_primer_pattern):
        
        '''
            The actual function responsible for generating the sample ID, being
            the filename into which successfully mapped sequences will be placed.
            
            Once reversed from the regex pattern, the F and R sequences are checked
            against the list of primers provided by the metadata file.
        '''
        
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
        
        
        fastq_string = ("@{0}\n{1}\n+{0}\n{2}\n".format(label+":sample_id_" + filename, seq, qual))
        try:
            record = SeqIO.read(StringIO(fastq_string), "fastq")
            return record
        except ValueError as self.e:
            self.logger.debug(e.args, label)
            return None   


def check_if_path_exists(path):
    
    '''
        Check if path exists. Not part of the main method.
    '''
    
    try: 
        print("- Checking the output path {0}".format(path))
        if os.path.exists(path):
            print("Path was found - {0}".format(path))
            return True
        else:
            print("{0} NOT found, creating output directory".format(path))
            try:
                os.makedirs(path)
                return True
            except OSError as e:
               return False
    except (OSError, IOError) as e:
        print("The filepath {0} does not exist, cannot be found, or cannot be read".format(e))

def random_char(length):
    
    '''
        Generate a string of random characters. Used for appending to the 
        temporary filename if the incoming sequence file is gz format.
    '''
                
    random_string = ""
    random_string = ''.join(random.choice(string.ascii_letters) for x in range(length))
    
    return random_string
        
def extract_file(gzip_path, CURRENT_PATH="./"):
    
    '''
        Extract any Gzip incoming file into a temporary file and return
        that temporary name to the main method.
    '''
    
    inF = gzip.open(gzip_path, 'rt')
    # uncompress the gzip_path INTO THE 's' variable
    s = inF.read()
    inF.close()

    # get gzip filename (without directories)
    gzip_fname = os.path.basename(gzip_path)
    # get original filename (remove 3 characters from the end: ".gz")
    fname = random_char(7) + gzip_fname[:-3]
    uncompressed_path = os.path.join(CURRENT_PATH, fname)

    # store uncompressed file data from 's' variable
    open(uncompressed_path, 'w').write(s)
    
    return fname
    
if __name__ == '__main__':
    
    '''
        The method that runs when the script is invoked from the commandline
    '''

    parser = argparse.ArgumentParser(description='A primer demultiplexer for Illumina NGS samples when the  \
    primer sequences are not present in the sequence header and the files were not demultiplexed on the machine \
    it was run on. It takes a metadata file (QIIME specific), a merged fastq file')

    parser.add_argument('-m', metavar='metadata file', required=True, type=argparse.FileType('r'))
    parser.add_argument('-f', metavar='sequence file - fastq only, OR gzipped', required=True, type=argparse.FileType('r'))
    parser.add_argument('-o', metavar='output directory', required=True ,action='store')
    try:
        results = parser.parse_args()

        if not (results.m or results.f or results.o):
            parser.error("You have to specify the -m and -f files, and -o output directory!")
        else:
            if check_if_path_exists(results.o) == False:
                sys.exit("Something went wrong creating the path.")
            if results.f.name[-3:] == '.gz':
                
                fname = extract_file(results.f.name)
                start_run = Demultiplex()
                start_run.run_demultiplex_and_trim(results, gzipFilename=fname)
                os.remove(fname)
            else:
                print("input sequence file - {0}".format(results.f.name))
                print("input metadata file - {0}".format(results.m.name))
                start_run = Demultiplex()
                start_run.run_demultiplex_and_trim(results)
    except IOError as e:
        parser.error(str(e))
