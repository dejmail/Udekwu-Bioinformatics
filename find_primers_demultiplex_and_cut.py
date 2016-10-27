import os
import sys

import logging
import logging.config
import yaml

from string import upper
from re import match
from re import compile
import argparse
import gzip
import random, string

from skbio.sequence import DNA
 
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Alphabet import IUPAC, generic_dna
        
from qiime.check_id_map import process_id_map
import subprocess

import progressbar

class Demultiplex(object):
    
    def __init__(self, opts, logger=None):
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
        
        self.record_buffer = {}
        self.iupac = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', \
                      'R': '[AG]', 'Y': '[CT]','S': '[GC]',\
                      'W': '[AT]', 'K': '[GT]', 'M': '[AC]',\
                      'B': '[CGT]','D': '[AGT]', 'H': '[ACT]',\
                      'V': '[ACG]', 'N': '[ACGT]'}
        self.search_length = int(opts.l)
        self.r1_filename = opts.f
        self.r2_filename = opts.r
        self.inverted_iupac = self.invert_regex_pattern(self.iupac)
        self.R1, self.r1_tot, self.R2, self.r2_tot = self.index_sequences([self.r1_filename, self.r2_filename])
            
    def number_of_lines_in_file(self, filename):
        
        import logging
        self.logger = logging.getLogger('_cntlne_')
            
        with open(filename) as f:
            for i, l in enumerate(f):
                pass
        self.logger.debug("{0} has {1} lines".format(filename, str(i+1)))
        return i + 1
        
        
    def index_sequences(self, filenames):
        
        import logging
        self.logger = logging.getLogger('_idxseq_')
        
        if ("R1" in filenames[0].name) and ("R2" in filenames[1].name):
            
            self.r1_number_of_lines = self.number_of_lines_in_file(filenames[0].name)
            self.logger.debug("thus {0} reads".format(self.r1_number_of_lines/4))
            self.r1_sequences = SeqIO.index(filenames[0].name, 'fastq', generic_dna)
            #self.logger.debug(self.r1_sequences)
            
            self.r2_number_of_lines = self.number_of_lines_in_file(filenames[1].name)
            self.logger.debug("thus {0} reads".format(self.r2_number_of_lines/4))
            self.r2_sequences = SeqIO.index(filenames[1].name, 'fastq', generic_dna)
            #self.logger.debug(self.r2_sequences)
            
            return self.r1_sequences, self.r1_number_of_lines,\
            self.r2_sequences, self.r2_number_of_lines
        else: 
            raise IOError("Can't find R1 or R2 in filenames")
            return None    
                
                
    def get_primers(self, header, mapping_data):
        
        """ Returns lists of forward/reverse primer regular expression
    
            header:  list of strings of header data.
            mapping_data:  list of lists of mapping data
    
            Will raise error if either the LinkerPrimerSequence or ReversePrimer fields
            are not present
        """
        import logging
        self.logger = logging.getLogger('_getprm_')
        
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
            # reverse primer were reverse complemented
            raw_reverse_primers.update([upper(str(DNA(primer))) for
                                            primer in line[rev_primer_ix].split(',')])
    
        if not raw_forward_primers:
            self.logger.critical("No forward primers detected in mapping file.")
            raise ValueError(("No forward primers detected in mapping file."))
            
        if not raw_reverse_primers:
            self.logger.critical("No reverse primers detected in mapping file.")
            raise ValueError(("No reverse primers detected in mapping file."))

     
        forward_primers = []
        forward_primers_rc = []
        reverse_primers = []
        reverse_primers_rc = []

        for curr_primer in raw_forward_primers:
    
            forward_primers.append(compile(''.join([self.iupac[symbol] for symbol in curr_primer[:self.search_length]])))
            forward_primers_rc.append(compile(''.join([self.iupac[symbol] for symbol in self.reverse_complement(curr_primer[:self.search_length])])))
        
        for curr_primer in raw_reverse_primers:
            reverse_primers.append(compile(''.join([self.iupac[symbol] for symbol in curr_primer[:self.search_length]])))
            reverse_primers_rc.append(compile(''.join([self.iupac[symbol] for symbol in self.reverse_complement(curr_primer[:self.search_length])])))
            
        return forward_primers, forward_primers_rc, reverse_primers, reverse_primers_rc
            
    def reverse_complement(self, sequence):
    # Sequences are reverse complemented when F and R reads were merged. So same has to happen to primers.
    
        rc_sequence = Seq(sequence, IUPAC.ambiguous_dna)
    
        return rc_sequence.reverse_complement()
    
    def write_fastq_sequence(self, fastq_list, filename):
        
        """ Takes in a filename and fastq entry as text and constructs
            a SeqIO object that is written to file. This helps ensure 
            that the output is correctly formatted when written to file.        
        """
        import logging
        self.logger = logging.getLogger('_fqwrte_')
        
        if len(filename) < 12:
            try:    
                r1_handle = filename + "_R1.fq"
                with open(r1_handle, 'a', 3290000) as self.out_seqs:
                    SeqIO.write(fastq_list[0], self.out_seqs, "fastq")
                    
                r2_handle = filename + "_R2.fq"
                with open(r2_handle, 'a', 3290000) as self.out_seqs:
                    SeqIO.write(fastq_list[1], self.out_seqs, "fastq")
                    
            except TypeError as e:
                self.quality_errors += 1
                self.logger.debug(fastq_list, e.args, e.message)
        else:
            with open(filename, 'a') as self.out_seqs:
                    SeqIO.write(fastq_list, self.out_seqs, "fastq")
#    def process_incoming_reads(self, sequence_file, orientation):
#        
#        for record in SeqIO.parse(handle=sequence_file,format="fastq-sanger"):
#            label = record.id
#            seq = record.seq
#            seq_rc = record.seq.reverse_complement()
#            qual = record.letter_annotations["phred_quality"]
#        
#            return None

    def determine_sequence_slices(self, sample_id, sample_primer_dict, search_result):
        
        try:
            r1_start_slice = search_result[0].get('start_position')
            r1_end_slice = len(sample_primer_dict.get(sample_id[0])[0])
    
            r2_start_slice = search_result[1].get('start_position')
            r2_end_slice = len(sample_primer_dict.get(sample_id[0])[1])
            
            return {'r1' : [r1_start_slice, r1_end_slice]}, {'r2' : [r2_start_slice, r2_end_slice]}
        
        except KeyError as e:
            print(e)
            return "Error with slice determination"        

    def regex_search_through_sequence(self, search_dict, regex_dict):
        
        ''' Takes in a regex pattern dictionary-list, and using the list values
        items searches through a tuple of strings until it either finds a match or not.
        
        Returns True/False (default False), the start and end position on the string and the
        regex pattern found and the tuple index indicating + or - orientation, and the 
        orientation of the search pattern used '''
        
        import logging
        self.logger = logging.getLogger('_regex__')
        
        pair_result = []
        for strand_key, sequence in search_dict.iteritems():
            for orient_key, pattern_list in regex_dict.items():
                for search_pattern in pattern_list:                
                    search_match = match(search_pattern, str(sequence.seq))
                    if search_match:
                        pair_result.append({'pattern_found' : True, 
                                            'pattern' : search_pattern.pattern,
                                            'start_position' : search_match.span()[0],
                                            'end_position' : search_match.span()[1],
                                            'index' : strand_key,
                                            'orient_key' : orient_key})
                        break
                    else:
                        continue
                    break
        else:
            try:
                if (pair_result == None) or (len(max(pair_result, key=len)) < 6):
                    pair_result.append({'pattern_found' : False,
                                        'index' : strand_key})    
            except ValueError:
                pair_result.append({'pattern_found' : False,
                                        'index' : strand_key})

        return pair_result

    def clip_primers_from_seq(self, search_result, primer_dict, seq_dict, sample_primer_dict, sample_id):
    
        # search_result - list of dictionaries
        # primer_dict - dict of lists
        # seq_dict - dictionary of seqio objects
        
        # loop through each sequence, whilst looping through each pattern_search
        # result. Match the keys and use the search result to clip the sequence
        # return a 
        
        # | search_result
        # 'pattern_found' : True, 
        # 'pattern' : search_pattern.pattern,
        # 'start_position' : search_match.span()[0],
        # 'end_position' : search_match.span()[1],
        # 'index' : strand_key,
        # 'orient_key' : orient_key})
        
        # | seq_dict
        # {'r1' : r1 sequence info, 'r2' : r2 sequence info}
        
        #(search_result, primer_dict, seq_dict, sample_primer_dict, sample_id)
        
        new_record = {}
        
        r1_slices, r2_slices = self.determine_sequence_slices(sample_id, sample_primer_dict, search_result)
        
        record_r1 = search_result[0]
        record_r2 = search_result[1]
            
        for seq_key, record in seq_dict.iteritems():
            
            if seq_key in record_r1.values():
                
                tmp_record = SeqRecord.SeqRecord(id=record.id,
                            seq=record.seq[r1_slices.get('r1')[0]:r1_slices.get('r1')[1]],
                            description=sample_id[0],
                            letter_annotations={'phred_quality' : record.letter_annotations['phred_quality'][r1_slices.get('r1')[0]:r1_slices.get('r1')[1]]})
                
                new_record.setdefault('r1', {}).append(tmp_record)
                del(tmp_record)
                
            elif seq_key in record_r2.values():
                
                tmp_record = SeqRecord.SeqRecord(id=record.id,
                            seq=record.seq[r2_slices.get('r2')[0]:r2_slices.get('r2')[1]],
                            description=sample_id[0],
                            letter_annotations={'phred_quality' : record.letter_annotations['phred_quality'][r2_slices.get('r2')[0]:r2_slices.get('r2')[1]]})
        
                new_record.setdefault('r2', {}).append(tmp_record)
                del(tmp_record)
                    
        return new_record
    
    def record_buffer_and_writer(self, record_dict):
        values = []
        
        for key, value in record_dict.items():
            self.record_buffer.setdefault(key, {}).append(value)
        
        
        for entry in self.record_buffer.values(): values.append(len(entry))
        
        print(sum(values))
        if (sum(values) >= 25) or (100 - sum(values) <= 25):
            print(sum(values))
            #read_key = r1, r2, discarded
            for pair_key, record in self.record_buffer.items():
                print(pair_key)
                #filename = record.description + "_" + read_key + ".fq"            
                #print(filename
            values=[]
            return "Buffer cleared"
        
        return "Buffer < 25"
            
    def run_demultiplex_and_trim(self, opts, **kwargs):
        
        """
            The main part of the script that pulls all the various 
            manipulations together. It takes arguments from the command
            line as well as **kwargs (currently only specifying gzip or not)
        
        """
        
        import logging
        self.logger = logging.getLogger('demultip')
        
        sample_primer_dict = {}
        file_list = []

        
        if not opts:
            sys.exit("command line options not getting to main method")
        

        metafile = opts.m
        output_directory = opts.o
        filename=""
        #check_both_orientations = opts.reverse_complement
        
        # extract .gz to temp file location
        if 'gzipFilename' in kwargs:
            self.logger.info("Incoming kwargs detected...gzip file?")
            #sequence_file = kwargs.get('gzipFilename')
        else:
            self.logger.info("No kwargs, normal Fastq file")
            #sequence_file = opts.f
        self.logger.info("processing {0} total sequences".format(str((self.r1_tot+self.r2_tot)/4)))
        self.logger.info("using the first {0} bases of primer in search".format(self.search_length))
        #self.logger.info("opt.f.name - {0}".format(sequence_file))
        #self.logger.info("sequences being read from {0}".format(sequence_file))
        #self.logger.info("checking in both sequence orientations = {0}".format(check_both_orientations))
    
        #extract the relevant data from the metadata file
        header, mapping_data, run_description, errors, warnings = process_id_map(metafile)
        self.logger.debug("csv mapping data from {0}...\n{1}".format(metafile, "\n".join([str(x) for x in mapping_data])))
        
        # get the primer search patterns
        #forward_primer_patterns, reverse_primer_patterns = self.get_primers(header, mapping_data)
        forward_primers, forward_primers_rc, reverse_primers, reverse_primers_rc = self.get_primers(header, mapping_data)
        self.primer_pattern_dict_list = {'fp' : forward_primers, 'fprc' : forward_primers_rc, 'rp' : reverse_primers, 'rprc' : reverse_primers_rc}
        
        
        self.logger.debug("forward_primer patterns\n{0}\n".format("\n".join([str(x.pattern) for x in self.primer_pattern_dict_list.get('fp')])))
        self.logger.debug("reverse_primers patterns\n{0}\n".format("\n".join([str(x.pattern) for x in self.primer_pattern_dict_list.get('rp')])))
        
        # replace colons with underscores in the sample_id names
        for samples in mapping_data:
            sample_primer_dict[samples[0].replace(":","_")] = (samples[2], samples[5])
            file_list.append(samples[0].replace(":","_"))
            
        self.logger.debug("sample_primer_dict...{0}".format("\n".join(x) for x in sample_primer_dict.items()))
        
        bar = progressbar.ProgressBar(max_value=(self.r1_tot+self.r2_tot)/4,redirect_stdout=True)
    
        for r1, r2 in zip(self.R1.values(), self.R2.values()):
            
            #self.pair_read = {'r1': r1, 'r2': r2}
            #self.logger.debug("Starting search through R1/R2 reads")
            #self.label = record.id
            #self.seq = record.seq
            #self.seq_rc = record.seq.reverse_complement()
            #self.qual = record.letter_annotations["phred_quality"]
            seq_dict = {'r1' : r1, 'r2' : r2}

            bar.update(self.total_seqs)        
            
            self.logger.debug("\nprocessing seq ID - R1 {0}... R2 {1}".format(r1.id, r2.id))
            self.logger.debug("R1 sequence - {0}".format(r1.seq))
            self.logger.debug("R2 sequence - {0}".format(r2.seq))

            self.sample_id = ""
            # because we process two sequences at a time (R1 and R2)
            self.total_seqs += 2
            
            # tuple with both orientations
            #self.curr_r1 = (r1.seq, r1.seq.reverse_complement())
            #self.curr_r2 = (r2.seq, r2.seq.reverse_complement())
            self.pair_read = {'r1' : r1, 'r2' :r2}
            #start_slice = 0
            #end_slice = -1
            self.f_primer_found = []
            self.r_primer_found = []
            
            self.logger.debug("Looking in pair read for patterns...")
            
            search_result = self.regex_search_through_sequence(seq_dict, self.primer_pattern_dict_list)
            #print(search_result)
            try:
                sample_id = self.get_sample_id_from_primer_sequence(sample_primer_dict, search_result[0].get('pattern'), search_result[1].get('pattern'))
            except IndexError as e:
                # sample is missing one or both the patterns keys
                print("Sample seq is missing a a pattern, {0} - {1}".format(sample_id, e))
                output = self.record_buffer_and_writer({'discarded' : seq_dict})
                del(output)
                continue
            try:
                new_seq = self.clip_primers_from_seq(search_result, self.primer_pattern_dict_list, seq_dict, sample_primer_dict, sample_id)
            except Exception as e:    
                print("attempt to get new sequence failed - e - {0}".format(e))
                output = self.record_buffer_and_writer({'discarded' : seq_dict})
                del(output)
                continue
        
            output = self.record_buffer_and_writer(new_seq)
            
            #primer_pair_search = self.regex_search_through_sequence(self.pair_read,  
            #                                                          self.primer_pattern_dict_list)
            #r2_primer_search = self.regex_search_through_sequence(str(self.curr_r1), 
            #                                                           self.primer_pattern_dict_list)
            #self.logger.debug("r1_forward_primer_search {0}...r1_reverse_primer_search {1}".format(r1_primer_search, r2_primer_search))
            
#            self.logger.debug("Looking in R2 for patterns...")
#            r2_forward_primer_search = self.regex_search_through_sequence(str(self.curr_r2), 
#                                                                       forward_primer_patterns)
#            r2_reverse_primer_search = self.regex_search_through_sequence(str(self.curr_r2),
#                                                                       reverse_primer_patterns)
#            self.logger.debug("r2_forward_primer_search {0}... r2_reverse_primer_search {1}".format(r2_forward_primer_search, r2_reverse_primer_search))
#            
#            self.logger.debug("regex result lengths -> {0}-{1}-{2}-{3}".format(len(r1_forward_primer_search),len(r1_reverse_primer_search),len(r2_forward_primer_search),len(r2_reverse_primer_search)))
#            #if len(r1_forward_primer_search) == 0: raise Exception ("Where is the regex result for self_f_primer_found ?")
            #if len(r1_reverse_primer_search) == 0: raise Exception ("Where is the regex result for self_r_primer_found ?")
 
#            
#            if primer_pair_search[0].get('pattern_found') == True:
#                self.logger.debug("F primer found between {0}-{1}".format(primer_pair_search[0].get('start'),primer_pair_search[0].get('end') ))
#                self.f_primer_found = primer_pair_search[0].get('pattern_found')
#            #    self.f_primer_found.append('r1f')
#            #    self.f_primer_found.append(r1_forward_primer_search)
#            
#            elif primer_pair_search[1].get('pattern_found') == True:
#                self.logger.debug("R primer found between {0}-{1}".format(primer_pair_search[1].get('start'),primer_pair_search[1].get('end')))
#                self.r_primer_found = primer_pair_search[1].get('pattern_found')
#                #self.f_primer_found.append('r2f')
#                #self.f_primer_found.append(r2_forward_primer_search)
#            else:
#                self.logger.debug("No search contig found for this sequence")
#                
                #self.f_primer_found = ['', {'pattern_found' : False}]
                
#            if r1_reverse_primer_search.get('pattern_found') == True:
#                self.logger.debug("found at tuple position {0}".format(r1_reverse_primer_search.get('index')))
#                self.r_primer_found.append('r1r')
#                self.r_primer_found.append(r1_reverse_primer_search)
#            elif r2_reverse_primer_search.get('pattern_found') == True:
#                self.logger.debug("found at tuple position {0}".format(r2_reverse_primer_search.get('index')))
#                self.r_primer_found.append('r2r')
#                self.r_primer_found.append(r2_reverse_primer_search)
#            else:
#                self.r_primer_found = ['', {'pattern_found' : False}]
#            
#            # can change this to a try catch
#            if len(self.f_primer_found) > 1 and len(self.r_primer_found) > 1:
#                self.logger.debug('F primer-{0}...R primer-{1}'.\
#                                format(self.f_primer_found[1].get('pattern_found'),
#                                       self.r_primer_found[1].get('pattern_found')))
#            
#            if (check_both_orientations == 'True') and not (f_primer_found == True) and not (r_primer_found == True):                    
#                #raw_input()
#                self.logger.debug('Looking for F primer in reverse complement')
#                forward_primer_search = self.regex_search_through_sequence(str(seq_rc), forward_primer_patterns)
#                
#                self.logger.debug('Looking for R primer in reverse complement')
#                reverse_primer_search = self.regex_search_through_sequence(str(seq_rc), 
#                                                                       reverse_primer_patterns)
#                f_primer_found = forward_primer_search.get('pattern_found')
#                r_primer_found = reverse_primer_search.get('pattern_found')
            #forward_primer_search= {'pattern_found' : False}
            
        
#            for curr_primer in forward_primers:
#                # regex search using curr_primer pattern
#                if curr_primer.search(str(seq)):
#                    self.logger.debug("found F primer pattern...{0}".format(curr_primer.pattern))        
#                    f_primer_seq = curr_primer
#                    start_slice = int(curr_primer.search(str(seq)).span()[1])
#                    self.f_count += 1
#                    f_primer_found = True


#            for curr_primer in reverse_primers:
#                if curr_primer.search(str(seq)):
#                    self.logger.debug("found R primer pattern...{0}".format(curr_primer.pattern))
#                    r_primer_seq = curr_primer
#                    # span() gives start and end coordinates
#                    end_slice = int(curr_primer.search(str(seq)).span()[0])
#                    self.r_count += 1
#                    r_primer_found = True
#            
        #  created a clipped sequence and quality score
        #curr_seq = seq[start_slice:end_slice]
        #curr_qual = qual[start_slice:end_slice]
        
        # can change this value so it can be specified from cmdline
        #if /len(curr_seq) < 1:
        #    self.no_seq_left += 1
        #    continue
            #self.logger.debug("self.f_primer_found {0}\nself.r_primer_found {1}".format(self.f_primer_found, self.r_primer_found))
#            if (self.f_primer_found == True) and (self.r_primer_found == True):
#        
#                #self.logger.debug("Both F and R primers found - \n R1 sequence {0}\n R2 sequences {1} \n forward pattern - {2} \n start_position - {3} \n end_position {4}\n reverse pattern - {3} \n reverse pattern - {4}\n".format(r1.seq, r2.seq, self.f_primer_found[1].get('pattern_found'), str(self.f_primer_found[1].get('start_position')),str(self.r_primer_found[1].get('end_position')),self.r_primer_found[1].get('pattern')))
#            #else: self.logger.debug("No primer match found, continuing...{0}".format(e))            
#            
#            # get filename from reversing the regex pattern
#            #self.logger.debug("sample_primer_dict - {0}".format(sample_primer_dict))
#
#            
#            #if len(self.f_primer_found) > 0 and len(self.r_primer_found) > 0:
#                #self.logger.debug("self.f_primer_found {0}... self.r_primer_found {1}".format(self.f_primer_found, 
#                #                                                                               self.r_primer_found))
#                #self.logger.debug("self.f_primer_found length = {0}, self.r_primer_found length {1}".format(str(len(self.f_primer_found)), 
#                 #                                                                                           str(len(self.r_primer_found))))
#                
#                self.sample_id = self.get_sample_id_from_primer_sequence(sample_primer_dict, self.primer_pair_search[0].get('pattern'), self.primer_pair_search[1].get('pattern'))
#                #filename = self.generate_filename(self.f_primer_found, self.r_primer_found, r1, r2, )
#            
#                self.logger.debug("assigning read {0} to sample id - {1}".format(r1.id, self.sample_id))
#
#                if self.sample_id == None:
#                    self.unmapped_count += 2
#                    
#                    self.logger.debug("Failed to get filename from primer sheet"+
#                                  "{0}\n{1}\n{2}\n{3}\n".format(r1.id, r2.id, 
#                                  self.f_primer_found[1].get('pattern'),
#                                  self.r_primer_found[1].get('pattern')))
#                    incorrect_primer_pairs_filename = os.path.join(output_directory,
#                                                                   "incorrect_primer_pairing.fq")
#                    
#                    with open(incorrect_primer_pairs_filename, 'a') as f:
#                        SeqIO.write([r1,r2], f, "fastq")
#                        
#                else:
#                    self.both_primers_count += 2
#                    
#                    # record is returned as a list
#                    record = self.return_fastq_seqio_object([r1, r2], filename)
#                    self.logger.debug("sample_id is...{0}".format(self.sample_id))
#                    self.write_fastq_sequence(record, os.path.join(output_directory, self.sample_id))  
#            
#            if (self.f_primer_found[1].get('pattern_found') == True) and\
#                (self.r_primer_found[1].get('pattern_found') == False):
#                self.logger.debug("Only forward primer found")
#                self.f_only_count += 1
#                fastq_entry = self.return_fastq_seqio_object([r1, r2], "Only_F")
#                self.write_fastq_sequence(fastq_entry, os.path.join(output_directory, "Only_F.fq"))
#                
#            if (self.f_primer_found[1].get('pattern_found') == False) and\
#                (self.r_primer_found[1].get('pattern_found') == True):
#                self.logger.debug("Only reverse primer found")
#                self.r_only_count += 1
#                fastq_entry = self.return_fastq_seqio_object([r1, r2],"Only_R")
#                self.write_fastq_sequence(fastq_entry, os.path.join(output_directory, "Only_R.fq"))
#            
#            if (self.f_primer_found[1].get('pattern_found') == False) and\
#                (self.f_primer_found[1].get('pattern_found') == False):
#                self.logger.debug("Neither F nor R primer found")
#                self.unmapped_count += 1
#                fastq_entry = self.return_fastq_seqio_object([r1, r2], "None")
#                self.write_fastq_sequence(fastq_entry, os.path.join(output_directory, "None.fq"))                
                
                    
        self.logger.info("__________________________")
        self.logger.info("Samples successfully mapped (F+R found): {0}".format(self.both_primers_count))        
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
            the sequence header. The pattern may also have ambiguities from wobble
            bases.
            
            This is accomplished by making a new IUPAC dictionary and simply 
            swapping the keys and values of the original iupac dictionary. 
        '''
        import logging
        self.logger = logging.getLogger('chgambig')
        
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
            against the list of primers provided by the metadata file. Using set
            means that only one sample_id is returned if both match as sets can
            only contain unique items, otherwise there are two items in the set.
        '''
        import logging
        self.logger = logging.getLogger('smplf4id')
        
        self.logger.debug("getting sample ID from the primer sequences")
        self.logger.debug("looking for patterns F - {0} and R - {1}".format(f_primer_pattern, r_primer_pattern))
        
        rectified_f_primer = self.replace_ambiguous_pattern_with_iupac_base(f_primer_pattern)
        
        rectified_r_primer = self.replace_ambiguous_pattern_with_iupac_base(r_primer_pattern)
        self.logger.debug("primer converted to F - {0} and R - {1}".format(rectified_f_primer, rectified_r_primer))
        
        values_set = set()
        
        for key, values in sample_id_primer_dict.items():
            
            if any(f_primer_pattern in f for f in values):
                values_set.add(key)
            elif any(r_primer_pattern in r for r in values):
                values_set.add(key)
                
            if len(values_set) == 2:
                return values_set
            else: 
                continue
                
        return list(values_set)
        

    def return_fastq_seqio_object(self, data, filename):
    
        '''Construct Biopython SeqIO object for each of the fastq reads
            in the incoming list '''
            
        import logging
        self.logger = logging.getLogger('fqseqIOr')
        
        record_list = []
        self.logger.debug("data coming in return_seqio object {0}".format(data))
        for items in data:
            
            record = SeqRecord.SeqRecord(id=items.id + 
                                        ":sample_id_" + 
                                        filename,
                                        seq=items.seq, 
                                        letter_annotations={'phred_quality' : items.letter_annotations['phred_quality']})
            record_list.append(record)
        
        return record_list
        
    def setup_logging(self, default_path='logging.yaml',
                  default_level=logging.DEBUG,
                  env_key='LOG_CFG'):
    
        """Setup logging configuration"""
    
        from yaml import safe_load
        
        path = default_path
        value = os.getenv(env_key, None)
        if value:
            path = value
        if os.path.exists(path):
            with open(path, 'rt') as f:
                config = safe_load(f.read())
            logging.config.dictConfig(config)
            
        else:
            logging.basicConfig(level=default_level)        
    


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
    import logging.config
    parser = argparse.ArgumentParser(description='A primer demultiplexer for Illumina NGS samples when the  \
    primer sequences are not present in the sequence header and the files were not demultiplexed on the machine \
    it was run on. It takes a metadata file (QIIME specific), a merged fastq file')

    parser.add_argument('--m', metavar='--> metadata file, QIIME formatted', required=True, type=argparse.FileType('r'))
    parser.add_argument('--f', metavar='--> sequence file - fastq or gzipped', required=True, type=argparse.FileType('r'))
    parser.add_argument('--r', metavar='--> reverse sequence file', required=True, type=argparse.FileType('r'))
    parser.add_argument('--l', metavar='--> length of primer to use in search', required=False, action='store', help='Length of primer in bp from 5 prime end for search')
    parser.add_argument('--t', metavar='--> trim barcode-primer from sequence', required=True, action='store', help='Whether to trim the sequencing used for demultiplexing or not')
    parser.add_argument('--o', metavar='--> output directory', required=True ,action='store')
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
                #print("R1 sequence file - {0}".format(results.f.name))
                #print("R2 sequence file - {0}".format(results.r.name))
                #print("input metadata file - {0}".format(results.m.name))
                start_run = Demultiplex(results)
                start_run.run_demultiplex_and_trim(results)
                
    except IOError as e:
        parser.error(str(e))
