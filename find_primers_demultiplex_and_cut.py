import os
import sys

import logging
import logging.config
import yaml

from string import upper, maketrans
from re import match, search
from re import compile
import argparse
import gzip
import random, string
import itertools

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
        self.processed_seqs = 0
        self.alternate_orientation = []
        
        self.record_buffer = {}
        self.iupac = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', \
                      'R': '[AG]', 'Y': '[CT]','S': '[GC]',\
                      'W': '[AT]', 'K': '[GT]', 'M': '[AC]',\
                      'B': '[CGT]','D': '[AGT]', 'H': '[ACT]',\
                      'V': '[ACG]', 'N': '[ACGT]'}
        self.search_length = int(opts.l)
        self.r1_filename = opts.f
        self.r2_filename = opts.r
        self.output_dir = opts.o
        self.inverted_iupac = self.invert_regex_pattern(self.iupac)
        self.R1, self.r1_tot, self.R2, self.r2_tot = self.index_sequence_files([self.r1_filename, self.r2_filename])
        self.starting_total = (self.r1_tot + self.r2_tot)/4
        self.buffer_limit = 5000
    
    
    def number_of_lines_in_file(self, filename):
        
        import logging
        self.logger = logging.getLogger('_cntlne_')
            
        with open(filename) as f:
            for i, l in enumerate(f):
                pass
        self.logger.debug("{0} has {1} lines".format(filename, str(i+1)))
        return i + 1
        
    
    def index_sequence_files(self, filenames):
        
        import logging
        self.logger = logging.getLogger('_idxseq_')
        
        self.logger.info("Indexing sequences....")
        if ("R1" in filenames[0].name.upper()) and ("R2" in filenames[1].name.upper()):
            
            self.r1_number_of_lines = self.number_of_lines_in_file(filenames[0].name)
            self.r1_sequences = SeqIO.index(filenames[0].name, 'fastq', generic_dna)
        
            self.r2_number_of_lines = self.number_of_lines_in_file(filenames[1].name)
            self.r2_sequences = SeqIO.index(filenames[1].name, 'fastq', generic_dna)
            
            return self.r1_sequences, self.r1_number_of_lines,\
            self.r2_sequences, self.r2_number_of_lines
        else: 
            raise IOError("Can't find R1 or R2 in filenames")
            return None    
                
    
    def create_primer_regex_patterns(self, header, mapping_data):
        
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
            raise ValueError("No forward primers detected in mapping file.")
            
        if not raw_reverse_primers:
            self.logger.critical("No reverse primers detected in mapping file.")
            raise ValueError("No reverse primers detected in mapping file.")

     
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
    
        rc_sequence = Seq(sequence, IUPAC.ambiguous_dna)
    
        return rc_sequence.reverse_complement()
    
    
    def determine_sequence_slices(self, sample_id, sample_primer_dict, search_result):
        
        '''Determine the coordinates where each sequence needs to be clipped 
           in order to remove the primer-barcode sequences
           
           sample_id is a tuple, so sampl_id[0] makes sure just the value is 
           returned.
           
           command like sample_primer_dict.get(sample_id[0])[0]) will for
           example get the first primer sequence from sample 1_2.        
           
        '''
        
        import logging
        self.logger = logging.getLogger('_slicer_')
        
        r1_seq_length = search_result[0].get('length')
        r2_seq_length = search_result[0].get('length')
        
        self.logger.debug("sample_id - {0}".format(sample_id))
        
        try:
            if len(search_result) == 4:
                
                self.logger.debug("pattern search contains {0} reads per primer pair".format(len(search_result)))
                r1_sublist = [num for num in search_result if num.get('index') == 'r1']
                for matches in r1_sublist:
                    if matches.get('start_position') < 40:
                        r1_start_slice = matches.get('start_position') + len(sample_primer_dict.get(sample_id[0])[0])
                
                    if matches.get('start_position') > r1_seq_length - 40:
                        r1_end_slice = r1_seq_length - len(sample_primer_dict.get(sample_id[0])[1])
                        
                        self.logger.debug("length of forward primer {0}".format(len(sample_primer_dict.get(sample_id[0])[0])))
                                    
                r2_sublist = [num for num in search_result if num.get('index') == 'r2']
                for matches in r2_sublist:
                
                    if matches.get('start_position') < 40:
                        r2_start_slice = matches.get('start_position') + len(sample_primer_dict.get(sample_id[0])[1])
            
                    if matches.get('start_position') > r1_seq_length-40:
                        r2_end_slice = r1_seq_length - len(sample_primer_dict.get(sample_id[0])[1])
                        self.logger.debug("length of reverse primer {0}".format(len(sample_primer_dict.get(sample_id[0])[1])))
    
                self.logger.debug("r1_start_slice - {0}, r1_end_slice {1}".format(r1_start_slice, r1_end_slice))
                self.logger.debug("r2_start_slice - {0}, r2_end_slice {1}".format(r2_start_slice, r2_end_slice))                
                self.logger.debug("sample_id[0][0] - {0}".format(sample_primer_dict.get(sample_id[0])[0]))
                self.logger.debug("sample_id[0][1] - {0}".format(sample_primer_dict.get(sample_id[0])[1]))
                
                return {'r1' : [r1_start_slice, r1_end_slice]}, {'r2' : [r2_start_slice, r2_end_slice]}
            elif len(search_result) == 3:
                self.logger.debug("pattern search contains {0} positive patterns - {1}".format(len(search_result), search_result))
                return None
                
            elif len(search_result) == 2:
                self.logger.debug("reads have one primer match each")
                r1_start_slice = search_result[0].get('start_position')
                r1_end_slice = len(sample_primer_dict.get(sample_id[0])[0])
                self.logger.debug("R1 start slice {0} - end slice {1}".format(r1_start_slice, r1_end_slice))
                
                r2_start_slice = search_result[1].get('start_position')
                r2_end_slice = len(sample_primer_dict.get(sample_id[0])[1])
                self.logger.debug("R2 start slice {0} - end slice {1}".format(r2_start_slice, r2_end_slice))
            
                self.logger.debug("sample_id[0][0] - {0}".format(sample_primer_dict.get(sample_id[0])[0]))
                self.logger.debug("sample_id[0][1] - {0}".format(sample_primer_dict.get(sample_id[0])[1]))
                
            return {'r1' : [r1_end_slice, r1_seq_length]}, {'r2' : [r2_end_slice, r2_seq_length]}
        
        except KeyError as e:
            self.logger.error("Error with slice determination...{0}...{1}".format(search_result, e))
            return "Error with slice determination"
    
    
    def regex_search_through_sequence(self, search_dict, regex_dict):
        
        ''' Takes in a regex pattern dictionary-list, and using the list values
        items searches through a tuple of strings until it either finds a match or not.
        
        Returns True/False (default False), the start and end position on the string and the
        regex pattern found and the tuple index indicating + or - orientation, and the 
        orientation of the search pattern used 
        
        Extra hits are removed if they are outside of the search range i.e.
        starting more than 40 bp from the start, or 40bp from the end. Anything
        else is discarded.
        
        '''
        
        import logging
        self.logger = logging.getLogger('_regex__')
        
        pair_result = []
        for strand_key, sequence in search_dict.iteritems():
            self.logger.debug("strand key {0}".format(strand_key))
            for orient_key, pattern_list in regex_dict.items():
                for search_pattern in pattern_list:                
                    search_match = search(search_pattern, str(sequence.seq))
                    if search_match:
                        # remove matches within the internal part of read
                        if (search_match.span()[0] < 40) or (search_match.span()[0] > len(sequence.seq) - 40):
                        
                            pair_result.append({'pattern_found' : True, 
                                                'pattern' : search_pattern.pattern,
                                                'start_position' : search_match.span()[0],
                                                'end_position' : search_match.span()[1],
                                                'length' : len(sequence.seq),
                                                'index' : strand_key,
                                                'orient_key' : orient_key})
                            #break
                    #else:
                    #    continue
                    #break
        else:
            try:
                if (pair_result == None) or (len(max(pair_result, key=len)) < 6):
                    pair_result.append({'pattern_found' : False,
                                        'index' : strand_key})    
            except ValueError:
                pair_result.append({'pattern_found' : False,
                                        'index' : strand_key})

        return pair_result
    

    def clip_primers_from_seq(self, search_result, primer_dict, pair_seq_dict, sample_primer_dict, sample_id):
        
        import logging
        self.logger = logging.getLogger('demultip')
        
        new_record = {}
        
        if len(sample_id) == 2:
            raise  ValueError ("Two sample IDs = mispriming")
        else:
            pass  

        r1_slices, r2_slices = self.determine_sequence_slices(sample_id, sample_primer_dict, search_result)
        
        sample_id = sample_id[0]
        
        if search_result[0].get('index') == 'r1':
            record_r1_regx_result = search_result[0]
            record_r2_regx_result = search_result[1]
        elif search_result[0].get('index') == 'r2':
            record_r1_regx_result = search_result[1]
            record_r2_regx_result = search_result[0]
            
        r1_orientation = record_r1_regx_result.get('orient_key')
        r2_orientation = record_r2_regx_result.get('orient_key')
        
        self.logger.debug("r1 orientation {0}, r2 orientation {1}".format(r1_orientation,
                                                                          r2_orientation))
        
        for seq_key, record in pair_seq_dict.iteritems():
            
            if seq_key in record_r1_regx_result.values():       
                
                self.logger.debug("attempt clip of r1 values - seq and qual")
                r1_seq = record.seq[r1_slices.get('r1')[0]:r1_slices.get('r1')[1]]
                r1_quality_scores=record.letter_annotations.get('phred_quality')[r1_slices.get('r1')[0]:r1_slices.get('r1')[1]]
                self.logger.debug("r1 seq clipped...{0} - r1 qual clipped...{1}".format(r1_seq[0:15], r1_quality_scores[0:10]))
                
                if r1_orientation != 'fp':
                    self.alternate_orientation.append([record.id, "\t", sample_id,"\n"])
                #    self.logger.debug("r1 read not in forward orientation...reverse complementing")
                #    r1_seq=self.reverse_complement(str(r1_seq))
                #    r1_quality_scores=r1_quality_scores[::-1]
                #    self.logger.debug("r1 rc sequence starting as...{0}".format(r1_seq[0:25]))
                #    self.logger.debug("r1 rc quality scores starting as...{0}".format(r1_quality_scores[0:15]))

                r1_tmp_record = SeqRecord.SeqRecord(id=record.id,
                            seq=str(r1_seq),
                            description=sample_id,
                            letter_annotations={'phred_quality' : r1_quality_scores})
                self.logger.debug("r1 unclipped length {0}, clipped length {1}".format(len(record.seq), len(r1_tmp_record.seq)))
                
                new_record.setdefault(sample_id +'_r1', []).append(r1_tmp_record)
                self.logger.debug("clipped r1_tmp_record seq...{0}".format(r1_tmp_record.seq[0:25]))
                
            if seq_key in record_r2_regx_result.values():
                
                self.logger.debug("attempt clip of r2 values - seq and qual")
                r2_seq = record.seq[r2_slices.get('r2')[0]:r2_slices.get('r2')[1]]
                r2_quality_scores=record.letter_annotations.get('phred_quality')[r2_slices.get('r2')[0]:r2_slices.get('r2')[1]]
                self.logger.debug("r2 seq clipped...{0} - r2 qual clipped...{1}".format(r2_seq[0:15], r2_quality_scores[0:10]))

                #if r2_orientation != 'rp':
                #    self.logger.debug("r2 read not in reverse orientation...reverse complementing")
                #    r2_seq=self.reverse_complement(str(r2_seq))
                #    r2_quality_scores=r2_quality_scores[::-1]
                #    self.logger.debug("r2 rc sequence starting as...{0}".format(r2_seq[0:25]))
                #    self.logger.debug("r2 rc quality scores starting as...{0}".format(r2_quality_scores[0:15]))
                
                r2_tmp_record = SeqRecord.SeqRecord(id=record.id,
                            seq=str(r2_seq),
                            description=sample_id,
                            letter_annotations={'phred_quality' : r2_quality_scores})
                self.logger.debug("r2 unclipped length {0}, clipped length {1}".format(len(record.seq), len(r2_tmp_record.seq)))
                new_record.setdefault(sample_id +'_r2', []).append(r2_tmp_record)

        self.logger.debug("returning {0} modified seqs".format(len(new_record)))
        
        return new_record

    def correct_orientation_of_reads(self, result_dict):
        
        '''
        
        The result dictionary coming in comes as a list of two dictionary
        items, hopefully with one r1 and one r2 read. The reads may not be
        in the r1=fp r2=rp orientation. It may be the reverse of that which
        is not helpful and needs to be reversed. Therefore we check for that
        and reverse the variables if need be.
        
        [
        {'pattern_found' : True, 
         'pattern' : search_pattern.pattern,
         'start_position' : search_match.span()[0],
         'end_position' : search_match.span()[1],
         'length' : len(sequence.seq),
         'index' : strand_key,
         'orient_key' : orient_key},
        {'pattern_found' : True, 
         'pattern' : search_pattern.pattern,
         'start_position' : search_match.span()[0],
         'end_position' : search_match.span()[1],
         'length' : len(sequence.seq),
         'index' : strand_key,
         'orient_key' : orient_key}
        ]
        
        '''
        
        self.logger = logging.getLogger('_corort_')
        
        if len(result_dict) == 2:
            orientation = {result_dict[0].get('index'): 
                            result_dict[0].get('orient_key'),
                           result_dict[1].get('index') :
                           result_dict[1].get('orient_key')}
            
            if not orientation.get('r1') == 'fp' and not orientation.get('r2') == 'rp':
                self.logger.debug("read is in alternate orientation, swapping around")
                
                result_dict[0]['index'], result_dict[1]['index'] = \
                           result_dict[1]['index'], result_dict[0]['index']
                                  
                return result_dict
            else:
                return result_dict
            
        else:
            return result_dict
              


    def record_buffer_and_writer(self, record_dict):
        
        import logging
        self.logger = logging.getLogger('_buffer_')
        
        buffer_count = []
        
        for key, value in record_dict.items():
            self.record_buffer.setdefault(key, []).append(value)
        
        for entry in self.record_buffer.values(): buffer_count.append(len(entry))
        
        eqn = sum(buffer_count) / float(self.buffer_limit)
        
        # not sure this is as clean as it needs to be
        if (eqn*100 > 99) or (self.starting_total == self.processed_seqs):
            self.logger.debug("Record buffer ({0}) full or analysis finished, writing {1} records to disk".format(self.buffer_limit, sum(buffer_count)))
            
            for pair_key, record in self.record_buffer.items():
                
                if not pair_key == 'discarded':
                    
                    for sample_entries in record:                        
                        for individual_records in sample_entries:
                            filename = os.path.join(self.output_dir, pair_key + ".fq")
                            with open(filename, "a") as output:
                                try:
                                    SeqIO.write(handle=output, sequences=individual_records, format='fastq')
                                except ValueError as e:
                                    self.logger.fatal("{0}".format(e))
                                    self.logger.fatal("{0}".format(individual_records.seq))
                                    raise ("problem writing to file, check record and SeqIO object")
                                
                if pair_key == 'discarded':
                    self.logger.debug("writing {0} discarded records to disk".format(len(record)))
                    for sample_entries in record:
                        for pair_end_key, pair_end_value in sample_entries.items():
                            filename = os.path.join(self.output_dir, pair_end_key + "_discarded" + ".fq")
                            with open(filename, "a") as output:
                                SeqIO.write(handle=output, sequences=pair_end_value, format='fastq')
                                
            #reset the count after the records written to disk                
            buffer_count=None
        
            return "cleared"
    
    
    def screen_read_pair_suitability(self, result):
        '''
        
        {
        'pattern' : search_pattern.pattern,
        'start_position' : search_match.span()[0],
        'end_position' : search_match.span()[1],
        'length' : len(sequence.seq),
        'index' : strand_key,
        'orient_key' : orient_key 
        }
        
        {
        'index': 'r1', 
        'pattern_found': True, 
        'orient_key': 'rp', 
        'pattern': 'ACTGACTGACTAC', 
        'end_position': 13, 
        'length': 301, 
        'start_position': 0
        }
        
        '''
        
        normal_combo = set(['fp', 'rp'])
        rc_combo =set(['fprc', 'rprc'])
        
#        if result[0].get('index') == 'r1':
#            orientation_set = set([result[0].get('orient_key'), 
#                                   result[1].get('orient_key')])
#        elif result[0].get('index') == 'r2':
        try:
            orientation_set = set([result[1].get('orient_key'),
                                   result[0].get('orient_key')])
        except IndexError as e:
            self.logger.debug("Too few primer orientations, failing sequence {0}".format(e))
            return "failed"
        except TypeError as e:
            self.logger.debug("Result is an integer")
        
        if not (len(set.intersection(normal_combo, orientation_set)) == 2 and
                  not len(set.intersection(rc_combo, orientation_set)) == 2):
               
               self.logger.debug("length of set intersection {0}".format(len(set.intersection(normal_combo, orientation_set))))
               self.logger.debug("length of rc set inersection {0}".format(len(set.intersection(rc_combo, orientation_set))))
               
               self.logger.debug("unaccepted primer orientations != 2")
               return "failed"
        else:
            self.logger.debug("accepted primer orientations == 2")
            return "proceed"
        
    
    def run_demultiplex_and_trim(self, opts, **kwargs):
        
        """
            The main part of the script that pulls all the various 
            manipulations together. It takes arguments from the command
            line as well as **kwargs (currently only specifying gzip or not)
        
        """
        
        import logging
        self.logger = logging.getLogger('demultip')
        
        sample_primer_dict = {}
        
        if not opts:
            sys.exit("command line options not getting to main method")
        
        metafile = opts.m
        
        # extract .gz to temp file location
        if 'gzipFilename' in kwargs:
            self.logger.info("Incoming kwargs detected...gzip file?")
            #sequence_file = kwargs.get('gzipFilename')
        else:
            self.logger.info("No kwargs, normal Fastq file")
            #sequence_file = opts.f
        self.logger.info("processing {0} total sequences".format(str((self.r1_tot+self.r2_tot)/4)))
        self.logger.info("using the first {0} bases of primer in search".format(self.search_length))

    
        #extract the relevant data from the metadata file, can maybe change this to non-qiime1
        self.logger.info("Getting header and mapping data...")
        header, mapping_data, run_description, errors, warnings = process_id_map(metafile)
        self.logger.debug("metadata headers {0}".format(header))
        self.logger.debug("csv mapping data from {0}...\n{1}".format(metafile, "\n".join([str(x) for x in mapping_data])))
        
        # get the primer regex search patterns
        self.logger.info("Generating regex search patterns...")
        forward_primers, forward_primers_rc, reverse_primers, reverse_primers_rc = self.create_primer_regex_patterns(header, mapping_data)
        self.primer_pattern_dict_list = {'fp' : forward_primers, 'fprc' : forward_primers_rc, 'rp' : reverse_primers, 'rprc' : reverse_primers_rc}
        
        
        self.logger.debug("forward_primer patterns\n{0}\n".format("\n".join([str(x.pattern) for x in self.primer_pattern_dict_list.get('fp')])))
        self.logger.debug("reverse_primers patterns\n{0}\n".format("\n".join([str(x.pattern) for x in self.primer_pattern_dict_list.get('rp')])))
        
        # replace all extra characters in header with underscore
        intab = '.-+|=:;,&$'
        outtab = '__________'
        trantab = maketrans(intab, outtab)
        
        for samples in mapping_data:
            try:
                sample_primer_dict[samples[header.index('SampleID')].translate(trantab)] = (samples[header.index('LinkerPrimerSequence')], samples[header.index('ReversePrimer')])
            except Exception as e:
                self.logger.error("Can not find {0} in header fields, please make sure metadata file has the required fields".format(e))
                        
        self.logger.debug("sample_primer_dict...{0}".format("\n".join(x) for x in sample_primer_dict.items()))
        self.logger.info("Starting demultiplex process...")
        
        bar = progressbar.ProgressBar(max_value=(self.r1_tot+self.r2_tot)/4,redirect_stdout=True)
    
        for r1, r2 in itertools.izip(self.R1.itervalues(), self.R2.itervalues()):
            #self.logger.debug("r1 {0}".format(r1))
            #self.logger.debug("r2 {0}".format(r1))

            pair_seq_dict = {'r1' : r1, 'r2' : r2}
            self.logger.debug("new read pair\n")
            self.logger.debug("processing new read pair {0}".format(pair_seq_dict.keys()))
            
            self.logger.debug("processing seq ID - R1 {0}... R2 {1}".format(r1.id, r2.id))
            self.logger.debug("R1 sequence - {0}...".format(r1.seq[0:50]))
            self.logger.debug("R2 sequence - {0}...".format(r2.seq[0:50]))

            self.sample_id = ""
            # because we process two sequences at a time (R1 and R2)
            self.processed_seqs += 2

            self.f_primer_found = []
            self.r_primer_found = []
            
            self.logger.debug("Looking in pair read for patterns...")
            
            search_result = self.regex_search_through_sequence(pair_seq_dict, self.primer_pattern_dict_list)
            #self.logger.debug("pre read correction search_result - {0}".format(search_result))
            #search_result = self.correct_orientation_of_reads(search_result)
            #self.logger.debug("post read correction search_result - {0}".format(search_result))
            
            try:
                if type(search_result) == list and len(search_result) > 1:
                    self.logger.debug("search result - {0}".format(search_result[0]))
                    self.logger.debug("search result - {0}".format(search_result[1]))
            except IndexError as e:
                self.logger.debug("search result - {0}".format(search_result))
                self.logger.debug("error in list index {0}".format(e))
            
            
            read_pair_proceed = self.screen_read_pair_suitability(search_result)
            
            self.logger.debug("proceed with read pair ? {0}".format(read_pair_proceed))
            
            if read_pair_proceed != 'failed':
                try:
                    sample_id = self.get_sample_id_from_primer_sequence(sample_primer_dict, 
                                                                        search_result[0].get('pattern'), 
                                                                        search_result[1].get('pattern'))
                    self.logger.debug("- R1 ID -> {0} & R2 ID -> {1} from sample {2}".format(r1.id, r2.id, sample_id))
                except IndexError as e:
                    # sample is missing one or both the patterns keys
                    self.logger.debug("Sample seq is missing a pattern, {0}- discarding read".format(e))
                    output = self.record_buffer_and_writer({'discarded' : pair_seq_dict})
                    self.unmapped_count += 2
                    continue
                try:
                    new_seq = self.clip_primers_from_seq(search_result, self.primer_pattern_dict_list, pair_seq_dict, sample_primer_dict, sample_id)
                    self.logger.debug("clipped read returned...{0} seqs".format(len(new_seq)))
                    output = self.record_buffer_and_writer(new_seq)
                    self.both_primers_count += 2
                except Exception as e:
                    self.logger.debug("attempt to clip sequence failed - errmsg - {0} - discarding read {1}".format(e, output))
                    output = self.record_buffer_and_writer({'discarded' : pair_seq_dict})
                    self.unmapped_count += 2
                    continue
            
                bar.update(self.processed_seqs) 
                
                if output == "cleared":
                    self.record_buffer = {}
                    self.logger.debug("buffer check {0}".format(self.record_buffer))

            elif read_pair_proceed == 'failed':
                self.unmapped_count += 2
                output = self.record_buffer_and_writer({'discarded' : pair_seq_dict})
                bar.update(self.processed_seqs) 
        
        self.logger.info("__________________________")
        self.logger.info("Samples successfully mapped (F+R found): {0}".format(self.both_primers_count))
        self.logger.info("Read pairs in alternate orientation - {0}".format(str(len(self.alternate_orientation))))
        self.logger.info("Sequences not mapped: {0}".format(self.unmapped_count))
        self.logger.info("Total sequences checked: {0}".format(self.processed_seqs))
    
        self.logger.info("writing alternate record IDs...")
        with open("alternate_orientation_records.txt", 'w') as f:
            for sequence_id in self.alternate_orientation:
                output_id = ''.join(sequence_id)
                f.write(output_id)            
    
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
        
        #self.logger.debug("reversing pattern to normal seq {0}".format(sequence))
        
        for key,value in self.inverted_iupac.items():
            if key in sequence:
                sequence = sequence.replace(key, value, 1)
            else:
                # loop does not continue if this is not here, strange?
                continue
        return str(sequence) 
    
    
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
    
    
    def get_sample_id_from_primer_sequence(self, sample_id_primer_dict, f_primer_pattern, r_primer_pattern):
        
        '''
            The actual function responsible for generating the sample ID, being
            the filename into which successfully mapped sequences will be placed.
            
            Once reversed from the regex pattern, the F and R sequences are checked
            against the list of primers provided by the metadata file. Using set
            means that only one sample_id is returned if both match as sets can
            only contain unique items, otherwise there are two items in the set,
            meaning a probably mispriming in which case the sequence is discarded.
        '''
        import logging
        self.logger = logging.getLogger('smplf4id')
        
        self.logger.debug("getting sample ID from the primer sequences")
        self.logger.debug("looking for patterns F - {0} and R - {1}".format(f_primer_pattern, r_primer_pattern))
        
        rectified_f_primer = self.replace_ambiguous_pattern_with_iupac_base(f_primer_pattern)
        
        rectified_r_primer = self.replace_ambiguous_pattern_with_iupac_base(r_primer_pattern)
        self.logger.debug("primer ambiguous bps converted to F - {0} and R - {1}".format(rectified_f_primer, rectified_r_primer))
        
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
        
    # this is not being used
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
            print("directory {0} NOT found, creating directory".format(path))
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
        elif (results.l < 8):
            parser.error("We like to use a minimum of 8bp for primer length")
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
