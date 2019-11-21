"""
This module keeps properties of annotation files and functions associated with those.
"""
import logging
import os

import plastid
from fivepseq.util.writers import FivePSeqOut
from preconditions import preconditions

from fivepseq import config


class Annotation:
    file_path = None
    gene_filter = None
    gene_filter_attribute = "gene_id"

    # the transcript_assembly_dict is a two-level dictionary.
    # the first level provides keys to gene-filters;
    # the second level provides keys to span-sizes
    transcript_assembly_dict = {}

    PROTEIN_CODING = "protein_coding"
    default_filter = PROTEIN_CODING

    gs_transcript_dict = None
    gs_transcriptInd_dict = None
    gs_geneID_attribute = None
    span_size = 0

    logger = logging.getLogger(config.FIVEPSEQ_LOGGER)

    @preconditions(lambda file_path: isinstance(file_path, str),
                   lambda transcript_assembly: isinstance(transcript_assembly, list))
    def __init__(self, transcript_assembly, file_path):

        """
        Initiates an Annotation instance with transcript_assembly list.
        Raises a ValueError if transcript assembly is None

        :param transcript_assembly:
        """

        self.file_path = file_path
        if transcript_assembly is None:
            error_message = "transcript_assembly is None. Cannot instantiate an Annotation object."
            self.logger.error(error_message)
            raise ValueError(error_message)
        self.transcript_count = len(transcript_assembly)

        self.transcript_assembly_dict.update({
            self.default_filter: {
                0: transcript_assembly
            }})  # the trnascript assembly returned by the reader is set for 0 span size and no filter

        self.logger.debug("Annotation object successfully created from file %s" % file_path)

    def set_permanent_gene_filter(self, gene_filter_file):
        """
        Provide the file where the names of genes to be filtered are present.
        The genes should be named according to the attribute in the gff file.
        The default "Name" attribute is used, unless provided otherwise.

        :param geneset_filter: the file containing genes separated by new lines
        :param attribute: the attribute name in which genes are presented
        :return:
        """

        self.gene_filter = self.read_gene_filter_file(gene_filter_file)
        self.apply_permanent_gene_filter(self.gene_filter, self.gene_filter_attribute)

    def apply_permanent_gene_filter(self, gene_filter, attribute):
        gene_filtered_assembly = []
        for transcript in self.get_transcript_assembly(span_size=0):
            attr_value = FivePSeqOut.get_transcript_attr(transcript, attribute)
            if attr_value in gene_filter:
                gene_filtered_assembly.append(transcript)
            # TODO this is not a universal solution, but when the transcripts have names with -1 in the end this works
            elif len(attr_value.split(":")) > 1 and attr_value.split(":")[
                1] in gene_filter:  # gene_id filtering results in IDs in the form gene:xxx
                gene_filtered_assembly.append(transcript)
            elif attr_value.split("-")[0] in gene_filter:
                gene_filtered_assembly.append(transcript)
            elif attr_value.split(".")[0] in gene_filter:
                gene_filtered_assembly.append(transcript)

        if len(gene_filtered_assembly) == 0:
            raise Exception("None of the genes in the geneset filter were present in the annotation file")

        #NOTE the gene filter is permanently applied in the beginning,
        # thus span sizes other than 0 are not taken into consideration
        self.transcript_assembly_dict.update({self.default_filter: {0: gene_filtered_assembly}})

    def read_gene_filter_file(self, gene_filter_file):

        if not os.path.exists(gene_filter_file):
            raise Exception("The gene set file %s does not exist" % gene_filter_file)

        gene_filter = []
        with open(gene_filter_file) as file:
            self.gene_filter_attribute = file.readline().rstrip("\n\r")
            self.logger.info("Gene filter attribute set as: %s" % self.gene_filter_attribute)

            line = file.readline()
            count = 1
            while line:
                # TODO check if this can work with spaces, as they are common
                # if " " in line or "\t" in line:
                #    raise Exception("The gene set file %s should not contain spaces or tabs. Found one in line %d"
                #                    % (gene_filter_file, count))

                gene_filter.append(line.rstrip("\n\r"))
                line = file.readline()
                count += 1

        return gene_filter

    def set_default_span_size(self, span_size):
        """
        Sets the default span size to be used with transcript assembly retrieval with no span size specifications.
        """

        self.span_size = span_size

    def get_transcript_assembly(self, span_size=None, geneset_filter=None):
        """
        Returns transcripts if a transcript assembly already exists with provided span_size and filter.
        Creates and returns a new transcript assembly otherwise.

        :param span_size:  int span_size - defines the number of nucleotides to span
        :param geneset_filter: filter name to be applied on transcripts
        :return: a vector of transcript or an iterable of type Transcript
        """

        if span_size is None:
            span_size = self.span_size

        if geneset_filter is None:
            geneset_filter = self.default_filter

        if geneset_filter not in self.transcript_assembly_dict:
            raise Exception("Provided geneset filter %s is not processed yet" % geneset_filter)

        if span_size in self.transcript_assembly_dict[geneset_filter]:
            return self.transcript_assembly_dict[geneset_filter][span_size]

        return self.generate_transcript_assembly(span_size, geneset_filter)

    @preconditions(lambda span_size: isinstance(span_size, int))
    def generate_transcript_assembly(self, span_size, geneset_filter=None):
        """
        Adds spanning regions to the start and the end of the transcript and yields the transcripts one-by-one.
        The original transcript assembly is left intact.

        :param span_size:  int span_size - defines the number of nucleotides to span
        :param geneset_filter: filter to apply on transcripts

        :return: transcript vector
        """

        # error check
        if geneset_filter is None:
            transcript_filter = self.default_filter

        if span_size < 0:
            error_message = "Negative span_size of %d provided: cannot span the transcripts with negative or " \
                            "zero values" % span_size
            self.logger.error(error_message)
            raise ValueError(error_message)
        self.logger.debug("Transcript span size: %d" % span_size)

        # info 
        self.logger.info(
            "Generating transcript assembly under filter: %s and span size: %d"
            % (geneset_filter, span_size))

        # generation
        this_transcript_assembly = []
        for transcript in self.get_transcript_assembly(span_size=0, geneset_filter=geneset_filter):
            # apply filter
#            if transcript_filter == self.FORWARD_STRAND_FILTER:
#                if transcript.strand != '+':
#                    continue

#            elif transcript_filter == self.REVERSE_STRAND_FILTER:
#                if transcript.strand != '-':
#                    continue

            # span
            spanned_transcript = self.span_transcript(transcript, span_size)

            # append 
            this_transcript_assembly.append(spanned_transcript)

        # add assembly to dictionary
        self.transcript_assembly_dict[geneset_filter].update({span_size: this_transcript_assembly})

        # info
        self.logger.info("Transcript assembly generated")

        return this_transcript_assembly

    @preconditions(lambda span_size: isinstance(span_size, int),
                   lambda span_size: span_size >= 0)
    def span_transcript(self, transcript, span_size):
        """
        Spans the transcript symmetrically by the given span size.
        Spanning may result in going beyond boundaries of current genome, with negative and out of range coordinates,
        which should be handled explicitly in other functions.

        :param transcript:
        :param span_size:
        :return:
        """
        span = transcript.spanning_segment

        new_start = span.start - span_size

        start_span = plastid.GenomicSegment(
            span.chrom, new_start, span.start, span.strand)

        end_span = plastid.GenomicSegment(
            span.chrom, span.end, span.end + span_size, span.strand)

        spanned_transcript = transcript.__deepcopy__(None)
        spanned_transcript.add_segments(start_span, end_span)

        return spanned_transcript

    def store_gene_sets(self, geneset_file):
        """
        Provide the file where the names of genes and their mappings to one or several gene sets.
        The file should be a tab delimited file, with the first column containing gene IDs and the second one - gene sets.
        The header of the first column should be equal to the attribute in the gff file, which has been used for gene ID specification.
        The header of the second column should be named "geneset"

        GFF_ATTRIBUTE_NAME(e.g.: gene_id)  geneset
        gene1   GS1
        gene2   GS1
        gene3   GS2
        gene4   GS2


        :param geneset_file: a tab-delimited file containing gene-geneset mapping
        :return: stores and returns a {geneset,[transcript]} dictionary
        """

        # check file
        if not os.path.exists(geneset_file):
            raise Exception("The gene set file %s does not exist" % geneset_file)

        # read genesets, store in GS:[geneIDs] dictionary
        gs_dict = {}
        geneIDs = []
        with open(geneset_file) as file:
            header = file.readline()
            tokens = header.split("\t")
            if len(tokens) != 2:
                raise Exception("The geneset file should have exactly two columns. Found %d instead."
                                % len(tokens))
            attribute = tokens[0]
            self.gs_geneID_attribute = attribute
            self.logger.info("The attribute %s will be used to read gff for setting the geneset dictionary"
                             % attribute)
            line = file.readline()
            count = 1
            while line:
                # if " " in line:
                #    raise Exception("The geneset file %s should not contain spaces. Found one in line %d"
                #                    % (geneset_file, count))
                tokens = line.split('\t')
                if len(tokens) != 2:
                    raise Exception(
                        "The geneset file should have exactly two columns. Found %d in line %d instead."
                        % (len(tokens), count))
                geneID = tokens[0]
                geneIDs.append(geneID)
                geneset = tokens[1].rstrip("\n\r")

                if gs_dict.has_key(geneset):
                    gs_dict[geneset].append(geneID)
                else:
                    gs_dict.update({geneset: [geneID]})

                line = file.readline()
                count += 1

        # check if the gene IDs are in the default transcript assembly
        # for those that are, create geneID:transcript dictionary

        geneID_transcript_dict = {}
        geneID_transcriptInd_dict = {}
        transcript_ind = 0
        for transcript in self.get_transcript_assembly(span_size=0):
            attr_value = FivePSeqOut.get_transcript_attr(transcript, attribute)
            # TODO this is not a universal solution, but when the transcripts have names with -1 in the end this works
            geneID = None
            if attr_value in geneIDs:
                geneID = geneIDs[geneIDs.index(attr_value)]
            elif len(attr_value.split(":")) > 1 and attr_value.split(":")[1] in geneIDs:
                geneID = geneIDs[geneIDs.index(attr_value.split(":")[1])]
            # attr_value.split("-")[0] in geneIDs or \
            # attr_value.split(".")[0] in geneIDs:
            if geneID is not None:
                geneID_transcript_dict.update({geneID: transcript})
                geneID_transcriptInd_dict.update({geneID: transcript_ind})

            transcript_ind = transcript_ind + 1

        if len(geneID_transcript_dict) == 0:
            raise Exception("None of the genes in the geneset file were present in the annotation file")

        # with those geneIDs that mapped to actual transcripts,
        # store a {GS: [transcripts]} dictionary
        #TODO remove if fine
        #gs_transcript_dict = {}
        gs_transcriptInd_dict = {}
        for gs in gs_dict.keys():
            #gs_transcript_dict.update({gs: []})
            self.transcript_assembly_dict.update({gs:{0:[]}})
            gs_transcriptInd_dict.update({gs: []})
            for geneID in gs_dict[gs]:
                if geneID_transcript_dict.has_key(geneID):
                    #TODO remove if fine
                    #gs_transcript_dict[gs].append(geneID_transcript_dict[geneID])
                    self.transcript_assembly_dict[gs][0].append(geneID_transcript_dict[geneID])
                    gs_transcriptInd_dict[gs].append(geneID_transcriptInd_dict[geneID])

        #TODO remove if fine
        #self.gs_transcript_dict = gs_transcript_dict
        self.gs_transcriptInd_dict = gs_transcriptInd_dict

        geneset_summary = ""
        for gs in gs_transcriptInd_dict.keys():
            geneset_summary += gs + "\n"

        self.logger.info("Genesets processed. %d out of %d unique geneIDs matched corresponding transcripts."
                         % (len(set(geneID_transcript_dict.keys())), len(set(geneIDs))))
        self.logger.info("Found genesets:\n%s\n" % geneset_summary)


    def remove_geneset_filter(self):
        #TODO remove if fine
        #self.transcript_assembly_dict.update({self.protein_coding_filter: {0: self.transcript_assembly_default}})
        self.default_filter = self.PROTEIN_CODING

    def apply_geneset_filter(self, gs):
        if not self.gs_transcriptInd_dict.has_key(gs):
            raise Exception("The annotation does not have a filter named %s" % gs)

        self.default_filter = gs

        #TODO remove if fine
        """
        gene_filtered_assembly = self.gs_transcript_dict[gs]

        if len(gene_filtered_assembly) == 0:
            raise Exception("No genes remain after applying the geneset filter %s " % gs)
        self.transcript_assembly_default = self.transcript_assembly_dict.get(self.protein_coding_filter). \
            get(0)
        self.transcript_assembly_dict.update({self.protein_coding_filter: {0: gene_filtered_assembly}})
        """