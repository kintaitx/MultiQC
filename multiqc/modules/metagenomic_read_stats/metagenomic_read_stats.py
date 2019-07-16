from __future__ import absolute_import

from .metagenomic_read_stats import MultiqcModule


#!/usr/bin/env python

""" MultiQC module to parse output from the human read removal stage of Kintai's metagenomics pipeline """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """ Metagenomic read stats """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Read Count Stats', anchor='ReadCountStats',
        info="Results of human read removal")

        # Find and load any Kallisto reports
        self.counts_data = dict()
        for f in self.find_log_files('metagenomic_read_stats', filehandles=True):
            log.info
            self.parse_stats_log(f)

        # Filter to strip out ignored sample names
        self.counts_data = self.ignore_samples(self.counts_data)

        if len(self.counts_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.counts_data)))

        # Write parsed report data to a file
        self.write_data_file(self.counts_data, 'multiqc_read_counts_stats')

        # Basic Stats Table
        self.metagenomic_read_stats_general_stats_table()

        # Alignment Rate Plot
        self.add_section( plot = self.counts_bar_plot() )


    def parse_stats_log(self, f):
        s_name = total_reads = paligned_reads = fraglength = None
        #get the sample name from the file name
        s_name = self.clean_s_name(f['fn'], f['root'])
        if s_name in self.counts_data:
            log.debug("Duplicate sample name found")
        self.counts_data[s_name] =  {}
        headers = f['f'].readline().split("\t")
        values = f['f'].readline().split("\t")
        for i in range(0, len(headers)):
            self.counts_data[s_name][headers[i]] = float(values[i])
        if 'total_reads' in self.counts_data[s_name] and 'non_human_reads' in self.counts_data[s_name]:
            self.counts_data[s_name]['human_reads'] = self.counts_data[s_name]['total_reads'] - \
                                                      self.counts_data[s_name]['non_human_reads']
        if 'fraction_human' in self.counts_data[s_name]:
            self.counts_data[s_name]['fraction_non_human'] = 1.0 - self.counts_data[s_name]['fraction_human']
        return




    def metagenomic_read_stats_general_stats_table(self):
        """ Take the parsed stats from the Kallisto report and add it to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['non_human_reads'] = {
            'title': 'Non Human Read Count',
            'description': 'Number of non-human reads',
        }
        headers['human_reads'] = {
            'title': 'Human Read Count',
            'description': 'Number of Human Reads',
        }
        headers['fraction_non_human'] = {
            'title': 'Fraction Human Reads',
            'description': 'Fraction of reads mapping to the human genome',
        }
        self.general_stats_addcols(self.counts_data, headers)




    def counts_bar_plot(self):
        """ Make the HighCharts HTML to plot the alignment rates """

        # Specify the order of the different possible categories
        keys['human_reads'] = {'color': '#437bb1', 'name': 'Human Reads'}
        keys['non_human_reads'] = {'color': '#7cb5ec', 'name': 'Mapped to multiple loci'}

        # Config for the plot
        pconfig = {
            'id': 'Read_Count_Plot',
            'title': 'Mapping Counts',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        return bargraph.plot(self.counts_data, keys, pconfig)

