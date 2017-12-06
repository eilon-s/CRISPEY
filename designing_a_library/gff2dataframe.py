#!/usr/bin/env python2

# https://gist.github.com/slowkow/8101481

from collections import defaultdict
import gzip
import pandas as pd
import re
from roman import fromRoman

GTF_HEADER  = ['chrom', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame'] # ,'attributes'
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')


def gff2dataframe(filename):
    """Open an optionally gzipped GTF file and return a pandas.DataFrame.
    """
    # Each column is a list stored as a value in this dict.
    result = defaultdict(list)

    for i, line in enumerate(lines(filename)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines.
            if key not in result:
                result[key] = [None] * i

        # Ensure this row has some value for each column.
        for key in result.keys():
            result[key].append(line.get(key, None))

    return pd.DataFrame(result)


def lines(filename):
    """Open an optionally gzipped GTF file and generate a dict for each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open

    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                yield parse(line)


def parse(line):
    """Parse a single GTF line and return a dict.
    """
    result = {}

    fields = line.rstrip().split('\t')

    for i, col in enumerate(GTF_HEADER):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

    for i, info in enumerate(infos, 1):
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value.
        if value:
            result[key] = _get_value(value)

    return result


def _get_value(value):
    if not value:
        return None

    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return value
  


def sgd_gff2dataframe(input_annot_gff_filename, ok_features_list, include_dubious=False):
    """
    loads sgd annotation gff file to dataframe\
   leave only features of the types in the input list
    filters mitochondrial and 2-micron genes
    """
    gff_df = gff2dataframe(input_annot_gff_filename)

    gff_df['Name'] = gff_df['Name'].str.strip()

    #gff_df = gff_df[gff_df['feature'].isin(['gene','CDS'])  ]
    # just coding sequence
    gff_df = gff_df[gff_df['feature'].isin(ok_features_list)  ]   # ['CDS']

    # filtering dubious and not in standard chromosomes
    if include_dubious:
      idx =  ~gff_df['chrom'].isin(['chrMito', '2-micron'])
    else:
      idx = gff_df['Note'].apply(lambda x: not str(x).startswith('Dubious')) & ~gff_df['chrom'].isin(['chrMito', '2-micron'])
    gff_df = gff_df[idx]

    # fixing chromosome to be chrXX where XX is integer 
    gff_df['chrom'] = gff_df['chrom'].apply(lambda x: 'chr' + str(fromRoman(str(x).replace('chr',''))).zfill(2) )

    # removing "SGD:" from the id column
    gff_df['dbxref'] = gff_df['dbxref'].apply(lambda x: str(x).replace('SGD:',''))
    gff_df['start'] = gff_df['start'].apply(int)
    

    
    return(gff_df)


 

 
 