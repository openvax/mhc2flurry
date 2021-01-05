# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
import os
import sys
import time

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from datacache import get_data_dir
from pyensembl import EnsemblRelease
import pandas as pd
import shellinford

from mhc2flurry.downloads import get_path


logger = logging.getLogger(__name__)


parser = ArgumentParser(
    description=__doc__,
    formatter_class=RawDescriptionHelpFormatter)

parser.add_argument(
    '--input-csv',
    default='',
    help='Path to CSV containing peptide DataFrame to annotate. Must contain "peptide" column')
parser.add_argument(
    '--output-csv',
    default='',
    help='Desired path to output CSV, will contain input dataframe with additional column denoting '
         'transcript IDs')
parser.add_argument(
    '--output-transcript-fasta',
    default='',
    help='Desired path to output transcript FASTA, will contain sequences for all transcripts '
         'found to match any sequence in the input')


def generate_all_transcript_sequences(genome):
    """
    Generates all coding transcript sequences from the input genome. Each sequence will have
    the structure: "id|sequence"

    Parameters
    ----------
    genome : pyensembl.EnsemblRelease
        Input genome to load for reference peptides
    """
    for t in genome.transcripts():
        if t.is_protein_coding:
            protein_sequence = t.protein_sequence
            transcript_id = t.id
            yield '%s|%s' % (transcript_id, protein_sequence)


def load_reference_peptides_index(genome, force_reload=False):
    """
    Loads the FM index containing reference peptides.

    Parameters
    ----------
    genome : pyensembl.EnsemblRelease
        Input genome to load for reference peptides

    force_reload : bool, optional
        If true, will recompute index for this genome even if it already exists.

    Returns
    -------
    fm : shellinford.FMIndex
        Index populated with reference peptides from the genome
    """
    cache_dir = get_data_dir()
    path = os.path.join(cache_dir, 'homo_sapiens_102_with_transcript_id.fm')

    if force_reload or not os.path.exists(path):
        start = time.time()
        logger.info('Building FM index at %s', path)
        fm = shellinford.FMIndex()
        fm.build(generate_all_transcript_sequences(genome), path)
        end = time.time()
        logger.info('Done building FM index, took %.1f minutes' % ((end - start)/60))
        return fm
    else:
        logger.info('Loading existing FM index from %s' % path)

    return shellinford.FMIndex(filename=path)


def get_matching_transcript_ids(peptide, fm):
    """
    Searches the input FM index for the input peptide. Returns a set of matching comma-separated
    transcript IDs.
    """
    docs = fm.search(peptide)
    ids = [doc.text.split('|')[0] for doc in docs]
    return ' '.join(ids)


def main(args_list=None):
    if args_list is None:
        args_list = sys.argv[1:]
    args = parser.parse_args(args_list)

    # create pyensembl genome object
    genome = EnsemblRelease(102)
    fm = load_reference_peptides_index(genome, force_reload=False)

    df = pd.read_csv(args.input_csv)
    if 'peptide' not in df.columns:
        raise ValueError('Input DataFrame must contain "peptide" column')

    logger.info('Annotating peptides...')
    start = time.time()
    df['transcripts'] = df.peptide.apply(lambda x: get_matching_transcript_ids(x, fm))
    end = time.time()
    logger.info('Done annotating peptides, took %.1f minutes' % ((end - start)/60))

    logger.info('Writing out to CSV...')
    start = time.time()
    df.to_csv(args.output_csv, index=False)
    end = time.time()
    logger.info('Done writing CSV, took %.1f minutes' % ((end - start)/60))

    # write out FASTA with all transcript sequences occurring here
    logger.info('Creating FASTA with transcript sequences...')
    all_transcripts = set()
    for transcripts in df['transcripts']:
        all_transcripts.update(transcripts.split(' '))

    with open(args.output_transcript_fasta, 'w') as f:
        for transcript in all_transcripts:
            doc = fm.search(transcript)
            if len(doc) != 1:
                logger.warning('More than one record matching transcript %s' % transcript)
            transcript_id, sequence = doc[0].text.split('|')
            f.write('>%s\n%s\n' % (transcript_id, sequence))
    logger.info('Done.')


if __name__ == '__main__':
    main()
