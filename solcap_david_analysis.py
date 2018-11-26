'''Joint analysis of the merged dataset from SolCAP

'''
from __future__ import division

import os.path

from os import listdir
from os.path import join as pjoin
from os.path import exists
from StringIO import StringIO
from csv import DictReader, DictWriter, Dialect, QUOTE_MINIMAL
from collections import Counter, defaultdict
import math
import random
from pprint import pprint

import numpy
from scipy.stats import distributions
from scipy.stats.stats import chisquare
from scipy.interpolate import interp2d
from scipy.spatial.distance import squareform, pdist

import simplekml

from pandas import Series, DataFrame, MultiIndex, read_csv

from solcap_snp_coding.snp_syns import get_snp_solcap_id
from solcap_snp_coding.snp_utils import (read_genotypes,
                                      write_genotypes_to_csv,
                                      change_genotypes_strand,
                                      change_ab_genotypes_strand,
                                      CUSTOMER_STRAND, DESIGN_STRAND,
                                      AB_CODING)

from variation.inout.dataframe import load_dataframe_csv
from variation.inout.genetic import (load_codominant_genetic_csv,
                                     load_markers_map, GenotypeCodec,
                                     MISSING_INT_GENOTYPE,
                                     get_codec_from_genotypes,
                                     write_genetic_csv)
from variation.matrixwithmeta import (PLOIDY, BITS_PER_ALLELE, ALLELE_CODING,
                                      MatrixWithMeta, load_matrix,
                                      INDIVIDUALS_IN_ROWS, LOCATION_COL,
                                      MOLECULE_COL, CLASSIFICATION_COL)
from variation.analyses.matrix_tools import (hstack_dframes,
                                             transpose_genotypes,
                                             filter_matrix_cols,
                                             filter_matrix_rows,
                                             sorted_marker_map,
                                             dframe_append_columns,
                                             vstack_dframes)
from variation.analyses.smartpca import do_smartpca
from variation.plot import (get_scatter_groups_from_matrix, scatter_groups,
                            geographic_map, stacked_bars, plot_matrix,
                            scatter_series, ColorMap, MarkerMap,
                            html_to_rgb_color)
from variation.analyses.phase import run_fast_phase, parse_fast_phase_results
from variation.analyses.pop_stats import (calculate_individual_heterozigosity,
                                          calculate_basic_stats,
                                          PopPolymorphism, _PopCounts,
                                          generate_plot_groups_for_rarefaction,
                                          AlleleFreqsDistribs)
from variation.analyses.adegenet import create_file_snp_fomat
from variation.analyses.finestructure import create_chromopainter_infiles
from variation.analyses.geography import lat_to_deg, lon_to_deg, split_in_grids
from variation.analyses.distance import (PopDistances,
                                         pairwise_individual_genetic_distance,
                                         pairwise_geographic_distance,
                                         correlation_distance_matrices)
from variation.analyses.trees import calc_nj
from variation.analyses.structure import (do_structure_for_ks,
                                          get_default_params,
                                          create_param_files,
                                          create_genotype_file,
                                          read_structure_results,
                                          sort_individual_classification)
from variation.analyses.adze import do_rarefaction
from variation.analyses.multivariate import do_pca
from variation.analyses.linkage_desequilibrium import calculate_unphased_LD
from accessions.abucode import get_comav_acc_codes

from decoder.decoder import (plot_color_chromosomes, compare_groups,
                             comparison_results_to_igv)


# pylint: disable=C0111

ANALYSIS_DIR = os.path.expanduser('~/analyses/solcap_david')
RESULTS_DIR = pjoin(ANALYSIS_DIR, 'results')
ORIG_DATA_DIR = pjoin(ANALYSIS_DIR, 'original_data')

WORLD_CLIMATE_IMAGE = pjoin(ORIG_DATA_DIR, 'World_Koppen_Map.recortado.png')
WORLD_CLIMATE_IMAGE = pjoin(ORIG_DATA_DIR, 'World_Koppen_Map.flipped.png')

SOLCAP_DESIGN_FPATH = pjoin(ORIG_DATA_DIR,
'solcap_web.tomato_ng_snp_10k_infinium_annotation_v2.10_FINAL_w_solcap_ids.csv')

BAUCHET_PIMPI_DIR = pjoin(ORIG_DATA_DIR, 'gbauchet')
BAUCHET_PIMPI_GENOTYPE_FPATH = pjoin(BAUCHET_PIMPI_DIR,
                                    'GBauchet_genotypes_AB.csv')
BAUCHET_PIMPI_PASSPORTS_FPATH = pjoin(BAUCHET_PIMPI_DIR,
                                    'GBauchet_passports.csv')
BAUCHET_PIMPI_GENOTYPE_FPATH_DESIGN = pjoin(BAUCHET_PIMPI_DIR,
                                            'GBauchet_passports.DESIGN.csv')

SMALL_CERAS_DIR = pjoin(ORIG_DATA_DIR, 'small_ceras')
SMALL_CERAS_GENOTYPE_FPATH_ORIG = pjoin(SMALL_CERAS_DIR,
                                        'FullDataTable_AB_upv.csv')
SMALL_CERAS_GENOTYPE_FPATH_DESIGN = pjoin(SMALL_CERAS_DIR,
                                          'small_ceras_genotypes.DESIGN.csv')
SMALL_CERAS_PASSPORTS_FPATH = pjoin(SMALL_CERAS_DIR,
                                    'small_ceras.PASSPORTS.csv')

INRA_ORIG_DATA_DIR = pjoin(ORIG_DATA_DIR, 'inra')
INRA_GENOTYPE_FPATH = pjoin(INRA_ORIG_DATA_DIR, 'data_SOLCAP_INRA.csv')
INRA_GENOTYPE_FPATH_DESIGN = pjoin(INRA_ORIG_DATA_DIR,
                                   'data_SOLCAP_INRA.DESIGN.csv')
INRA_PASSPORTS_FPATH = pjoin(INRA_ORIG_DATA_DIR,
                             'passport_data_SOLCAP_INRA.csv')

COMAV_ORIG_DATA_DIR = pjoin(ORIG_DATA_DIR, 'plos')
COMAV_GENOTYPE_FPATH = pjoin(COMAV_ORIG_DATA_DIR,
                             'genotypes_all_experiments.csv')
COMAV_PASSPORTS_FPATH = pjoin(COMAV_ORIG_DATA_DIR,
                              'pasportes_banco.r069_r068.csv')
PLOS_PASSPORTS_FPATH = pjoin(COMAV_ORIG_DATA_DIR,
                             'journal.pone.0048198.s003.csv')

MERGED_SOLCAP_DESIGN_NON_FILTERED_GENOTYPES = pjoin(RESULTS_DIR,
                             'merged_non_filtered_design_solcap_genotypes.csv')

KOENING_RESULT_DIR = pjoin(RESULTS_DIR, 'koenig')
KOENING_ORIG_DIR = pjoin(ORIG_DATA_DIR, 'koenig')
KOENING_GENOTYPES_FPATH = pjoin(KOENING_RESULT_DIR,
                                'koenig_solcap_genotypes.csv')
KOENIG_PASSPORTS_FPATH = pjoin(KOENING_ORIG_DIR, 'passport_data_koenig.csv')

ACC_SYNS_FPATH = pjoin(COMAV_ORIG_DATA_DIR,
                       'solcap_acc_samples_synonyms.MODIFIED.csv')
DAVID_ORIG_DATA_DIR = pjoin(ORIG_DATA_DIR, 'david')
DAVID_GENOTYPE_FPATH = pjoin(DAVID_ORIG_DATA_DIR,
                             'genotypes.journal.pone.0045520.s004.csv')
DAVID_PASSPORTS_FPATH = pjoin(DAVID_ORIG_DATA_DIR,
                            'passports_journal.pone.0045520.s003.MODIFIED.csv')
# DAVID_PASSPORTS_FPATH = pjoin(DAVID_ORIG_DATA_DIR,
#                               'passports_journal.pone.0045520.s003.csv')

DAVID2_ORIG_DATA_DIR = pjoin(ORIG_DATA_DIR, 'david2')
DAVID2_GENOTYPE_FPATH = pjoin(DAVID2_ORIG_DATA_DIR, 'david2_genotypes.csv')
DAVID2_PASSPORTS_FPATH = pjoin(DAVID2_ORIG_DATA_DIR, 'david2_passports.csv')
DAVID2_GENOTYPE_FPATH_DESING = pjoin(DAVID2_ORIG_DATA_DIR,
                                     'david2_genotypes.DESIGN.csv')


SHAPE_DIR = pjoin(ORIG_DATA_DIR, 'shape_genotypes')
FAS_GENOTYPES_FPATH = pjoin(SHAPE_DIR, 'fas_genotypes.csv')
FW32_GENOTYPES_FPATH = pjoin(SHAPE_DIR, 'fw32_genotypes.csv')
FW113_GENOTYPES_FPATH = pjoin(SHAPE_DIR, 'fw113_genotypes.csv')
FW22_GENOTYPES_FPATH = pjoin(SHAPE_DIR, 'fw22_genotypes.csv')
LC_HAPLOTYPES_FPATH = pjoin(SHAPE_DIR, 'lc_haplotypes.csv')
OVATE_GENOTYPES_FPATH = pjoin(SHAPE_DIR, 'ovate_genotypes.csv')
SUN_GENOTYPES_FPATH = pjoin(SHAPE_DIR, 'results_sun.csv')

ESTHER_SHAPE_GENOTYPES_FPATH = pjoin(SHAPE_DIR,
                                     'Esther_fruit_genotyped_accs.csv')
ESTHER_SHAPE_GENOTYPES_SYNONYMS_FPATH = pjoin(SHAPE_DIR, 'Esther_synonyms.csv')
ESTHER_2_SHAPE_GENOTYPES_FPATH = pjoin(SHAPE_DIR,
                                       'Esther_2_fruit_genotypes.csv')
INRA_SHAPE_GENOTYPES_FPATH = pjoin(SHAPE_DIR, 'INRA_fruit_genotypes.csv')
INRA_OVATE_GENOTYPES_FPATH = pjoin(SHAPE_DIR, 'ovate_genotypes_inra.csv')

PROCESSED_DATA_DIR = pjoin(RESULTS_DIR, 'processed_data')
# PROCESSED_DATA_DIR = pjoin(RESULTS_DIR, 'processed_data_2')
SAVED_MERGED_FILTERED_GENOTYPES_FPATH = pjoin(PROCESSED_DATA_DIR,
                                        'merged_filtered_genotypes.p')
SAVED_PHYSICAL_MAP_FPATH = pjoin(PROCESSED_DATA_DIR, 'physical_map.p')
SAVED_GENETIC_MAP_FPATH = pjoin(PROCESSED_DATA_DIR, 'genetic_map.p')
SAVED_MERGED_PASSPORTS_FPATH = pjoin(PROCESSED_DATA_DIR, 'merged_passports.p')

RECOMBINATION_MAP_DIR = pjoin(RESULTS_DIR, 'chrom_recombination')
INTERPOLED_RECOM_MAP = pjoin(RECOMBINATION_MAP_DIR,
                             'interpolated_genetic_map.csv')

LONG_ACCESSION_NAMES = {'Brandywine (Sudduth/Quisenberry)': 'Brandywine',
                      '1091-Chonto 21 (mataverde (3-21-2))': '1091-Chonto 21',
                        'Pomodoro Superselezione di Marmande': 'Pomodoro'}
SYMBOLS_TO_REPLACE_IN_ACCESSION_NAMES = {'(': '_', ')': '_',
                                         ' ': '_'}

PASSPORT_COUNTRY_SYNS = {'Bolivia': 'BOL', 'Brazil': 'BRA', 'CANADA': 'CAN',
                         'Colombia': 'COL', 'Costa Rica': 'CRI', 'Cuba': 'CUB',
                         'Czechoslovakia': 'TCH', 'Ecuador': 'ECU',
                         'El Salvador': 'SLV', 'England': 'ENG',
                         'France': 'FRA', 'Honduras': 'HND',
                         'Indonesia': 'IDN', 'Italy': 'ITA', 'Mexico': 'MEX',
                         'Netherlands': 'NLD', 'Peru': 'PER', 'Poland': 'POL',
                         'Russia': 'RUS', 'Spain': 'ESP',
                         'United States': 'USA'}

PASSPORT_SPECIES_SYNS = {'S. lycopersicum ': 'tomato',
                         'S. lycopersicum': 'tomato',
                         'S. lycopersicon': 'tomato',
                         'lycopersicum': 'tomato',
                         'S. lycopersicum \'cerasiforme\'': 'ceras',
                         'S. pimpinellifolium': 'pimpi',
                         'S. pimpinellifolium ': 'pimpi',
                         ' S. pimpinellifolium ': 'pimpi',
                         'pimpinellifolium': 'pimpi',
                         'S. pimpinellifolium x S. lycopersicum': 'pimpixtom',
                         'S. penellii x S. lycopersicum': 'penelliixtom',
                         'cerasiforme': 'ceras'}

# these samples are duplicated between comav's and david's datasets:
REPEATED_SAMPLES_COMAV = {'UC-82': 'UC-82_comav', 'LA1545': 'LA1545_comav',
                    'LA1549': 'LA1549_comav', 'LA2312': 'LA2312_comav',
                    'LA2675': 'LA2675_comav', 'LA1547': 'LA1547_comav'}

REPEATED_SAMPLES_BAUCHET = {'LA0373': 'LA0373_inra',
                            'LA2181': 'LA2181_inra',
                            'LA1582': 'LA1582_inra',
                            'LA1301': 'LA1301_inra',
                            'LA1388': 'LA1388_inra',
                            'LA1547': 'LA1547_inra',
                            'LA1933': 'LA1933_inra'}

REPEATED_SAMPLES_KOENIG = {'LA1589': 'LA1589_koenig'}

REPEATED_SAMPLES_DAVID2 = {'RioGrande': 'RioGrande_david',
                           'Yellow_Stuffer': 'Yellow_Stuffer_david',
                           'LA2561': 'LA2561_david_new'}

SMALL_CERAS_CODES = {'CH2': 'ECU996',
                     'CH3': 'ECU1015',
                     'CH4': 'ECU1018',
                     'CH5': 'ECU1019',
                     'CH6': 'ECU1023',
                     'CH7': 'ECU1024',
                     'CH8': 'ECU1067',
                     'CH1': 'ECU1104',
                     'CH9': 'ECU1105',
                     'CH10': 'ECU1107',
                     'CH11': 'ECU1121',
                     'CH12': 'ECU1142'}

GROUPING_ORDER = {'species': ['SC', 'SP', 'SLC', 'SLL', 'mixture'],
         'group1': ['SC', 'SP_Peru', 'SP_Montane', 'SP_Ecuador',
                    'SP_mixture', 'SLC_Ecuador', 'SLC_Peru', 'SLC_SP_Peru',
                    'SLC_non_Andean', 'SLC_1', 'SLC_mixture', 'SLC_LA2792',
                    'SLC_Wva106', 'SLC_LA2135', 'SLL_vintage', 'SLL_fresh',
                    'SLL_processing_1', 'SLL_processing_2', 'SLL_1', 'mixture'],
         'group2': ['SC',
                    'SP_Peru_1', 'SP_Peru_2', 'SP_Peru_3', 'SP_Peru_4',
                    'SP_Peru_5', 'SP_Peru_6', 'SP_Peru_7', 'SP_Peru_8',
                    'SP_Peru_9', 'SP_Montane_1', 'SP_Montane_2',
                    'SP_Ecuador_1', 'SP_Ecuador_2', 'SP_Ecuador_3',
                    'SP_non_Andean', 'SP_mixture', 'SLC_SP_Peru',
                    'SLC_Ecuador_1', 'SLC_Ecuador_2', 'SLC_Ecuador_3', 'SLC_vintage',
                    'SLC_Peru_1', 'SLC_Peru_2', 'SLC_Peru_3',
                    'SLC_Colombia', 'SLC_Costa_Rica', 'SLC_Mesoamerica',
                    'SLC_Sinaloa', 'SLC_Asia', 'SLC_world', 'SLC_1',
                    'SLC_mixture', 'SLC_LA2792', 'SLC_LA2135', 'SLC_Wva106',
                    'SLL_Mesoamerica', 'SLL_vintage_1', 'SLL_early_breed',
                    'SLL_vintage_2', 'SLL_vintage/fresh', 'SLL_fresh_1',
                    'SLL_fresh_2', 'SLL_processing_1_1', 'SLL_processing_1_2',
                    'SLL_processing_1_3', 'SLL_processing_2', 'SLL_1',
                    'mixture']}

COMAV_SHAPE_ALLELES = {'fas': {'335': 'A', '466': 'B'},
                       'lc': {'395A': 'A', '395B': 'A', '533A': 'B',
                              '533B': 'B'},
                       'sun': {'wt': 'A', 'weird': None, 'nd': None,
                               'wt_pimpi': 'B', 'raro4': 'C', 'raro3': 'D'},
                       'sun_dup': {'A': 'A', None: None, 'B': 'A', 'C': 'B',
                                   'D': 'C'}}

INRA_SHAPE_ALLELES = {'fas': {'wt': 'A', 'fas': 'B', 'He': 'AB'},
                      'fw22': {'T': 'A', 'C': 'B'},
                      'fw32': {'G': 'A', 'A': 'B'},
                      'lc': {'TA': 'A', 'CG': 'B', 'YR': 'AB', 'NN': 'nan',
                             'CA': 'B'},
                      'ovate': {'wt': 'A', 'pst': 'B', '/': 'nan'}}


def _fix_acc_name(acc):
    if isinstance(acc, str):
        acc = acc.replace('+AC0', '')
        for long_, short in LONG_ACCESSION_NAMES.viewitems():
            acc = acc.replace(long_, short)
        for long_, short in SYMBOLS_TO_REPLACE_IN_ACCESSION_NAMES.viewitems():
            acc = acc.replace(long_, short)
    return acc


def read_physical_map():
    physical_map_fhand = StringIO()

    physical_map_header = ['chromosome', 'position']
    physical_map_fhand.write('\t'.join(physical_map_header))
    physical_map_fhand.write('\n')

    snp_id_col = 0
    chrom_col = 8
    physical_col = 9

    for index, line in enumerate(open(SOLCAP_DESIGN_FPATH)):
        if index < 9:
            continue
        line = line.strip().split('\t')
        snp_name = get_snp_solcap_id(line[snp_id_col])
        if line[chrom_col] == 'unknown':
            chrom = ''
        else:
            chrom = line[chrom_col].replace('SL2.40ch', '')
        if line[physical_col] == 'unknown':
            position = ''
        else:
            position = line[physical_col]
        map_line = '\t'.join([snp_name, chrom, position]) + '\n'
        physical_map_fhand.write(map_line)

    physical_map_fhand.seek(0)
    physical_map = load_markers_map(physical_map_fhand, sep='\t',
                                    molecule_col='chromosome',
                                    location_col='position')
    physical_map = MatrixWithMeta(physical_map.data.dropna(axis=0, how='any'),
                                  physical_map.meta)
    physical_map = sorted_marker_map(physical_map)

    return physical_map


def read_david_genotype_file(genotype_codec):
    first_acc_col = 7
    first_genetic_map_col = 5
    genotypes_fhand = StringIO()
    genetic_map_fhand = StringIO()
    for index, line in enumerate(open(DAVID_GENOTYPE_FPATH)):
        line = line.strip()
        if index == 0:
            continue
        elif index == 1:
            accessions = line.split('\t')[first_acc_col:]
            accessions = [_fix_acc_name(acc) for acc in accessions]
            genotypes_fhand.write('\t'.join(accessions))
            genotypes_fhand.write('\n')

            genetic_map_fhand.write('\t'.join(['chromosome', 'position']))
            genetic_map_fhand.write('\n')

        elif index > 2:
            items = line.split('\t')
            snp_genos = items[first_acc_col:]
            snp_genos = ['' if '-' in geno else geno for geno in snp_genos]
            snp_name = get_snp_solcap_id(items[0])
            genotypes_fhand.write(snp_name + '\t')
            genotypes_fhand.write('\t'.join(snp_genos))
            genotypes_fhand.write('\n')

            snp_genetic_map = items[first_genetic_map_col:
                                    first_acc_col]
            genetic_map_fhand.write(snp_name + '\t')
            map_line = '\t'.join(snp_genetic_map).replace(',', '.') + '\n'
            map_line = map_line.replace('Unknown', '')
            genetic_map_fhand.write(map_line)

    genetic_map_fhand.seek(0)
    genetic_map = load_markers_map(genetic_map_fhand, sep='\t',
                                   molecule_col='chromosome',
                                   location_col='position')
    genetic_map = MatrixWithMeta(genetic_map.data.dropna(axis=0, how='any'),
                                 genetic_map.meta)
    genetic_map = sorted_marker_map(genetic_map)

    genotypes_fhand.seek(0)
    genotypes = load_codominant_genetic_csv(genotypes_fhand, sep='\t',
                                            individuals_in_rows=False,
                                            genotype_codec=genotype_codec)

    # filter accs with too many missing data
    too_many_empty_genotypes = create_empty_genotypes_filter()
    genotypes = filter_matrix_cols(genotypes, too_many_empty_genotypes)

    return genotypes, genetic_map


def create_code_conversor(rosetta_code_table, from_field, to_field=None):
    def _convert_code(code):
        for sample, codes in rosetta_code_table.viewitems():
            if codes[from_field] == code:
                if to_field:
                    converted_code = codes[to_field]
                else:
                    converted_code = sample
                return converted_code
        msg = code + ' was not found in the table code field: ' + from_field
        raise ValueError(msg)
    return _convert_code


def read_solcap_sample_codes(preferred_key_columns=None):
    fhand = open(ACC_SYNS_FPATH)
    if preferred_key_columns is None:
        preferred_key_columns = ['our_comav_id']
    # every line is a row in our SolCap experiment
    rows = [row for row in DictReader(fhand, delimiter='\t')]

    # Now we give every sample in the experiment an id
    table = {}
    id_counts = Counter()
    for row in rows:
        code = row['preferred_id']
        if not code:
            continue
        preferred_id = None
        for column in preferred_key_columns:
            preferred_id = row[column]
            if preferred_id:
                break

        if not preferred_id:
            msg = 'No id in the columns given as preferred',
            msg += 'for row: ' + str(row)
            raise RuntimeError(msg)

        if id_counts[preferred_id]:
            id_ = preferred_id + '_' + str(id_counts[preferred_id] + 1)
        else:
            id_ = preferred_id
        id_counts[preferred_id] += 1
        table[id_] = row
    return table


def read_comav_genotype_file():
    rosetta_accs_table = read_solcap_sample_codes()

    convert_from_german_to_std = create_code_conversor(rosetta_accs_table,
                                                      from_field='german_code')

    # from customer strand to design strand
    #assert (guess_snp_coding(read_genotypes(open(COMAV_GENOTYPE_FPATH))) ==
    #                                                              'customer')
    customer_strand_genotypes = read_genotypes(open(COMAV_GENOTYPE_FPATH))
    design_genotypes = change_genotypes_strand(customer_strand_genotypes,
                                               input_strand=CUSTOMER_STRAND,
                                               output_strand=DESIGN_STRAND)
    design_genotypes_fhand = StringIO()
    write_genotypes_to_csv(design_genotypes, design_genotypes_fhand)
    design_genotypes_fhand.seek(0)

    # read the genotypes and modify the accession names
    genotypes = load_codominant_genetic_csv(design_genotypes_fhand,
                                            index_col=0, missing_tag='failed',
                                            individuals_in_rows=False)

    # filter accs with too many missing data
    too_many_empty_genotypes = create_empty_genotypes_filter()
    genotypes = filter_matrix_cols(genotypes, too_many_empty_genotypes)

    cols = genotypes.data.columns
    cols = [convert_from_german_to_std(acc_sample) for acc_sample in cols]
    genotypes.data.columns = cols

    # some comav accessions were not included in the Plos study
    plos_passports = load_dataframe_csv(open(PLOS_PASSPORTS_FPATH),
                                         index_col=0, sep='\t')
    code_conversor = create_code_conversor(rosetta_accs_table,
                                           from_field='preferred_id')
    rows = [code_conversor(acc_sample) for acc_sample in plos_passports.data.index]
    samples_in_plos_publication = rows
    samples_not_in_plos = set(genotypes.data.columns).difference(samples_in_plos_publication)
    samples_in_plos = set(genotypes.data.columns).difference(samples_not_in_plos)

    genotypes_in_plos = MatrixWithMeta(genotypes.data[list(samples_in_plos)],
                                       genotypes.meta)
    genotypes_not_in_plos = MatrixWithMeta(genotypes.data[list(samples_not_in_plos)],
                                           genotypes.meta)

    return genotypes_in_plos, genotypes_not_in_plos, rosetta_accs_table


def extract_small_ceras_genotypes():
    reader = DictReader(open(SMALL_CERAS_GENOTYPE_FPATH_ORIG), delimiter='\t')
    snp_fields = filter(lambda x: 'GType' in x, reader.fieldnames)
    accessions = map(lambda x: SMALL_CERAS_CODES[x.replace('.GType', '')],
                     snp_fields)

    ch_to_acc = dict(zip(snp_fields, accessions))
    genotypes_fhand = StringIO()
    sep = '\t'
    fields = ['maker']
    fields.extend(accessions)
    writer = DictWriter(genotypes_fhand, fieldnames=fields, delimiter=sep)
    writer.writeheader()
    for snp in reader:
        snp_name = get_snp_solcap_id(snp['Name'])
        snp_genotypes = {ch_to_acc[f]: '' if ('-' in snp[f] or 'NC' == snp[f]) else snp[f] for f in snp_fields}
        snp_genotypes['maker'] = snp_name
        writer.writerow(snp_genotypes)
    genotypes_fhand.seek(0)

    ab_genotypes = read_genotypes(genotypes_fhand)
    design_genotypes = change_ab_genotypes_strand(ab_genotypes,
                                                  input_strand=AB_CODING,
                                                  output_strand=DESIGN_STRAND)

    write_genotypes_to_csv(design_genotypes,
                           open(SMALL_CERAS_GENOTYPE_FPATH_DESIGN, 'w'))


def read_small_ceras_genotypes(genotype_codec):
    fhand = open(SMALL_CERAS_GENOTYPE_FPATH_DESIGN)
    genotypes = load_codominant_genetic_csv(fhand, individuals_in_rows=False,
                                            sep=',', index_col=0,
                                            genotype_codec=genotype_codec)
    return genotypes


def read_koenig_genotypes(genotype_codec, marker_names):
    fhand = open(KOENING_GENOTYPES_FPATH)
    genotypes = load_codominant_genetic_csv(fhand, individuals_in_rows=False,
                                            sep=',', index_col=0,
                                            genotype_codec=genotype_codec)
    # There are lots of markers not present in koenig and we don't wont
    # NANs later because their floats, so we fill the missing data with 0
    genos = genotypes.data.reindex(index=marker_names, fill_value=0)
    genotypes = MatrixWithMeta(genos, metadata=genotypes.meta)
    return genotypes


def extract_inra_genotypes():
    reader = DictReader(open(INRA_GENOTYPE_FPATH), delimiter='\t')
    snp_fields = filter(lambda x: 'CR' in x, reader.fieldnames)

    genotypes_fhand = StringIO()
    sep = '\t'
    fields = ['maker']
    fields.extend(snp_fields)

    writer = DictWriter(genotypes_fhand, fieldnames=fields, delimiter=sep)
    writer.writeheader()
    for snp in reader:
        snp_name = get_snp_solcap_id(snp['Name'])
        snp_genotypes = {f: '' if ('-' in snp[f] or 'NC' == snp[f]) else snp[f] for f in snp_fields}
        snp_genotypes['maker'] = snp_name
        writer.writerow(snp_genotypes)
    genotypes_fhand.seek(0)

    ab_genotypes = read_genotypes(genotypes_fhand)

    design_genotypes = change_ab_genotypes_strand(ab_genotypes,
                                                  input_strand=AB_CODING,
                                                  output_strand=DESIGN_STRAND)

    write_genotypes_to_csv(design_genotypes,
                           open(INRA_GENOTYPE_FPATH_DESIGN, 'w'))


def read_inra_genotypes(genotype_codec):
    genotypes_fhand = open(INRA_GENOTYPE_FPATH_DESIGN)
    genotypes = load_codominant_genetic_csv(genotypes_fhand,
                                            individuals_in_rows=False,
                                            sep=',', index_col=0,
                                            genotype_codec=genotype_codec)

    return genotypes


def extract_bauchet_pimpi_genotypes():
    fhand = open(BAUCHET_PIMPI_GENOTYPE_FPATH)
    genotypes = load_dataframe_csv(fhand, index_col=0, sep=';')
    genotypes = genotypes.data.dropna()
    snp_ids = [get_snp_solcap_id(snp_id) for snp_id in genotypes.columns]
    genotypes.columns = snp_ids
    genotypes.replace('NC', '', inplace=True)

    genotypes_dict = genotypes.to_dict()

    genotypes_fhand = StringIO()
    fields = ['marker']
    fields.extend(genotypes.index)
    writer = DictWriter(genotypes_fhand, fieldnames=fields, delimiter='\t')
    writer.writeheader()
    for snp, geno in genotypes_dict.viewitems():
        genotype = geno
        genotype['marker'] = snp
        writer.writerow(genotype)
    genotypes_fhand.seek(0)

    class TabDialect(Dialect):
        delimiter = '\t'
        quoting = QUOTE_MINIMAL
        quotechar = '"'
        doublequote = True
        skipinitialspace = False
        lineterminator = '\r\n'
        quoting = QUOTE_MINIMAL

    genotypes = read_genotypes(genotypes_fhand, TabDialect)

    design_genotypes = change_ab_genotypes_strand(genotypes,
                                                  input_strand=AB_CODING,
                                                  output_strand=DESIGN_STRAND)

    write_genotypes_to_csv(design_genotypes,
                           open(BAUCHET_PIMPI_GENOTYPE_FPATH_DESIGN, 'w'))


def read_bauchet_pimpi_genotypes(genotype_codec):
    genotypes_fhand = open(BAUCHET_PIMPI_GENOTYPE_FPATH_DESIGN)
    genotypes = load_codominant_genetic_csv(genotypes_fhand,
                                            individuals_in_rows=True,
                                            index_col=0, sep=',', ploidy=2,
                                            missing_tag='',
                                            genotype_codec=genotype_codec)
    return genotypes


def extract_david2_genotypes():
    fhand = open(DAVID2_GENOTYPE_FPATH)

    class TabDialect(Dialect):
        delimiter = '\t'
        quoting = QUOTE_MINIMAL
        quotechar = '"'
        doublequote = True
        skipinitialspace = False
        lineterminator = '\r\n'
        quoting = QUOTE_MINIMAL

    genotypes = read_genotypes(fhand, TabDialect)
    write_genotypes_to_csv(genotypes,
                           open(DAVID2_GENOTYPE_FPATH_DESING, 'w'))


def read_david2_genotypes(genotype_codec):
    genotypes = load_codominant_genetic_csv(DAVID2_GENOTYPE_FPATH_DESING,
                                            individuals_in_rows=False,
                                            index_col=0,
                                            genotype_codec=genotype_codec)
    return genotypes


def merge_genotypes(plos_genotypes, new_comav_genotypes, david_genotypes,
                    small_ceras_genotypes, inra_genotypes, bauchet_genotypes,
                    david2_genotypes, koening_genotypes):
    if plos_genotypes.meta[INDIVIDUALS_IN_ROWS]:
        raise ValueError('merge genotypes only prepared for indis in cols')

    n_plos_accs = len(plos_genotypes.data.columns)
    n_new_comav_accs = len(new_comav_genotypes.data.columns)
    n_david_accs = len(david_genotypes.data.columns)
    n_inra_accs = len(inra_genotypes.data.columns)
    n_small_ceras_accs = len(small_ceras_genotypes.data.columns)
    n_bauchet_accs = len(bauchet_genotypes.data.columns)
    n_david2_accs = len(david2_genotypes.data.columns)
    n_koening_accs = len(koening_genotypes.data.columns)

    genotypes = hstack_dframes(plos_genotypes.data, david_genotypes.data,
                               lsuffix='_comav')
    genotypes = hstack_dframes(genotypes, new_comav_genotypes.data,
                               rsuffix='_comav')
    genotypes = hstack_dframes(genotypes, small_ceras_genotypes.data,
                               rsuffix='_comav')
    genotypes = hstack_dframes(genotypes, inra_genotypes.data, rsuffix='_inra')
    genotypes = hstack_dframes(genotypes, bauchet_genotypes.data,
                               rsuffix='_inra')
    genotypes = hstack_dframes(genotypes, david2_genotypes.data,
                               rsuffix='_david')
    genotypes = hstack_dframes(genotypes, koening_genotypes.data,
                               rsuffix='_koenig')

    provenances = n_plos_accs * ['comav'] + n_david_accs * ['david'] + \
                  n_new_comav_accs * ['comav_new'] + \
                  n_small_ceras_accs * ['comav_new'] + \
                  n_inra_accs * ['inra'] + n_bauchet_accs * ['inra'] + \
                  n_david2_accs * ['david_new'] + n_koening_accs * ['koenig']
    provenances = DataFrame({'provenance':
                              Series(provenances, index=genotypes.columns)})
    provenances = MatrixWithMeta(provenances)
    genotypes = MatrixWithMeta(genotypes, plos_genotypes.meta)

    return genotypes, provenances


def do_multivariate(genotypes, acc_info, marker_col, color_col,
                    color_definitions, marker_definitions,
                    labels=False):
    pca_dir = pjoin(RESULTS_DIR, 'pca')
    if not exists(pca_dir):
        os.mkdir(pca_dir)

    smartpca_params = {'num_chromosomes': 12}
    pca = do_smartpca(genotypes, params=smartpca_params, max_princomps=10)
    pca_plot_fhand12 = open(os.path.join(pca_dir, 'snp_pca_1_2.svg'), 'w')
    pca_plot_fhand13 = open(os.path.join(pca_dir, 'snp_pca_1_3.svg'), 'w')

    if marker_col == color_col:
        classification = [acc_info.data[marker_col],]
    else:
        classification = [acc_info.data[marker_col],
                          acc_info.data[color_col]]
    pca_projections = dframe_append_columns(pca['projections'].data,
                                            classification)
    fhand = open(pjoin(pca_dir, 'projections.csv'), 'w')
    pca_projections.to_csv(fhand)
    print 'pca written'
    print 'axis\t variance (%)\n', pca['var_percentages']
    plot_groups = get_scatter_groups_from_matrix(pca_projections,
                                                 x_col='pca-01',
                                                 y_col='pca-02',
                                                 marker_col=marker_col,
                                                 color_col=color_col,
                                    color_mapper=color_definitions[color_col],
                                 marker_mapper=marker_definitions[marker_col])
    scatter_groups(plot_groups, pca_plot_fhand12, labels=labels,
                   alpha=0.7, marker_size=150)
    plot_groups = get_scatter_groups_from_matrix(pca_projections,
                                                 x_col='pca-01',
                                                 y_col='pca-03',
                                                 marker_col=marker_col,
                                                 color_col=color_col,
                                    color_mapper=color_definitions[color_col],
                                  marker_mapper=marker_definitions[marker_col])
    scatter_groups(plot_groups, pca_plot_fhand13, labels=labels, alpha=0.9,
                   marker_size=230)


def filter_incongruent_snps(physical_map, genetic_map):
    '''It removes the SNPs that are not in the same chromosome in the physical
     and genetic map
     '''
    gen_chrom_col = genetic_map.meta[MOLECULE_COL]
    phy_chrom_col = physical_map.meta[MOLECULE_COL]

    filter_snp = []
    for snp in genetic_map.data.index:
        gen_chrom = genetic_map.data.loc[snp, gen_chrom_col]
        if snp in physical_map.data.index:
            phy_chrom = physical_map.data.loc[snp, phy_chrom_col]
        else:
            continue
        if gen_chrom != phy_chrom:
            filter_snp.append(snp)

    physical_map_data = physical_map.data.drop(filter_snp)
    genetic_map_data = genetic_map.data.drop(filter_snp)

    physical_map = MatrixWithMeta(physical_map_data, physical_map.meta)
    genetic_map = MatrixWithMeta(genetic_map_data, genetic_map.meta)

    return physical_map, genetic_map


def create_empty_genotypes_filter(empty_percent_threshold=10):
    empty_per_thresh = empty_percent_threshold / 100

    def _too_many_empty_genotypes(series):
        if len(list(filter(bool, series))) / len(series) > empty_per_thresh:
            return True
        else:
            return False
    return _too_many_empty_genotypes


def create_empty_genotypes_in_accs_filter(accs):

    def _empty_genotypes(series):
        if len(list(filter(bool, series.reindex(accs)))) > 1:
            return True
        else:
            return False
    return _empty_genotypes


def create_monomorphic_marker_check(genotype_codec, maf_threshold=0.95):
    def _monomorphic_marker_check(marker_genotypes):
        genotypes = [genotype_codec.decode_to_ints(genotype) for genotype in marker_genotypes if genotype != MISSING_INT_GENOTYPE and genotype]
        counts = Counter(allele for geno in genotypes for allele in geno)
        assert 0 < len(counts) <= 2
        counts = list(counts.viewvalues())
        if max(counts) / sum(counts) < maf_threshold:
            return True
        else:
            return False
    return _monomorphic_marker_check


def filter_spreaded_markers_in_markers(genotypes, physical_map,
                                       distance=0.1):
    if genotypes.meta[INDIVIDUALS_IN_ROWS]:
        markers = list(genotypes.data.columns)
    else:
        markers = list(genotypes.data.index)

    location_col = physical_map.meta[LOCATION_COL]
    molecule_col = physical_map.meta[MOLECULE_COL]

    map_dframe = physical_map.data
    selected_markers = []
    prev_location = None
    prev_molecule = None
    for marker in markers:
        try:
            location = map_dframe.get_value(marker, location_col)
            molecule = map_dframe.get_value(marker, molecule_col)
        except KeyError:
            continue
        if (prev_molecule is None or
            prev_molecule != molecule or
            (prev_molecule == molecule and location - prev_location > distance)):
            selected_markers.append(marker)
        prev_location = location
        prev_molecule = molecule

    if genotypes.meta[INDIVIDUALS_IN_ROWS]:
        dframe = genotypes.data[selected_markers]
    else:
        dframe = genotypes.data.ix[selected_markers]

    genotypes = MatrixWithMeta(dframe, genotypes.meta)
    return genotypes


def filter_monomorphic_markers(genotypes, codec, maf_threshold=0.95):
    monomorphic_marker = create_monomorphic_marker_check(genotype_codec=codec,
                                                   maf_threshold=maf_threshold)
    if genotypes.meta[INDIVIDUALS_IN_ROWS]:
        genotypes = filter_matrix_cols(genotypes, monomorphic_marker)
    else:
        genotypes = filter_matrix_rows(genotypes, monomorphic_marker)
    return genotypes


def filter_missing_koenig_genotypes(genotypes, acc_info):
    koenig_info = filter_matrix_rows(acc_info, {'provenance': 'koenig'})
    accs = koenig_info.data.index

    empty_genotypes = create_empty_genotypes_in_accs_filter(accs)

    if genotypes.meta[INDIVIDUALS_IN_ROWS]:
        genotypes = filter_matrix_cols(genotypes, empty_genotypes)
    else:
        genotypes = filter_matrix_rows(genotypes, empty_genotypes)
    return genotypes


def filter_missing_data_markers(genotypes):
    too_many_empty_genotypes = create_empty_genotypes_filter()

    if genotypes.meta[INDIVIDUALS_IN_ROWS]:
        genotypes = filter_matrix_cols(genotypes, too_many_empty_genotypes)
    else:
        genotypes = filter_matrix_rows(genotypes, too_many_empty_genotypes)
    return genotypes


def filter_samples(genotypes, acc_info, condition, reverse=False):
    filtered_sample_matrix = filter_matrix_rows(acc_info, condition)

    filtered_samples = list(filtered_sample_matrix.data.index)
    if genotypes.meta[INDIVIDUALS_IN_ROWS]:
        genotypes = filter_matrix_rows(genotypes, filtered_samples,
                                       reverse=reverse)
    else:
        genotypes = filter_matrix_cols(genotypes, filtered_samples,
                                       reverse=reverse)
    return genotypes


def filter_matrix(matrix, acc_info, condition, reverse=False):
    filtered_sample_matrix = filter_matrix_rows(acc_info, condition)
    filtered_samples = list(filtered_sample_matrix.data.index)

    matrix = filter_matrix_rows(matrix, filtered_samples,
                                       reverse=reverse)
    matrix = filter_matrix_cols(matrix, filtered_samples,
                                       reverse=reverse)
    return matrix


def read_comav_sample_passports(comav_acc_syns):
    passports = load_dataframe_csv(open(COMAV_PASSPORTS_FPATH), sep='\t')
    tpassports = passports.data.T
    passports = MatrixWithMeta(tpassports)

    # we need the passports for the samples genotyped in our SOLcap
    # experiment. Some accessions were genotyped twice, so we have to
    # duplicate those passports and give them the corresponding sample codes
    pass_samples = {}
    for sample, syns in comav_acc_syns.viewitems():
        sample_id = syns['german_code']
        if not sample_id:
            continue
        passports_for_sample = filter_matrix_cols(passports,
                                            criteria={'sample_ID': sample_id})

        n_passports = len(passports_for_sample.data.columns)
        if not n_passports:
            print 'no pass for ' + sample_id
        elif n_passports == 1:
            col_name = list(passports_for_sample.data.columns)[0]
            sample_pass = passports_for_sample.data[col_name]
        else:
            raise RuntimeError('Several passports for one sample: ' + sample)
        pass_samples[sample] = sample_pass

    passports = DataFrame(pass_samples)
    tpassports = passports.T

    columns_to_del = ['ID', 'sample_ID', 'group2',
                      'peticion', 'anotaciones morfologia maria jose',
                      'var_ini', 'sp_fin', 'var_fin', 'sp_ini', 'entregas',
                      'Nombre local', 'fail_kk', 'poblacion',
                      'Especie_nuestra', 'Obs_descripcion',
                      'Obs_internas banco', 'uso', '28_Observaciones',
                      'GENERO', 'sample_ID', 'clasificacion_morfologica_paper',
                      'obs_origen', 'clasificacion_morfologica', 'Unnamed: 7',
                      'clasificacion_banco_original']

    tpassports = tpassports.drop(columns_to_del, axis=1)
    passports = MatrixWithMeta(tpassports)

    return passports


def read_david_sample_passports(david_acc_syns):
    passports = load_dataframe_csv(open(DAVID_PASSPORTS_FPATH), sep='\t',
                                   index_col=1, header=1)
    accs = passports.data.index
    passports_index = [_fix_acc_name(acc) for acc in accs]
    for acc in passports_index:
        assert acc in david_acc_syns
    passports.data.index = passports_index

    columns_to_del = ['Donor Source', 'missclassified', 'Comments_TGRC',
                      'Photos', 'Male', 'Female',
                      ' IBC', 'Other e.g. BC', 'Donor Institution', 'group1',
                      'Name of Donor', 'Synonyms', 'SolCAP T+AC0-number']
    passports_data = passports.data.drop(columns_to_del, axis=1)

    return MatrixWithMeta(passports_data)


def read_classification(classification_fhand):
    passports = load_dataframe_csv(classification_fhand, sep='\t',
                                   index_col=0, header=0)
    accs = passports.data.index
    passports_index = [_fix_acc_name(acc) for acc in accs]
    passports.data.index = passports_index
    return passports


def read_esther_shape_genotypes():
    fhand = open(ESTHER_SHAPE_GENOTYPES_FPATH)
    esther_shape_genotypes = load_dataframe_csv(fhand,
                                                index_col=0, header=0,
                                                sep='\t')
    esther_genotype_accs = esther_shape_genotypes.data.index
    shape_geno = esther_genotype_accs

    cols = ['sun_dup' if col == 'sun' else col for col in esther_shape_genotypes.data.columns]
    esther_shape_genotypes.data.columns = cols

    fhand = ESTHER_SHAPE_GENOTYPES_SYNONYMS_FPATH
    esther_synonims = load_dataframe_csv(fhand, index_col=0, header=0,
                                                sep='\t')
    esther_synoms_accs = esther_synonims.data.index

    corrected_accs = []

    for acc in esther_genotype_accs:
        if acc in esther_synoms_accs:
            corrected_accs.append(esther_synonims.data.ix[acc]['our_id'])
        else:
            corrected_accs.append(acc)
    esther_shape_genotypes.data.index = corrected_accs

    def convert_esther_geno_coding(genotype):
        if not isinstance(genotype, str):
            genotype = '%.0f' % genotype
        if genotype == '1':
            return ('B', 'B')
        elif genotype == '2':
            return ('A', 'B')
        elif genotype == '3':
            return ('A', 'A')
        else:
            return (float('nan'), float('nan'))

    shape_genotypes = {}

    for gene, gene_genotypes in esther_shape_genotypes.data.iteritems():
        gene_genotypes = gene_genotypes.dropna()
        accs = gene_genotypes.index
        genos = [convert_esther_geno_coding(geno) for geno in gene_genotypes]
        df = DataFrame(genos, columns=['allele1', 'allele2'],
                       index=accs).dropna()
        shape_genotypes[gene] = MatrixWithMeta(df)

    return shape_genotypes


def read_esther_2_shape_genotypes():
    fhand = open(ESTHER_2_SHAPE_GENOTYPES_FPATH)
    esther2_shape_genotypes = load_dataframe_csv(fhand, index_col=1, sep='\t')
    shape_genes = ['FW2.2', 'FW3.2', 'FW11.3', 'LC', 'FAS',
                   'OVATE', 'SUN']

    esther2_shape_genotypes = esther2_shape_genotypes.data[shape_genes]
    accs = [_fix_acc_name(acc) for acc in esther2_shape_genotypes.index]
    esther2_shape_genotypes.index = accs

    def convert_esther_geno_coding(genotype):
        if not isinstance(genotype, str):
            genotype = '%.0f' % genotype
        if genotype == '1':
            return ('B', 'B')
        elif genotype == '2':
            return ('A', 'B')
        elif genotype == '3':
            return ('A', 'A')
        else:
            return (float('nan'), float('nan'))

    esther2_shape_genotypes = esther2_shape_genotypes.applymap(convert_esther_geno_coding)

    shape_genes = ['fw22', 'fw32', 'fw113', 'lc', 'fas', 'ovate', 'sun_dup']
    esther2_shape_genotypes.columns = shape_genes

    shape_genotypes = {}
    for gene, gene_genotypes in esther2_shape_genotypes.iteritems():
        gene_genotypes = gene_genotypes.dropna()
        df = DataFrame(list(gene_genotypes), columns=['allele1', 'allele2'],
                       index=gene_genotypes.index).dropna()
        shape_genotypes[gene] = MatrixWithMeta(df)

    return shape_genotypes


def read_inra_shape_files():
    fhand = open(INRA_SHAPE_GENOTYPES_FPATH)
    inra_shape_genotypes = load_dataframe_csv(fhand, index_col=0, sep=',')
    inra_shape_genotypes.data.replace('NA', numpy.float('nan'), inplace=True)

    inra_shape_genotypes.data['lc'] = inra_shape_genotypes.data['lc1'] + \
                                      inra_shape_genotypes.data['lc2']

    shape_genes = ['lc', 'fw22', 'fw32', 'fas']
    inra_shape_genotypes = inra_shape_genotypes.data[shape_genes]

    fhand = open(INRA_OVATE_GENOTYPES_FPATH)
    inra_ovate_genotypes = load_dataframe_csv(fhand, index_col=0, sep='\t')
    inra_shape_genotypes['ovate'] = inra_ovate_genotypes.data

    def convert_inra_geno_coding(gene, genotype):
        genotype = INRA_SHAPE_ALLELES[gene][genotype]
        if len(genotype) == 1 and genotype == 'B':
            return ('B', 'B')
        elif len(genotype) == 1 and genotype == 'A':
            return ('A', 'A')
        elif len(genotype) == 2:
            return ('A', 'B')
        elif genotype == 'nan':
            return (numpy.float('nan'), numpy.float('nan'))
        else:
            raise ValueError('Genotype ' + genotype + ' is not in the list')

    shape_genotypes = {}
    for gene, gene_genotypes in inra_shape_genotypes.iteritems():
        gene_genotypes = gene_genotypes.dropna()
        accs = gene_genotypes.index
        genos = [convert_inra_geno_coding(gene, geno) for geno in gene_genotypes]
        df = DataFrame(genos, columns=['allele1', 'allele2'],
                       index=accs).dropna()
        shape_genotypes[gene] = MatrixWithMeta(df)

    return shape_genotypes


def read_comav_shape_files():

    fw32_genotypes = load_dataframe_csv(FW32_GENOTYPES_FPATH, sep=',',
                                        index_col=0)
    fw113_genotypes = load_dataframe_csv(FW113_GENOTYPES_FPATH, sep=',',
                                         index_col=0)
    lc_haplotypes = load_dataframe_csv(LC_HAPLOTYPES_FPATH, sep=',',
                                       index_col=0)
    fas_genotypes = load_dataframe_csv(FAS_GENOTYPES_FPATH, sep=',',
                                       index_col=0)
    fw22_genotypes = load_dataframe_csv(FW22_GENOTYPES_FPATH, sep=',',
                                       index_col=0)
    ovate_genotypes = load_dataframe_csv(OVATE_GENOTYPES_FPATH, sep='\t',
                                       index_col=0)
    sun_genotypes = load_dataframe_csv(SUN_GENOTYPES_FPATH, sep='\t',
                                       index_col=0)

    rosetta_accs_table = read_solcap_sample_codes()

    convert_acc_id_to_std1 = create_code_conversor(rosetta_accs_table,
                                                  from_field='preferred_id')
    convert_acc_id_to_std2 = create_code_conversor(rosetta_accs_table,
                                                  from_field='our_comav_id')

    def convert_acc_id_to_std(acc_id):
        try:
            new_id = convert_acc_id_to_std1(acc_id)
        except ValueError:
            new_id = convert_acc_id_to_std2(acc_id)
        return new_id

    def codify_alleles(shape_gene, genotypes):
        if shape_gene not in COMAV_SHAPE_ALLELES.viewkeys():
            return genotypes

        def code(value):
            return COMAV_SHAPE_ALLELES[shape_gene][value]

        genotypes = genotypes.data
        if any(genotypes.dtypes == numpy.int64):
            genotypes = genotypes.astype(str)
        genotypes = genotypes.applymap(code).dropna()
        return MatrixWithMeta(genotypes)

    accs = fw32_genotypes.data.index
    fw32_index = [convert_acc_id_to_std(acc) for acc in accs]
    fw32_genotypes.data.index = fw32_index

    accs = fw113_genotypes.data.index
    fw113_index = [convert_acc_id_to_std(acc) for acc in accs]
    fw113_genotypes.data.index = fw113_index

    accs = fw22_genotypes.data.index
    fw22_index = [convert_acc_id_to_std(acc) for acc in accs]
    fw22_genotypes.data.index = fw22_index

    accs = fas_genotypes.data.index
    fas_index = [convert_acc_id_to_std(acc) for acc in accs]
    fas_genotypes.data.index = fas_index
    fas_genotypes = codify_alleles('fas', fas_genotypes)

    accs = lc_haplotypes.data.index
    lc_index = [convert_acc_id_to_std(acc) for acc in accs]
    lc_haplotypes.data.index = lc_index
    lc_haplotypes = codify_alleles('lc', lc_haplotypes)

    accs = ovate_genotypes.data.index
    ovate_index = [convert_acc_id_to_std(acc) for acc in accs]
    ovate_genotypes.data.index = ovate_index

    accs = sun_genotypes.data.index
    sun_index = [convert_acc_id_to_std(acc) for acc in accs]
    sun_genotypes.data.index = sun_index
    sun_genotypes = codify_alleles('sun', sun_genotypes)

    sun_dup_genotypes = codify_alleles('sun_dup', sun_genotypes)

    return {'fw113': fw113_genotypes, 'lc': lc_haplotypes,
            'fas': fas_genotypes, 'fw32': fw32_genotypes,
            'fw22': fw22_genotypes, 'ovate': ovate_genotypes,
            'sun': sun_genotypes, 'sun_dup': sun_dup_genotypes}


def _check_sample_name(shape_genotypes, acc_syns):

    for gene in shape_genotypes:
        sample_names = []
        for sample in shape_genotypes[gene].data.index:
            if sample in acc_syns:
                sample_names.append(acc_syns[sample])
            else:
                sample_names.append(sample)
        shape_genotypes[gene].data.index = sample_names

    return shape_genotypes


def merge_shape_genotypes(acc_syns):
    esther_shape_genotypes = read_esther_shape_genotypes()
    esther_shape_genotypes = _check_sample_name(esther_shape_genotypes,
                                                acc_syns)

    comav_shape_genotypes = read_comav_shape_files()
    comav_shape_genotypes = _check_sample_name(comav_shape_genotypes, acc_syns)

    print 'comav sun accs', comav_shape_genotypes['sun'].data.index
    print 'alleles', set(comav_shape_genotypes['sun'].data['allele1'])


    inra_shape_genotypes = read_inra_shape_files()
    inra_shape_genotypes = _check_sample_name(inra_shape_genotypes, acc_syns)

    esther2_shape_genotyes = read_esther_2_shape_genotypes()
    esther2_shape_genotyes = _check_sample_name(esther2_shape_genotyes,
                                                acc_syns)

    genes = list(esther_shape_genotypes.viewkeys())
    genes.append('sun')
    merged_shape_genotypes = {}
    for gene in genes:
        if gene in esther_shape_genotypes.viewkeys():
            esther_genos = esther_shape_genotypes[gene].data
            merged_genos = esther_genos
        else:
            merged_genos = None

        if gene in comav_shape_genotypes.viewkeys():
            comav_genos = comav_shape_genotypes[gene].data

            if merged_genos is not None:
                merged_accs = merged_genos.index
                comav_acc_to_remove = set()
                merged_acc_to_remove = set()
                for acc in merged_accs:
                    if acc in comav_genos.index:
                        if all(merged_genos.ix[acc] != comav_genos.ix[acc]):
                            print 'Shape genotype for', acc, 'in gene', gene, \
                                  'do not match between Esther and COMAV'
                            merged_acc_to_remove.add(acc)
                            comav_acc_to_remove.add(acc)
                        else:
                            merged_acc_to_remove.add(acc)

                comav_genos = comav_genos.drop(comav_acc_to_remove)
                merged_genos = merged_genos.drop(merged_acc_to_remove)
                merged_genos = vstack_dframes(merged_genos, comav_genos)
            else:
                merged_genos = comav_genos

        if gene in inra_shape_genotypes:
            inra_genos = inra_shape_genotypes[gene].data
            merged_accs = merged_genos.index
            inra_acc_to_remove = set()
            merged_acc_to_remove = set()
            for acc in merged_accs:
                if acc in inra_genos.index:
                    if all(merged_genos.ix[acc] != inra_genos.ix[acc]):
                        print 'Shape genotype for', acc, 'in gene', gene, \
                              'do not match between Esther or COMAV and INRA'
                        inra_acc_to_remove.add(acc)
                        merged_acc_to_remove.add(acc)
                    else:
                        inra_acc_to_remove.add(acc)
            inra_genos = inra_genos.drop(inra_acc_to_remove)
            merged_genos = merged_genos.drop(merged_acc_to_remove)
            merged_genos = vstack_dframes(merged_genos, inra_genos)

        if gene in esther2_shape_genotyes:
            esther2_genos = esther2_shape_genotyes[gene].data
            merged_accs = merged_genos.index
            esther2_acc_to_remove = set()
            merged_acc_to_remove = set()
            for acc in merged_accs:
                if acc in esther2_genos.index:
                    if all(merged_genos.ix[acc] != esther2_genos.ix[acc]):
                        print 'Shape genotype for', acc, 'in gene', gene, \
                       'do not match between Esther, COMAV or INRA and Eshter2'
                        esther2_acc_to_remove.add(acc)
                        merged_acc_to_remove.add(acc)
                    else:
                        esther2_acc_to_remove.add(acc)
            esther2_genos = esther2_genos.drop(esther2_acc_to_remove)
            merged_genos = merged_genos.drop(merged_acc_to_remove)
            merged_genos = vstack_dframes(merged_genos, esther2_genos)

        merged_shape_genotypes[gene] = MatrixWithMeta(merged_genos)

    return merged_shape_genotypes


def create_syns_replace(syns_dict):
    def _syns_replace(item):
        if item in syns_dict.keys():
            return syns_dict[item]
        else:
            return item
    return _syns_replace


def merge_passports(comav_passports, david_passports):
    column_syns = {'Pais': 'Origin_Country',
                   'Origin Country': 'Origin_Country',
                   'Scientific Name': 'final_species',
                   'Lugar de colecta': 'Collection_site',
                   'latitud': 'Latitude', 'longitud': 'Longitude',
                   'altitud': 'Elevation'}

    comav_passports.data.rename(columns=column_syns, inplace=True)
    david_passports.data.rename(columns=column_syns, inplace=True)

    comav_passports.data.rename(index=REPEATED_SAMPLES_COMAV, inplace=True)

    inra_passports = load_dataframe_csv(open(INRA_PASSPORTS_FPATH), sep='\t',
                                        index_col=0, header=0)
    column_syns = {'Species': 'final_species'}
    inra_passports.data.rename(columns=column_syns, inplace=True)

    small_ceras_pass_fhand = open(SMALL_CERAS_PASSPORTS_FPATH)
    small_ceras_passports = load_dataframe_csv(small_ceras_pass_fhand,
                                               sep='\t', index_col=0, header=0)
    column_syns = {'species': 'final_species', 'latitud': 'Latitude',
                   'longitud': 'Longitude', 'altitud': 'Elevation',
                   'Pais': 'Origin_Country'}
    small_ceras_passports.data.rename(columns=column_syns, inplace=True)

    bauchet_passports = load_dataframe_csv(open(BAUCHET_PIMPI_PASSPORTS_FPATH),
                                         header=0, sep=';', index_col=0)
    bauchet_passports.data.replace('#N/A', float('nan'), inplace=True)
    column_syns = {'species': 'final_species', 'Altitude': 'Elevation',
                   'Country': 'Origin_Country', 'Coll.site': 'Collection_site',
                   'passport_taxonomy': 'final_species',
                   'Province/Department': 'Origin State/Province'}
    bauchet_passports.data.rename(columns=column_syns, inplace=True)
    bauchet_passports.data.rename(index=REPEATED_SAMPLES_BAUCHET, inplace=True)

    columns_to_delete = ['Assoc. biota', 'Coll.no.', 'Coll.notes', 'Habitat',
                         'Elevation_source', 'Coordenates_source', 'Pop.size',
                         'Sample_size', 'Site_details', 'Veg. type']

    bauchet_passports_data = bauchet_passports.data.drop(columns_to_delete,
                                                         axis=1)
    bauchet_passports = MatrixWithMeta(bauchet_passports_data)

    david2_passports = load_dataframe_csv(open(DAVID2_PASSPORTS_FPATH),
                                          index_col=0, header=0, sep='\t')
    david2_passports.data.rename(index=REPEATED_SAMPLES_DAVID2, inplace=True)

    koenig_passports = load_dataframe_csv(open(KOENIG_PASSPORTS_FPATH),
                                          sep='\t', index_col=0, header=0)
    column_syns = {'Species': 'final_species'}
    koenig_passports.data.rename(columns=column_syns, inplace=True)
    koenig_passports.data.rename(index=REPEATED_SAMPLES_KOENIG, inplace=True)

    passports = vstack_dframes(comav_passports.data, david_passports.data)
    passports = vstack_dframes(passports, inra_passports.data)
    passports = vstack_dframes(passports, small_ceras_passports.data)
    passports = vstack_dframes(passports, bauchet_passports.data)
    passports = vstack_dframes(passports, david2_passports.data)
    passports = vstack_dframes(passports, koenig_passports.data)

    country_syns = create_syns_replace(PASSPORT_COUNTRY_SYNS)
    passports['Origin_Country'] = passports['Origin_Country'].map(country_syns)

    species_syns = create_syns_replace(PASSPORT_SPECIES_SYNS)
    passports['final_species'] = passports['final_species'].map(species_syns)

    passports = passports.fillna('None')

    return MatrixWithMeta(passports)


def merge_filter_data(filter_missing_data=True, filter_monomorphic=True,
                      filter_close_data=True,
                      filter_missing_koenig=False):

    interpolated_recomb_map = load_markers_map(open(INTERPOLED_RECOM_MAP),
                                               molecule_col='chromosome',
                                               location_col='genetic_position',
                                               index_col=0, sep='\t')

    physical_map = load_markers_map(open(INTERPOLED_RECOM_MAP),
                                    molecule_col='chromosome',
                                    location_col='position',
                                    index_col=0, sep='\t')

    plos_genotypes, new_comav_genotypes, comav_acc_syns = read_comav_genotype_file()
    comav_sample_passports = read_comav_sample_passports(comav_acc_syns)
    meta = plos_genotypes.meta
    genotype_codec = GenotypeCodec(ploidy=meta[PLOIDY],
                                   bits_per_allele=meta[BITS_PER_ALLELE],
                                   alleles_coding=meta[ALLELE_CODING])

    koening_genotypes = read_koenig_genotypes(genotype_codec,
                                        marker_names=plos_genotypes.data.index)

    extract_small_ceras_genotypes()
    small_ceras_genotypes = read_small_ceras_genotypes(genotype_codec)

    extract_inra_genotypes()
    inra_genotypes = read_inra_genotypes(genotype_codec)

    extract_bauchet_pimpi_genotypes()
    pimpi_genotypes = read_bauchet_pimpi_genotypes(genotype_codec)

    extract_david2_genotypes()
    david2_genotypes = read_david2_genotypes(genotype_codec)

    david_genotypes, genetic_map = read_david_genotype_file(genotype_codec)
    david_acc_syns = david_genotypes.data.columns
    david_sample_passports = read_david_sample_passports(david_acc_syns)

    genotypes, provenances = merge_genotypes(plos_genotypes,
                                             new_comav_genotypes,
                                             david_genotypes,
                                             small_ceras_genotypes,
                                             inra_genotypes, pimpi_genotypes,
                                             david2_genotypes,
                                             koening_genotypes)
    # write_merged_non_filtered_genotypes
    fhand = open(MERGED_SOLCAP_DESIGN_NON_FILTERED_GENOTYPES, 'w')
    write_genetic_csv(fhand, genotypes)

    genotypes = transpose_genotypes(genotypes)

    passports = merge_passports(comav_sample_passports, david_sample_passports)

    acc_info = provenances.data.join(passports.data, how='inner')

    acc_info = MatrixWithMeta(acc_info)

    print 'No marker filtering'
    print genotypes
    print '-' * 50

    if filter_missing_koenig:
        print 'Missing koenig genotypes filtering'
        genotypes = filter_missing_koenig_genotypes(genotypes, acc_info)
        print genotypes

        print '-' * 50

    if filter_missing_data:
        # filter markers with too much missing data
        print 'Missing data marker filtering'
        genotypes = filter_missing_data_markers(genotypes)
        print genotypes
        print '-' * 50

    if filter_monomorphic:
        print 'Monomorphic marker filtering'
        genotypes = filter_monomorphic_markers(genotypes, genotype_codec,
                                               maf_threshold=0.95)
        print genotypes
        print '-' * 50

    markers = genotypes.data.columns
    sorted_markers = interpolated_recomb_map.data.index
    right_markers = [mar for mar in sorted_markers if mar in markers]
    genotypes = MatrixWithMeta(genotypes.data[right_markers],
                                   genotypes.meta)
    if filter_close_data:

        print 'Markers too close filtered'
        genotypes = filter_spreaded_markers_in_markers(genotypes,
                                                       interpolated_recomb_map)
        print genotypes
        print '-' * 50

    if not exists(RESULTS_DIR):
        os.mkdir(RESULTS_DIR)

    genotypes.save(open(SAVED_MERGED_FILTERED_GENOTYPES_FPATH, 'w'))
    physical_map.save(open(SAVED_PHYSICAL_MAP_FPATH, 'w'))
    interpolated_recomb_map.save(open(SAVED_GENETIC_MAP_FPATH, 'w'))
    acc_info.save(open(SAVED_MERGED_PASSPORTS_FPATH, 'w'))


def load_data_matrices():
    'Loads the genotype, physical map, genetic map and passports matrices'
    genotypes = load_matrix(open(SAVED_MERGED_FILTERED_GENOTYPES_FPATH))
    physical_map = load_matrix(open(SAVED_PHYSICAL_MAP_FPATH))
    genetic_map = load_matrix(open(SAVED_GENETIC_MAP_FPATH))
    passports = load_matrix(open(SAVED_MERGED_PASSPORTS_FPATH))
    return genotypes, physical_map, genetic_map, passports


def _calculate_binomial_confidence_intervals(counts, total_counts):

    df11 = 2 * (total_counts - counts + 1)
    df12 = 2 * counts
    df21 = 2 * (counts + 1)
    df22 = 2 * (total_counts - counts)

    f_distrib_1 = distributions.f.ppf(0.975, df11, df12)
    f_distrib_2 = distributions.f.ppf(0.975, df21, df22)

    if not math.isnan(f_distrib_1):
        lower_ci = counts / ((total_counts - counts + 1) * f_distrib_1 + counts)
    else:
        lower_ci = 0

    if not math.isnan(f_distrib_2):
        denom = total_counts - counts + (counts + 1) * f_distrib_2
        upper_ci = (counts + 1) * f_distrib_2 / denom
    else:
        upper_ci = 1

    return [lower_ci, upper_ci]


def _calculate_allele_freqs(genotypes, classification, binomial_ci=False,
                            min_num_alleles=2):
    freqs = {}
    pop_size = {}

    if binomial_ci:
        conf_intervals = []
        multiindex = []
    else:
        conf_intervals = None

    accs_with_shape = genotypes.index
    accs_for_analysis = set(classification.index).intersection(accs_with_shape)
    accs_for_analysis = list(accs_for_analysis)
    classification = classification.ix[accs_for_analysis]

    populations = classification.groupby(classification.values).groups
    for pop, accs in populations.viewitems():
        alleles = []
        for acc in accs:
            acc_genotype = list(genotypes.ix[acc])
            alleles.extend(acc_genotype)
        al_counts = Counter(alleles)
        tot_alleles = sum(al_counts.values())
        if tot_alleles < min_num_alleles:
            continue
        else:
            al_freqs = {}
            for allele, count in al_counts.viewitems():
                al_freqs[allele] = count / tot_alleles
                if binomial_ci:
                    ci = _calculate_binomial_confidence_intervals(count,
                                                                  tot_alleles)
                    conf_intervals.append(ci)
                    multiindex.append((pop, allele))

            if binomial_ci and len(al_freqs) == 1:
                lower_ci, upper_ci = conf_intervals[-1]
                conf_intervals.append([1 - upper_ci, 1 - lower_ci])

                allele = multiindex[-1][1]
                allele = 'AB'.replace(allele, '')
                multiindex.append((pop, allele))

        freqs[pop] = al_freqs
        pop_size[pop] = len(alleles) / 2

    if binomial_ci:
        conf_intervals = DataFrame(conf_intervals)
        multiindex = MultiIndex.from_tuples(multiindex)
        conf_intervals.index = multiindex
        conf_intervals.columns = ['lower_ci', 'upper_ci']

    return DataFrame(freqs).T.fillna(0), Series(pop_size), conf_intervals


def plot_shape_histogram(dir_path, shape_genotypes, classification,
                         binomial_ci=False, min_num_alleles=12):

    for marker in shape_genotypes:
        classified_inds = classification.data.index
        genotyped_inds = shape_genotypes[marker].data.index
        acc_analysis = set(classified_inds).intersection(genotyped_inds)
        shape_genotype = shape_genotypes[marker].data.ix[acc_analysis]
        shape_genotype = shape_genotype.dropna()

        for group in ['species', 'group1', 'group2']:
            shape_dir_grp = pjoin(dir_path, group)
            if not os.path.exists(shape_dir_grp):
                os.mkdir(shape_dir_grp)
            plot_classification = classification.data[group]

            freqs, pop_size, ci = _calculate_allele_freqs(shape_genotype,
                                                          plot_classification,
                                                       binomial_ci=binomial_ci,
                                               min_num_alleles=min_num_alleles)

            freqs = freqs.reindex(GROUPING_ORDER[group]).dropna()

            fhand = open(pjoin(shape_dir_grp, 'histogram.' + marker + '.' +
                               group + '.svg'), 'w')
            args = [freqs, fhand]
            kwargs = {}
            if ci is not None:
                upper_ci = ci['upper_ci'].unstack(level=1)
                upper_ci = upper_ci.reindex(GROUPING_ORDER[group]).dropna()
                upper_ci['B'] = 0
                lower_ci = ci['lower_ci'].unstack(level=1)
                lower_ci = lower_ci.reindex(GROUPING_ORDER[group]).dropna()
                lower_ci['B'] = 0
                kwargs['y_error_bars_lower'] = freqs - lower_ci
                kwargs['y_error_bars_upper'] = upper_ci - freqs

            stacked_bars(*args, **kwargs)


            pop_size = pop_size.reindex(GROUPING_ORDER[group]).dropna()
            freqs['pop_size'] = pop_size
            fhand = open(pjoin(shape_dir_grp, 'pop_size.' + marker + '.' +
                               group + '.csv'), 'w')
            freqs.to_csv(fhand)


def read_coancestry_matrix():

    chromopaint_dir = pjoin(RESULTS_DIR, 'chromopainter')
    ind_order_fhand = open(pjoin(chromopaint_dir, "Individual_names.csv"))

    coancestry_matrix_fhand = open(pjoin(chromopaint_dir,
                                        'combined_output.chunkcounts.out'))
    coancestry_matrix = load_dataframe_csv(coancestry_matrix_fhand,
                                           index_col=0, header=1, sep=' ',
                                           meta={INDIVIDUALS_IN_ROWS: True})

    ind_codes = load_dataframe_csv(ind_order_fhand, index_col=0, header=0)
    coancestry_matrix.data.index = ind_codes.data['Ind_name']
    coancestry_matrix.data.columns = ind_codes.data['Ind_name']

    return coancestry_matrix


def plot_coancestry(matrix, classification):
    matrix = matrix.data

    ordering_criterum = 'group2'
    grps = classification.data.groupby(ordering_criterum).groups
    ordered_accs = []
    for grp in GROUPING_ORDER[ordering_criterum]:
        ordered_accs.extend(grps[grp])
    matrix = matrix.reindex(ordered_accs).dropna()
    ordered_accs = [acc for acc in ordered_accs if acc in matrix.columns]
    matrix = matrix.reindex_axis(ordered_accs, axis=1)

    fhand = open(pjoin(RESULTS_DIR, 'coancestry_matrix_indi.png'), 'w')
    plot_matrix(matrix, fhand)


def plot_grouped_coancestry(matrix, classification, groupby):

    ids = matrix.data.index

    classification = classification.data.ix[ids]
    groups = classification.groupby(groupby).groups

    means = {}

    for row_pop, row_accs in groups.viewitems():
        nrow_accs = len(row_accs)
        row_sum = matrix.data.ix[row_accs].sum()
        col_pop_means = {}
        for col_pop, col_accs in groups.viewitems():
            ncol_accs = len(col_accs)
            if row_accs == col_accs:
                ntot = nrow_accs * ncol_accs - nrow_accs
            else:
                ntot = nrow_accs * ncol_accs

            pop_mean = row_sum[col_accs].sum() / ntot
            col_pop_means[col_pop] = pop_mean
        mean = Series(col_pop_means)
        means[row_pop] = mean

    means = DataFrame(means).T.fillna(0)

    means = means.reindex(GROUPING_ORDER[groupby]).dropna()

    groups = [grp for grp in GROUPING_ORDER[groupby] if grp in means.columns]
    means = means.reindex_axis(groups, axis=1)

    fhand = open(pjoin(RESULTS_DIR, 'coancestry_means_' + groupby + '.svg'),
                 'w')
    stacked_bars(means, fhand)

    fhand = open(pjoin(RESULTS_DIR, 'coancestry_matrix_' + groupby + '.svg'),
                 'w')
    plot_matrix(means, fhand)


GROUP_COLORS = {}


def _get_color_for_group2(group2):

    if group2 in GROUP_COLORS:
        return GROUP_COLORS[group2]

    hex_ = random.choice(['4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd',
                          'e', 'f'])
    color = ''.join([random.choice(hex_) for i in range(6)])
    color = 'ff' + color

    GROUP_COLORS[group2] = color
    return color


def write_kml_file(acc_info, kml_fpath):
    kml = simplekml.Kml()
    kml.document.name = 'Tomato joint variation analysis'

    species_styles = {
           'SLC': {'icon': 'http://maps.google.com/mapfiles/kml/paddle/C.png'},
           'SLL': {'icon': 'http://maps.google.com/mapfiles/kml/paddle/L.png'},
           'SP': {'icon': 'http://maps.google.com/mapfiles/kml/paddle/P.png'},
       'mixture': {'icon': 'http://maps.google.com/mapfiles/kml/paddle/M.png'},
          'wild': {'icon': 'http://maps.google.com/mapfiles/kml/paddle/W.png'},
          'SC': {'icon': 'http://maps.google.com/mapfiles/kml/paddle/D.png'},
          'SG': {'icon': 'http://maps.google.com/mapfiles/kml/paddle/G.png'},
          'Schm': {'icon': 'http://maps.google.com/mapfiles/kml/paddle/H.png'},
          'SN': {'icon': 'http://maps.google.com/mapfiles/kml/paddle/N.png'},
                    }

    folders = {}
    for acc_name, acc_pass in acc_info.data.iterrows():
        species = acc_pass['species']
        group1 = acc_pass['group1']
        group2 = acc_pass['group2']
        passport = acc_pass['Passport_Classification']

        if acc_pass['Longitude'] == 'None':
            continue
        provenance = acc_pass['provenance']
        lon = lon_to_deg(acc_pass['Longitude'])
        lat = lat_to_deg(acc_pass['Latitude'])

        if species not in folders:
            folder = kml.newfolder(name=species)
            folders[species] = {'folder': folder, 'subfolders': {}}
        if group1 not in folders[species]['subfolders']:
            parent_folder = folders[species]['folder']
            folder = parent_folder.newfolder(name=group1)
            folders[species]['subfolders'][group1] = {'folder': folder,
                                                      'subfolders': {}}
        if group2 not in folders[species]['subfolders'][group1]['subfolders']:
            parent_folder = folders[species]['subfolders'][group1]['folder']
            folder = parent_folder.newfolder(name=group2)
            style = simplekml.StyleMap()
            icon_url = species_styles[species]['icon']
            color = color = 'ff88ff00'
            colormode = 'normal'
            style.normalstyle.iconstyle.icon.href = icon_url
            style.highlightstyle.iconstyle.icon.href = icon_url
            style.highlightstyle.iconstyle.scale = 1.1
            style.normalstyle.labelstyle.scale = 0
            style.normalstyle.iconstyle.color = color
            style.normalstyle.iconstyle.colormode = colormode
            style.highlightstyle.iconstyle.color = color
            style.highlightstyle.iconstyle.colormode = colormode

            folders[species]['subfolders'][group1]['subfolders'][group2] = {'folder': folder,
                                                                            'style': style}

        folder = folders[species]['subfolders'][group1]['subfolders'][group2]['folder']
        style = folders[species]['subfolders'][group1]['subfolders'][group2]['style']

        description = '<p>Species: <strong>' + species + '</strong></p>'
        description += '<p>Group: <strong>' + group1 + '</strong></p>'
        description += '<p>Subgroup: <strong>' + group2 + '</strong></p>'
        description += '<p>Passport group: <strong>' + passport + '</strong></p>'

        pnt = folder.newpoint(name=acc_name, coords=[(lon, lat)],
                              description=description)
        pnt.stylemap = style

    kml.save(kml_fpath)


def calculate_dist_distributions(genotypes, acc_info, dists_fhand):
    genotype_accs = genotypes.data.index
    acc_info_accs = acc_info.data.index
    accs_for_analysis = set(acc_info_accs).intersection(genotype_accs)
    genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                               genotypes.meta)
    acc_info = MatrixWithMeta(acc_info.data.ix[accs_for_analysis],
                              acc_info.meta)

    dists = pairwise_individual_genetic_distance(genotypes)

    dist_per_group = {}
    dist_per_group['all'] = squareform(dists, force='tovector')

    for grouping in ['species', 'group1', 'group2', 'accession']:
        group_dists = []
        for group, accs in acc_info.data.groupby(grouping).groups.items():
            in_group_dists = dists.ix[accs][accs]
            in_group_dists = squareform(in_group_dists, force='tovector')
            group_dists.extend(in_group_dists)
        dist_per_group[grouping] = group_dists

    # All should have the same name of items to be printed
    max_n_dists = max(len(dists) for dists in dist_per_group.values())
    dists_per_group = {}
    for group, dists in dist_per_group.items():
        to_fill = max_n_dists - len(dists)
        dists = list(dists)
        dists.extend([float('nan')] * to_fill)
        dists_per_group[group] = dists
    dists = DataFrame(dists_per_group)
    dists = dists[['all', 'species', 'group1', 'group2', 'accession']]
    dists.to_csv(dists_fhand)


def load_color_marker_definitions():
    figure_coding_fhand = open(pjoin(ANALYSIS_DIR, 'color_map.csv'))
    figure_coding_fhand.readline()
    color_definitions = defaultdict(dict)
    marker_definitions = defaultdict(dict)
    marker_definitions2 = defaultdict(dict)
    for line in figure_coding_fhand.readlines():
        line = line.strip().split('\t')
        color_definitions[line[0]].update({line[1]: line[2]})
        marker_definitions[line[0]].update({line[1]: line[3]})
        marker_definitions2[line[0]].update({line[1]: line[4]})

    color_definitions = {group: ColorMap(colors) for group, colors in color_definitions.items()}
    marker_definitions = {group: MarkerMap(colors) for group, colors in marker_definitions.items()}
    marker_definitions2 = {group: MarkerMap(colors) for group, colors in marker_definitions2.items()}
    marker_definitions2[line[0]].update({line[1]: line[4]})

    return color_definitions, marker_definitions, marker_definitions2


def filter_out_of_geographic_region(df, limits):

    keep_accs = []
    for acc in df.index:
        if ((limits[0] < df.ix[acc]['Longitude'] < limits[1]) and
            (limits[2] < df.ix[acc]['Latitude'] < limits[3])):
            keep_accs.append(acc)
    return df.ix[keep_accs]


def write_distances_output(fhand, genetic_dist_matrix, geographic_dist_matrix):
    line = 'ind1\tind2\tgenetic\tdistance\n'
    fhand.write(line)
    n_indiv = genetic_dist_matrix.shape[0]
    for i in xrange(n_indiv):
        for j in xrange(i + 1, n_indiv):
            line = '%s\t%s\t%s\t' % (i, j, genetic_dist_matrix.iloc[i, j])
            dist = geographic_dist_matrix.iloc[i, j]
            if dist == 0:
                dist = 0.000001
            line += '%s\n' % (dist)
            fhand.write(line)
    fhand.flush()


def write_arlequin(fpath, genotypes, inner_classification,
                   outter_classification=None):

    pop_class = inner_classification.groupby(inner_classification).groups
    if outter_classification is not None:
        struct_class = outter_classification.groupby(outter_classification).groups

    arp_file = open(fpath, 'w')

    arp_header = '''[Profile]
Title="Input Arlequin"
GenotypicData=1
GameticPhase=0
RecessiveData=0
DataType=STANDARD
LocusSeparator=WHITESPACE
MissingData='0'
'''
    arp_header += 'NbSamples=%s\n' % len(pop_class)
    arp_file.write(arp_header)

    sample_block = '''[Data]
[[Samples]]
'''
    codec = get_codec_from_genotypes(genotypes)

    for pop, indvs in pop_class.viewitems():
        sample_block += 'SampleName="%s"\n' % pop
        sample_block += 'SampleSize=%s\n' % len(indvs)
        sample_block += 'SampleData=    {\n'
        for indv in indvs:
            sample_block += '%s\t1\t' % indv
            indv_geno = genotypes.data.ix[indv]
            hap1, hap2 = '', ''
            for locus in indv_geno:
                decoded = codec.decode_to_genotype(locus)
                if decoded is None:
                    decoded = ('0', '0')
                hap1 += decoded[0] + ' '
                hap2 += decoded[1] + ' '
            sample_block += hap1.strip() + '\n'
            sample_block += '\t\t' + hap2.strip() + '\n'
        sample_block += '}\n'
    arp_file.write(sample_block)

    if outter_classification is not None:
        structure_block = '''[[Structure]]

\tStructureName="Structure"
'''
        structure_block += '\tNbGroups=%s\n\n' % len(struct_class)

        for group, indvs in struct_class.viewitems():
            pops = list(set(inner_classification.ix[indvs]))
            structure_block += '\tGroup={\n'
            for pop in pops:
                structure_block += '\t\t"%s"\n' % pop
            structure_block += '\t}\n\n'
        arp_file.write(structure_block)
    arp_file.close()


def _nexus_two_state_format(indv_genotype, codec, coding):
    formated_geno = []
    #It keeps the first allele that founds as a reference allele
    if coding is None:
        coding = [snp[0] if snp is not None else None for snp in indv_genotype]
    for index in range(len(indv_genotype)):
        ref = coding[index]
        snp = indv_genotype[index]
        if (ref is None) and (snp is not None):
            coding[index] = snp[0]
            ref = coding[index]
        if snp is None:
            formated_geno.append('?')
        else:
            if (snp[0] == snp[1]) and (snp[0] == ref):
                formated_geno.append('0')
            elif (snp[0] != snp[1]) and (snp[0] == ref or snp[1] == ref):
                formated_geno.append('1')
            elif (snp[0] == snp[1]) and (snp[0] != ref):
                formated_geno.append('2')
            else:
                msg = 'Genotype %s cannot be converted ' % (snp, )
                raise ValueError(msg)
    return formated_geno, coding


def write_nexus_file(fpath, genotypes, two_state_format=True):

    fhand = open(fpath, 'w')

    fhand.write('#NEXUS\n')
    data_block = 'Begin data;\n'
    data_block += '\tdimensions ntax=%s nchar=%s;\n' % (len(genotypes.data.index),
                                                   len(genotypes.data.columns))
    data_block += '\tformat datatype=standard symbols="012" missing=?;\nMatrix\n'
    fhand.write(data_block)

    codec = get_codec_from_genotypes(genotypes)
    coding = None
    for acc in genotypes.data.index:
        fhand.write(acc + '\t')
        indv_genotype = list(genotypes.data.ix[acc].apply(codec.decode_to_ints))
        if two_state_format:
            formated_geno, coding = _nexus_two_state_format(indv_genotype,
                                                            codec, coding)
        fhand.write(('').join(formated_geno) + '\n')

    data_block = ';\nEND;\n'
    fhand.write(data_block)
    fhand.close()


def generate_acc_syns(acc_info):

    acc_syns = {}
    for sample in acc_info.data.index:
        acc_name = acc_info.data.ix[sample]['accession']
        acc_syns[sample] = sample
        acc_syns[acc_name] = sample
        if ' ' in sample:
            acc_syns[sample] = sample
            acc_syns[sample.replace(' ', '_')] = sample
            acc_syns[sample.replace(' ', '')] = sample
            acc_syns[sample.replace(' ', '-')] = sample
        if '_' in sample:
            acc_syns[sample] = sample
            acc_syns[sample.replace('_', ' ')] = sample
            acc_syns[sample.replace('_', '')] = sample
            acc_syns[sample.replace('_', '-')] = sample
        if '-' in sample:
            acc_syns[sample] = sample
            acc_syns[sample.replace('-', ' ')] = sample
            acc_syns[sample.replace('-', '')] = sample
            acc_syns[sample.replace('-', '_')] = sample
        if  sample == 'Pomodoro':
            acc_syns[sample + ('_')] = sample
        if  sample == 'Brandywine':
            acc_syns[sample + ('_')] = sample
        if  'Chonto' in sample:
            acc_syns['1091-Chonto_21_'] = sample

        if ' ' in acc_name:
            acc_syns[acc_name] = sample
            acc_syns[acc_name.replace(' ', '_')] = sample
            acc_syns[acc_name.replace(' ', '')] = sample
            acc_syns[acc_name.replace(' ', '-')] = sample
        if '_' in acc_name:
            acc_syns[acc_name] = sample
            acc_syns[acc_name.replace('_', ' ')] = sample
            acc_syns[acc_name.replace('_', '')] = sample
            acc_syns[acc_name.replace('_', '-')] = sample
        if '-' in acc_name:
            acc_syns[acc_name] = sample
            acc_syns[acc_name.replace('-', ' ')] = sample
            acc_syns[acc_name.replace('-', '')] = sample
            acc_syns[acc_name.replace('-', '_')] = sample

    return acc_syns


def write_gephi(fpath, network):
    fhand = open(fpath, 'w')
    line = '''<gexf version="1.1">
<graph defaultedgetype="undirected" idtype="string" type="static">
'''
    fhand.write(line)
    fhand.write('<nodes count="%s">\n' % len(network.index))
    chrom_names = {}
    for indx, chrom_name in enumerate(network.index):
        chrom_names[chrom_name] = indx
        fhand.write('<node id="%.1f" label="%s_%s"/>\n' % (indx, chrom_name[0], chrom_name[1]))
    fhand.write('</nodes>\n')
    n_edges = 0
    lines = ''
    for chrom_name, neighbours in network['neighbours'].iteritems():
        for neighbour in neighbours:
            lines += '<edge id="%s" source="%.1f" target="%.1f"/>\n' % (n_edges, chrom_names[neighbour], chrom_names[chrom_name])
            n_edges += 1
    fhand.write('<edges count="%s">\n' % n_edges)
    fhand.write(lines)
    fhand.write('</edges>\n</graph>\n</gexf>')
    fhand.close()


def write_merged_passports(acc_info):
    out_dir = pjoin(RESULTS_DIR, 'merged_data')
    out_fpath = pjoin(out_dir, 'merged_passports.csv')
    acc_info.data.to_csv(out_fpath, sep='\t', index_label='sample')


def main():

    if not exists(PROCESSED_DATA_DIR):
        os.mkdir(PROCESSED_DATA_DIR)

    filter_missing_data = True
    filter_monomorphic = False
    filter_close_data = True
    filter_missing_koenig_genotypes = False
    one_sample_per_accesion = True

    color_definitions, marker_definitions, marker_definitions2 = load_color_marker_definitions()

    if len(listdir(PROCESSED_DATA_DIR)) < 4:
        merge_filter_data(filter_missing_data=filter_missing_data,
                          filter_monomorphic=filter_monomorphic,
                          filter_close_data=filter_close_data,
                         filter_missing_koenig=filter_missing_koenig_genotypes)

    genotypes, physical_map, genetic_map, acc_info = load_data_matrices()
    #open('results/markers.txt', 'w').write('\n'.join(genotypes.data.columns))

    basic_stats = False
    global_genotypes_stats = False
    check_heterozigotes = False
    pop_stats = False
    pca_analysis = False
    phase = False
    finestructure = False
    adegenet = False
    plot_geomap = True
    shape = False
    shape_table_by_sample = False
    shape_haplotypes = False
    plot_coancestry_matrix = False
    population_distance = False
    rarefaction = False
    allele_freq_distribution = False
    kml = False
    stats_per_grid = False
    plot_stats_per_grid = False
    adze = False
    pairwise_adze = False
    do_structure = False
    structure_results = False
    dist_distributions = False
    indiv_genetic_tree = False
    isolation_by_distance = False
    climatic_correlation = False
    climatic_geographic_correlation = False
    arlequin = False
    snapp = False
    freq_by_locus = False
    genalex = False
    genepop_file = False
    bayescan_file = False
    tassel_file = False
    tassel_file_from_phased = False
    LD = False
    write_passports = False

    # filter samples
    filter_condition = {'group1': ['failed_geno']}
    genotypes = filter_samples(genotypes, acc_info, filter_condition,
                               reverse=True)

    # classification:
    if True:

        del acc_info.data['species']
        del acc_info.data['group1']

        classification_fpath = pjoin(RESULTS_DIR, 'classification',
                                     'classification.csv')
        classification = read_classification(open(classification_fpath))

        acc_info_data = dframe_append_columns(acc_info.data,
                                       [classification.data['accession']])
        acc_info_data = dframe_append_columns(acc_info_data,
                                       [classification.data['species']])
        acc_info_data = dframe_append_columns(acc_info_data,
                                      [classification.data['group1']])
        acc_info_data = dframe_append_columns(acc_info_data,
                                      [classification.data['group2']])
        acc_info_data = dframe_append_columns(acc_info_data,
                                      [classification.data['passport_group']])
        acc_info_data = dframe_append_columns(acc_info_data,
                                 [classification.data['heterozigosity_index']])
        acc_info_data = dframe_append_columns(acc_info_data,
                                 [classification.data['Geographic_Region']])
        acc_info_data = dframe_append_columns(acc_info_data,
                              [classification.data['Passport_Classification']])
        acc_info_data = dframe_append_columns(acc_info_data,
                                        [classification.data['Rarefaction_1']])
        acc_info_data = dframe_append_columns(acc_info_data,
                                     [classification.data['Passport_Species']])
        acc_info_data = dframe_append_columns(acc_info_data,
                                     [classification.data['original_bank']])
        acc_info_data = dframe_append_columns(acc_info_data,
                                     [classification.data['source_bank']])
        acc_info_data = dframe_append_columns(acc_info_data,
                                     [classification.data['%Het']])
        acc_info_data = dframe_append_columns(acc_info_data,
                                     [classification.data['collection_source']])
        acc_info_data = dframe_append_columns(acc_info_data,
                                     [classification.data['Collection_year']])
        acc_info_data = dframe_append_columns(acc_info_data,
                                     [classification.data['Country']])


        # Add selection for Esther 2014 project
        project_esther_dir = pjoin(ORIG_DATA_DIR, 'project_esther_2014')
        selection_fpath = pjoin(project_esther_dir,
                                'datos_integrados_por_mj_20140604.csv')
        # there are rows with the same acc name
        uniq_fhand = StringIO()
        accs = set()
        for line in open(selection_fpath):
            items = line.split('\t')

            acc = items[0] if not(items[1]) else items[1]
            if not acc:
                continue
            if acc in accs:
                continue
            accs.add(acc)
            items[0] = acc
            uniq_fhand.write('\t'.join(items))
        uniq_fhand.seek(0)
        classification = read_classification(uniq_fhand)
        classification = classification.data.reindex(acc_info_data.index)

        #print set(classification['Seleccion_20140325'])
        selection_dict = {float('nan'): 'no_select', 'X': 'selected',
                          'x': 'selected', '1': 'selected', 1: 'selected'}
        selection = classification['Selecc. 2014-06-04'].map(selection_dict)
        selection.name = 'selection_esther'
        acc_info_data = dframe_append_columns(acc_info_data, [selection])
        #print classification.columns
        s2_val = classification['S2 seedValencia']
        s2_esther = classification['sent Esther']
        s2_tunel = classification['Reprod tunel parcela 2013']
        s2s = {}
        for acc in s2_val.index:
            if str(s2_val[acc]) != 'nan' and s2_val[acc]:
                s2 = 's2_val' + str(int(s2_val[acc]))
            elif str(s2_esther[acc]) != 'nan' and s2_esther[acc]:
                s2 = 's2_Esther'
            elif str(s2_tunel[acc]) != 'nan' and s2_tunel[acc]:
                s2 = 's2_tunel'
            else:
                s2 = 'no_s2'
            s2s[acc] = s2

        s2 = Series(s2s, name='s2')

        acc_info_data = dframe_append_columns(acc_info_data, [s2])

        acc_info = MatrixWithMeta(acc_info_data, acc_info.meta)

    print 'Excluded samples filtered'
    species = ['SG', 'SP', 'SLL', 'SLC', 'mixture', 'SN', 'Schm', 'SC']
    #species = ['SP']
    filter_condition = {'species': species}
    genotypes = filter_samples(genotypes, acc_info, filter_condition,
                               reverse=False)

    classification2 = ['SLL_fresh', 'SLL_processing_1', 'SLL_processing_2',
                       'SLL1', 'SLC_Wva106', 'SLC_LA2135',
                       'SLC_LA2792', 'SLC_1', 'SLL_1']
    filter_condition = {'group1': classification2}
    #genotypes = filter_samples(genotypes, acc_info, filter_condition,
    #                           reverse=True)

    classification3 = ['SLC_Asia', 'SLC_Colombia', 'SLC_Costa_Rica',
                       'SLC_world']
    filter_condition = {'group2': classification3}
    #genotypes = filter_samples(genotypes, acc_info, filter_condition,
    #                           reverse=True)

    classificationpass = []
    filter_condition = {'Rarefaction_1': classificationpass}
#     genotypes = filter_samples(genotypes, acc_info, filter_condition,
#                                 reverse=False)

    samples = []
    samples = list(set(genotypes.data.index).intersection(samples))
    if samples:
        genotypes = MatrixWithMeta(genotypes.data.drop(samples),
                                   metadata=genotypes.meta)

    genotypes_accs = genotypes.data.index
    acc_info_accs = acc_info.data.index
    accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)

    print accs_for_analysis

    genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                               genotypes.meta)
    acc_info = MatrixWithMeta(acc_info.data.ix[accs_for_analysis],
                              acc_info.meta)

    if True:
        print 'Accessions with samples that do not match'
        groups_by_acc = acc_info.data.groupby('accession').groups
        repeated_accs = {acc: samples for acc, samples in groups_by_acc.viewitems() if len(samples) > 1}
        for acc, samples in repeated_accs.items():
            acc_classes = [(s, acc_info.data.get_value(s, 'group2')) for s in samples]
            if len(set([class_[1] for class_ in acc_classes])) > 1:
                print acc, acc_classes

    if write_passports:
        write_merged_passports(acc_info)

    if one_sample_per_accesion:
        # We have to keep one sample per accession
        repeated_accs = acc_info.data.groupby('accession').groups
        samples_to_remove = []
        for acc, samples in repeated_accs.items():
            if len(samples) == 1:
                continue
            random.shuffle(samples)
            if 'koenig' in samples[0]:
                samples = samples[1:] + [samples[0]]
            samples_to_remove.extend(samples[1:])
        acc_info = MatrixWithMeta(acc_info.data.drop(samples_to_remove),
                                  acc_info.meta)

    if basic_stats:
        print 'Provenance'
        print pprint(dict(Counter(acc_info.data.provenance)))
        print 'total: ', sum(Counter(acc_info.data.provenance).values())
        print 'Passport species'
        counts = dict(Counter(acc_info.data['Passport_Species']))
        print pprint(counts)
        print 'total: ', sum(counts.values())
        print 'Genetic species'
        counts = dict(Counter(acc_info.data['species']))
        print pprint(counts)
        print 'total: ', sum(counts.values())

    if global_genotypes_stats:

        maf_threshold = 0.95
        missing_threshold = 0.1

        n_missing_markers = 0
        n_variable_markers = 0

        for marker, genos in genotypes.data.iteritems():
            counts = Counter(genos)
            total_counts = sum(counts.values())
            if 0 in counts:
                if counts[0] / total_counts > missing_threshold:
                    n_missing_markers += 1
                n_missing = counts[0]
                del counts[0]
                total_counts -= n_missing
            n_freq_geno = max(counts.values())
            if n_freq_geno / total_counts < maf_threshold:
                n_variable_markers += 1

        line = 'Numbers of markers with more that'
        line += ' a %s%% of missing_markers: %s' % (missing_threshold * 100,
                                                n_missing_markers)
        print line
        line = 'Numbers of markers with a MAF lower than'
        line += ' %s%%: %s' % (maf_threshold * 100, n_variable_markers)
        print line

    if not filter_missing_koenig_genotypes:
        koenig_samples = ['LA1589_koenig', 'LA0530', 'SLL_koenig']
        koenig_samples = set(koenig_samples).intersection(genotypes.data.index)
        if koenig_samples:
            genotypes = MatrixWithMeta(genotypes.data.drop(koenig_samples),
                                       metadata=genotypes.meta)

    if pop_stats:
        stats_dir = pjoin(RESULTS_DIR, 'pop_stats')
        if not exists(stats_dir):
            os.mkdir(stats_dir)

        for group in ['species', 'group1', 'group2']:
            stats_fhand = open(pjoin(stats_dir, 'pop_stats.' + group +
                                     '.csv'), 'w')
            classification = acc_info.data[group].ix[genotypes.data.index]
            tgenotypes = transpose_genotypes(genotypes)
            stats = calculate_basic_stats(tgenotypes, classification,
                                          min_num_individuals=5)
            stats = stats.dropna(axis=1)
            sorted_group_labels = [group for group in GROUPING_ORDER[group] if group in stats.columns]
            stats[sorted_group_labels].T.to_csv(stats_fhand)
            stats_fhand.close()

    if adegenet:
        adegenet_dir_path = pjoin(RESULTS_DIR, 'adegenet')
        if not exists(adegenet_dir_path):
            os.mkdir(adegenet_dir_path)

        geno_fhand = open(pjoin(adegenet_dir_path, 'adegenet_geno.snp'), 'w')
        geo_fhand = open(pjoin(adegenet_dir_path, 'adegenet_geo.csv'), 'w')
        pop_group = 'group2'

        acc_info_geo = acc_info.data.ix[genotypes.data.index]
        acc_info_geo = acc_info_geo[['Longitude', 'Latitude']]
        acc_info_geo = acc_info_geo.replace('None', float('nan'))
        acc_info_geo = acc_info_geo.dropna(how='any')
        acc_info_geo['Longitude'] = acc_info_geo['Longitude'].apply(lon_to_deg)
        acc_info_geo['Latitude'] = acc_info_geo['Latitude'].apply(lat_to_deg)
        acc_info_geo = acc_info_geo.dropna(how='any')

        indv_with_geo_info = acc_info_geo.index

        if genotypes.meta[INDIVIDUALS_IN_ROWS]:
            genotypes = MatrixWithMeta(genotypes.data.ix[indv_with_geo_info],
                                       genotypes.meta)
        else:
            genotypes = MatrixWithMeta(genotypes.data[indv_with_geo_info],
                                       genotypes.meta)
        indv_classification = acc_info_data.ix[indv_with_geo_info][pop_group]

        create_file_snp_fomat(geno_fhand, genotypes, indv_classification)
        acc_info_geo.to_csv(geo_fhand)
        locus_info = open(pjoin(adegenet_dir_path, 'adegenet_locus.csv'), 'w')

        if genotypes.meta[INDIVIDUALS_IN_ROWS]:
            strng = '\n'.join(genotypes.data.columns)
        else:
            strng = '\n'.join(genotypes.data.index)
        locus_info.write(strng)

        geno_fhand.close()

    if phase:
        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)
        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)

        classification = acc_info.data['species']
        phase_dir_path = pjoin(RESULTS_DIR, 'PHASE')
        run_fast_phase(genotypes, physical_map, dir_path=phase_dir_path,
                       individual_classification=classification,
                       parameters=['-H200', '-i'])

    if finestructure:
        interpoled_map_fhand = open(INTERPOLED_RECOM_MAP)
        interpolated_gen_map = load_markers_map(interpoled_map_fhand,
                                                molecule_col='chromosome',
                                               location_col='genetic_position',
                                                index_col=0, sep='\t')
        interpoled_map_fhand = open(INTERPOLED_RECOM_MAP)
        interpolated_physical_map = load_markers_map(interpoled_map_fhand,
                                                   molecule_col='chromosome',
                                                   location_col='position',
                                                   index_col=0, sep='\t')
        genetic_position = interpolated_gen_map.data['genetic_position']
        gen_dist = genetic_position.map(lambda x: x / 100)
        interpolated_gen_map.data['genetic_position'] = gen_dist

        phase_dir_path = pjoin(RESULTS_DIR, 'PHASE')
        chromopainter_dir_path = pjoin(RESULTS_DIR, 'chromopainter')
        if not exists(chromopainter_dir_path):
            os.mkdir(chromopainter_dir_path)

        individual_info_fhand = open(pjoin(chromopainter_dir_path,
                                           'individual_names.txt'), 'w')
        individual_info = False
        haplotypes = parse_fast_phase_results(phase_dir_path)

        hap_type = 'individual_haplotypes'
        #hap_type = 'switch_haplotypes'

        for chrom in haplotypes[hap_type].keys():
            chrom_dir_path = pjoin(chromopainter_dir_path, chrom)
            if not exists(chrom_dir_path):
                os.mkdir(chrom_dir_path)

            chrom_haplotypes = haplotypes[hap_type][chrom]
            indv = acc_info.data.index
            chrom_haplotypes = MatrixWithMeta(chrom_haplotypes.data.ix[indv],
                                              metadata=chrom_haplotypes.meta)

            if individual_info is False:
                indiv = chrom_haplotypes.data.index.get_level_values(0)
                indiv = [indiv[i] for i in range(len(indiv)) if i % 2 != 0]
                individual_info_fhand.write('\n'.join(indiv))
                individual_info_fhand.close()
                individual_info = True

            create_chromopainter_infiles(chrom_haplotypes,
                                         interpolated_physical_map,
                                         interpolated_gen_map,
                                         chrom_dir_path)

    if check_heterozigotes:
        stats_dir = pjoin(RESULTS_DIR, 'pop_stats')
        if not exists(stats_dir):
            os.mkdir(stats_dir)

        hets = calculate_individual_heterozigosity(genotypes)
        hets.to_csv(pjoin(stats_dir, 'individual_heterozygosity.csv'))


    labels = True
#    marker_col = 'provenance'
    #marker_col = 'passport_group'
    marker_col = 'group1'
    #marker_col = 'species'
    #marker_col = 's2'
#    color_col = 'species'
    color_col = 'group2'
#    color_col = 'Geographic_Region'
    #color_col = 'selection_esther'

    if pca_analysis:
        do_multivariate(genotypes, acc_info=acc_info,
                        marker_col=marker_col, color_col=color_col,
                        labels=labels, color_definitions=color_definitions,
                        marker_definitions=marker_definitions)

    if plot_geomap:

        geo_dir = pjoin(RESULTS_DIR, 'geomap')

        if not exists(geo_dir):
            os.mkdir(geo_dir)

        individuals = genotypes.data.index
        cols = ['accession', 'species', 'Longitude', 'Latitude']
        acc_info_geo = acc_info.data.ix[individuals][cols].replace('None', '')
        acc_info_geo['Longitude'] = acc_info_geo['Longitude'].map(lon_to_deg)
        acc_info_geo['Latitude'] = acc_info_geo['Latitude'].map(lat_to_deg)
        acc_info_geo[color_col] = acc_info.data[color_col]
        acc_info_geo[marker_col] = acc_info.data[marker_col]

        all_indiv_geo = acc_info_geo.index
        acc_info_geo = acc_info_geo.dropna()
        indiv_with_geo = acc_info_geo.index

        print 'Individuals lacking geographical information'
        print set(all_indiv_geo).difference(indiv_with_geo)

        background_image = WORLD_CLIMATE_IMAGE
        background_image = None

        #Andean
        plot_fhand = open(pjoin(geo_dir, 'map_andean.svg'), 'w')
        limits = (-82.0, -65, -18, 1.5)

        filtered_acc_info_geo = filter_out_of_geographic_region(acc_info_geo,
                                                                limits)

        plot_groups = get_scatter_groups_from_matrix(filtered_acc_info_geo,
                                                     x_col='Longitude',
                                                     y_col='Latitude',
                                                   marker_col=marker_col,
                                  marker_mapper=marker_definitions[marker_col],
                                                     color_col=color_col,
                                     color_mapper=color_definitions[color_col])

        geographic_map(plot_groups, plot_fhand, labels=labels,
                       background_image=background_image, limits=limits,
                       marker_size=65)

        #Ecuador
        plot_fhand = open(pjoin(geo_dir, 'map_ecuador.svg'), 'w')
        limits = (-81.2, -76, -7, 1.5)

        filtered_acc_info_geo = filter_out_of_geographic_region(acc_info_geo,
                                                                limits)

        plot_groups = get_scatter_groups_from_matrix(filtered_acc_info_geo,
                                                     x_col='Longitude',
                                                     y_col='Latitude',
                                                   marker_col=marker_col,
                                  marker_mapper=marker_definitions[marker_col],
                                                     color_col=color_col,
                                     color_mapper=color_definitions[color_col])

        geographic_map(plot_groups, plot_fhand, labels=labels,
                       background_image=background_image, limits=limits,
                       marker_size=85)

        #Peru
        plot_fhand = open(pjoin(geo_dir, 'map_peru.svg'), 'w')
        limits = (-82.0, -65, -18, -4)

        filtered_acc_info_geo = filter_out_of_geographic_region(acc_info_geo,
                                                                limits)

        plot_groups = get_scatter_groups_from_matrix(filtered_acc_info_geo,
                                                     x_col='Longitude',
                                                     y_col='Latitude',
                                                   marker_col=marker_col,
                                  marker_mapper=marker_definitions[marker_col],
                                                     color_col=color_col,
                                     color_mapper=color_definitions[color_col])

        geographic_map(plot_groups, plot_fhand, labels=labels,
                       background_image=background_image, limits=limits,
                       marker_size=85)

        #colombia and mesoamerica
        plot_fhand = open(pjoin(geo_dir, 'map_mesoamerica.svg'), 'w')
        #30.524413,-112.155763 1.318243,-74.714356
        limits = (-112, -73, 1.3, 30.5)

        filtered_acc_info_geo = filter_out_of_geographic_region(acc_info_geo,
                                                                limits)

        plot_groups = get_scatter_groups_from_matrix(filtered_acc_info_geo,
                                                     x_col='Longitude',
                                                     y_col='Latitude',
                                                   marker_col=marker_col,
                                  marker_mapper=marker_definitions[marker_col],
                                                     color_col=color_col,
                                     color_mapper=color_definitions[color_col])

        geographic_map(plot_groups, plot_fhand, labels=labels,
                       background_image=background_image, limits=limits,
                       marker_size=85)

        # all
        plot_fhand = open(pjoin(geo_dir, 'map_all.svg'), 'w')

        plot_groups = get_scatter_groups_from_matrix(acc_info_geo,
                                                     x_col='Longitude',
                                                     y_col='Latitude',
                                                   marker_col=marker_col,
                                marker_mapper=marker_definitions[marker_col],
                                                     color_col=color_col,
                                     color_mapper=color_definitions[color_col])

        geographic_map(plot_groups, plot_fhand, labels=labels,
                       background_image=background_image, marker_size=40)

    if shape:

        shape_dir = pjoin(RESULTS_DIR, 'shape', )

        if not exists(shape_dir):
            os.mkdir(shape_dir)

        acc_syns = generate_acc_syns(acc_info)
        shape_genotypes = merge_shape_genotypes(acc_syns)

        genes = ['lc', 'fw22', 'fw32', 'fw113', 'fas', 'ovate', 'sun', 'sun_dup']
        shape_geno_combined = DataFrame(None, index=acc_info.data.index,
                                        columns=genes)
        for gene in genes:
            shape_geno = shape_genotypes[gene].data
            acc_info_accs = acc_info.data.index
            common_accs = list(set(acc_info_accs).intersection(shape_geno.index))
            shape_geno = shape_geno.ix[common_accs]
            for acc, values in shape_geno.iterrows():
                if values[0] == values[1]:
                    shape_geno_combined.ix[acc][gene] = values[0]
                else:
                    shape_geno_combined.ix[acc][gene] = values[0] + '/' + values[1]

        fine_sample_names = []
        for sample in shape_geno_combined.index:
            try:
                fine_name = get_comav_acc_codes(sample)['publication']
                if fine_name == '':
                    fine_name = sample
                fine_sample_names.append(fine_name)
            except:
                fine_sample_names.append(sample)
        shape_geno_combined.index = fine_sample_names
        shape_geno_combined = shape_geno_combined.dropna(how='all')
        shape_geno_combined.to_csv(pjoin(shape_dir, 'shape_genotypes.csv'))

       # We have to keep one sample per accession
        for gene in shape_genotypes:
            repeated_accs = acc_info.data.groupby('accession').groups
            samples_to_remove = []
            for acc, samples in repeated_accs.items():
                if len(samples) > 1:
                    repeated_samples = [sample for sample in shape_genotypes[gene].data.index if sample in samples]
                    if repeated_samples > 1:
                        random.shuffle(samples)
                    samples_to_remove.extend(repeated_samples[1:])
            shape_genotypes[gene] = MatrixWithMeta(shape_genotypes[gene].data.drop(samples_to_remove),
                                                   shape_genotypes[gene].meta)

        #table of pairwise genes
        for sp in ['SP', 'SLC', 'SLL']:
            sp_acc = acc_info.data[acc_info.data['species'] == sp].index
            shape_genes = ['fw22', 'fw32', 'lc', 'ovate', 'fas', 'sun', 'sun_dup']
            shape_pairs = DataFrame(index=sp_acc, columns=shape_genes)
            counts = {}
            index_order = []
            for indx_gene1, gene1 in enumerate(shape_genes):
                if gene1 not in counts:
                    counts[gene1] = {}
                index_order.append(gene1)
                keep_accs = set(sp_acc).intersection(shape_genotypes[gene1].data.index)
                shape_gene1 = shape_genotypes[gene1].data.ix[list(keep_accs)]
                shape_pairs[gene1] = (shape_gene1['allele1'] == 'B') | (shape_gene1['allele2'] == 'B')

                if indx_gene1 == 0:
                    counts[gene1]['tot_counts'] = len(shape_gene1.index)
                    counts[gene1]['b_counts'] = sum(shape_pairs[gene1].dropna())

                for gene2 in shape_genes[indx_gene1 + 1:]:
                    pair = gene1 + '-' + gene2
                    index_order.append(pair)
                    if gene2 not in counts:
                        counts[gene2] = {}
                    counts[pair] = {}
                    keep_accs = set(sp_acc).intersection(shape_genotypes[gene2].data.index)
                    shape_gene2 = shape_genotypes[gene2].data.ix[list(keep_accs)]
                    shape_pairs[gene2] = (shape_gene2['allele1'] == 'B') | (shape_gene2['allele2'] == 'B')
                    counts[gene2]['tot_counts'] = len(shape_gene2.index)
                    counts[gene2]['b_counts'] = sum(shape_pairs[gene2].dropna())

                    counts[pair]['tot_counts'] = len(shape_pairs[[gene1, gene2]].dropna(how='any'))
                    shape_pair = shape_pairs[[gene1, gene2]].dropna(how='any')
                    counts[pair]['b_counts'] = sum((shape_pair[gene1] == True) & (shape_pair[gene2] == True))
                    b_gene1 = counts[gene1]['b_counts']
                    t_gene1 = counts[gene1]['tot_counts']
                    b_gene2 = counts[gene2]['b_counts']
                    t_gene2 = counts[gene2]['tot_counts']
                    b_pair = counts[pair]['b_counts']
                    t_pair = counts[pair]['tot_counts']
                    #counts[pair]['b_expected'] = b_gene1 / t_gene1 * b_gene2 / t_gene2 * t_pair

            fpath = pjoin(shape_dir, sp + '_pairwise_genes.csv')
            #DataFrame(counts, columns=index_order).transpose().to_csv(fpath)

        plot_shape_histogram(shape_dir, shape_genotypes, acc_info,
                             binomial_ci=True, min_num_alleles=12)

    if shape_table_by_sample:

        shape_dir = pjoin(RESULTS_DIR, 'shape_per_sample', )

        if not exists(shape_dir):
            os.mkdir(shape_dir)

        acc_syns = {}
        shape_genotypes = merge_shape_genotypes(acc_syns)

        genes = ['lc', 'fw22', 'fw32', 'fw113', 'fas', 'ovate', 'sun',
                 'sun_dup']

        combined_genos = {}
        for gene, genos in shape_genotypes.iteritems():
            if gene not in genes:
                continue
            combined_genos_values =[]
            combined_genos_index = []
            for acc, values in genos.data.iterrows():
                if values[0] == values[1]:
                    value = values[0]
                else:
                    value = values[0] + '/' + values[1]
                combined_genos_values.append(value)
                combined_genos_index.append(acc)
            combined_genos[gene] = Series(combined_genos_values, index=combined_genos_index)

        combined_genos['species'] = acc_info.data['species']
        combined_genos['group1'] = acc_info.data['group1']
        combined_genos['group2'] = acc_info.data['group2']

        combined_genos = DataFrame.from_dict(combined_genos)
        combined_genos.to_csv(pjoin(shape_dir, 'shape_genotypes.csv'))


    if shape_haplotypes:
        shape_haplo_dir = pjoin(RESULTS_DIR, 'shape_haplo', )

        if not exists(shape_haplo_dir):
            os.mkdir(shape_haplo_dir)

        flanking_snps_fpath = pjoin(shape_haplo_dir, 'flanking_snps1.csv')
        flanking_snps = load_dataframe_csv(flanking_snps_fpath, sep='\t',
                                        index_col=0)

        phase_dir_path = pjoin(RESULTS_DIR, 'PHASE')
        hap_type = 'individual_haplotypes'
        haplos = parse_fast_phase_results(phase_dir_path)[hap_type]

        shape_haplotypes = {}
        for marker in flanking_snps.data.index:
            chrom = flanking_snps.data.ix[marker]['chrom']
            snps = flanking_snps.data.ix[marker]
            snps = snps.dropna()
            snps = [snp for snp in snps if snp in haplos[chrom].data.columns]
            marker_hap = haplos[chrom].data[snps]
            marker_hap = marker_hap.apply(lambda x: ''.join(x[snps]), axis=1)
            marker_hap = marker_hap.unstack(level=1)
            shape_haplotypes[marker] = MatrixWithMeta(marker_hap)

        plot_shape_histogram(shape_haplo_dir, shape_haplotypes, acc_info,
                             binomial_ci=False, min_num_alleles=2)

    if plot_coancestry_matrix:

        coancestry_matrix = read_coancestry_matrix()

        species = ['SLC', 'SLL', 'SP']
        filter_condition = {'species': species}
#         coancestry_matrix = filter_matrix(coancestry_matrix, acc_info,
#                                           filter_condition, reverse=False)

        species = ['SP_SLC_inter_species_1', 'SLC_SP_inter_species', 'SLC_1',
                   'SP_1']
        filter_condition = {'group1': species}
#         coancestry_matrix = filter_matrix(coancestry_matrix, acc_info,
#                                           filter_condition, reverse=True)

        species = []
        filter_condition = {'group2': species}
#         coancestry_matrix = filter_matrix(coancestry_matrix, acc_info,
#                                           filter_condition, reverse=True)

        groupby = 'group2'

        plot_grouped_coancestry(coancestry_matrix, acc_info, groupby)
        plot_coancestry(coancestry_matrix, acc_info)

    if population_distance:

        distance_dir = pjoin(RESULTS_DIR, 'distance')
        if not exists(distance_dir):
            os.mkdir(distance_dir)

        genotype_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotype_accs)
        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)
        classification = acc_info.data['group2']

        dist = PopDistances(genotypes, classification,
                            conf_intervals_repeats=0, permutation_repeats=0,
                            min_num_individuals_per_pop=6)

        dist_matrix = dist.pairwise_distances('jost_dest')['distance_matrix']

        fhand = open(pjoin(distance_dir, 'group2_distance.csv'), 'w')
        dist_matrix.to_csv(fhand)
        fhand.close()

        boot_matrices = dist.jackknife_pairwise_distances('jost_dest',
                                                          num_repeats=1000)

        tree = calc_nj(dist_matrix, boot_matrices)
        tree_fpath = pjoin(distance_dir, 'all_tree_1000_6_indv_wo_mixture.newick')
        newick = tree.getNewick(with_distances=True, with_boot=True)
        open(tree_fpath, 'w').write(newick)

    if rarefaction:
        rarefaction_dir = pjoin(RESULTS_DIR, 'rarefaction', 'FM_PROC_VIN')
        if not exists(rarefaction_dir):
            os.mkdir(rarefaction_dir)

        genotype_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotype_accs)
        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)
        classification = acc_info.data['group1']

        poly = PopPolymorphism(genotypes, classification,
                               num_repeats=100)
        to_do = False
        to_plot = True
        median_fhand = pjoin(rarefaction_dir, 'median.csv')
        upper_fhand = pjoin(rarefaction_dir, 'upper_ci.csv')
        lower_fhand = pjoin(rarefaction_dir, 'lower_ci.csv')
        if to_do:
            rarefaction = poly.rarefaction()
            rarefaction['median'].to_csv(open(median_fhand, 'w'))
            rarefaction['upper_ci'].to_csv(open(upper_fhand, 'w'))
            rarefaction['lower_ci'].to_csv(open(lower_fhand, 'w'))
        if to_plot:
            hieral_classes = {}
            for acc_name, acc in acc_info.data.iterrows():
                species = acc['species']
                group1, group2 = acc['group1'], acc['group2']
                if species not in hieral_classes:
                    hieral_classes[species] = {}
                if group1 not in hieral_classes[species]:
                    hieral_classes[species][group1] = set()
                hieral_classes[species][group1].add(group2)
                if not group1.startswith(species):
                    print species, group1, group2, acc_name

            conf_intervals = False
            for species, groups1 in hieral_classes.viewitems():
                groups1 = groups1.keys()
                rarefaction = {}
                groups1 = ['Fresh_market', 'Processing', 'Vintage']
                rarefaction['median'] = read_csv(open(median_fhand),
                                                 index_col=0)
                rarefaction['median'] = rarefaction['median'][groups1]
                if conf_intervals:
                    rarefaction['upper_ci'] = read_csv(open(upper_fhand),
                                                       index_col=0)
                    rarefaction['lower_ci'] = read_csv(open(lower_fhand),
                                                       index_col=0)
                plot_groups = generate_plot_groups_for_rarefaction(
                                                         rarefaction['median'],
                                                         rarefaction['upper'],
                                                         rarefaction['lower'])
                fhand = open(pjoin(rarefaction_dir,
                                   'rarefaction' + species + '.svg'), 'w')
                scatter_groups(plot_groups, fhand, plot_lines=True)

    if allele_freq_distribution:
        allele_freq_dir = pjoin(RESULTS_DIR, 'allelic_freq')
        if not exists(allele_freq_dir):
            os.mkdir(allele_freq_dir)

        genotype_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotype_accs)
        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)

        classification = acc_info.data['group2']
        freqs = AlleleFreqsDistribs(genotypes, classification)
        hist = freqs.maf_histograms()
        fhand = open(pjoin(allele_freq_dir, 'hist_freqs.csv'), 'w')
        hist.to_csv(fhand)

        classification = acc_info.data['group2'].replace(['SLC_Peru_1_1',
                                                          'SP_Peru_1'],
                                                         'mixed')
        freqs = AlleleFreqsDistribs(genotypes, classification)
        hist = freqs.maf_histograms()
        fhand = open(pjoin(allele_freq_dir, 'hist_freqs_mixed.csv'), 'w')
        hist.to_csv(fhand)

    if kml:
        fpath = pjoin(RESULTS_DIR, 'SolCAPjoin.kml')
        write_kml_file(acc_info, fpath)

    grid_dir = pjoin(RESULTS_DIR, 'grid')
    n_bins = 5

    if stats_per_grid:
        #colombia and mesoamerica
        #limits = (-112, -73, 1.3, 30.5)
        #Peru
        #limits = (-82.0, -65, -18, -4)
        if not exists(grid_dir):
            os.mkdir(grid_dir)

        acc_info_geo = acc_info.data[['Longitude', 'Latitude']].replace('None',
                                                                        '')
        acc_info_geo['Longitude'] = acc_info_geo['Longitude'].map(lon_to_deg)
        acc_info_geo['Latitude'] = acc_info_geo['Latitude'].map(lat_to_deg)

        acc_info_geo = acc_info_geo.dropna()
        to_keep = set(acc_info_geo.index).intersection(genotypes.data.index)
        acc_info_geo = acc_info_geo.ix[list(to_keep)]
        gps = MatrixWithMeta(acc_info_geo)

        # andean and mesoamerica
        # lon_min=-112, lon_max=-65,
        # lat_min=-18, lat_max=30.5
        # andean
        # lon_min=-82, lon_max=-65,
        # lat_min=-18, lat_max=1.6
        #Ecuador and north peru
        # lon_min=-82, lon_max=-77.6,
        # lat_min=-6.5, lat_max=1.6

        bins = split_in_grids(gps, n_bins, lon_col='Longitude',
                              lat_col='Latitude',
                              lon_min=-82, lon_max=-77.6,
                              lat_min=-6.5, lat_max=1.6)
        # keep only the regions with enough accessions
        min_num_individuals = 5
        bins = {k: v for k, v in bins.items() if len(v['row_names']) > min_num_individuals}

        # transform to a classification series
        pops = []
        accs = []
        geo_locs = {}
        for index, (geo_pop, v) in enumerate(bins.items()):
            accs.extend(v['row_names'])
            pops.extend([index] * len(v['row_names']))
            geo_locs[index] = v['mid_point']
        geo_class = Series(pops, index=accs)
        geo_locs = DataFrame(geo_locs, index=['Longitude', 'Latitude'])

        to_keep = list(set(geo_class.index).union(genotypes.data.index))
        genos = genotypes.data.ix[to_keep]
        genos = genos.dropna()
        geo_class = geo_class[to_keep]

        genotypes = MatrixWithMeta(genos, genotypes.meta)

        tgenotypes = transpose_genotypes(genotypes)
        stats = calculate_basic_stats(tgenotypes, geo_class,
                                      min_num_individuals=5)
        stats.to_csv(pjoin(grid_dir, 'pop_stas.csv'))
        geo_locs.to_csv(pjoin(grid_dir, 'pop_locs.csv'))

    if plot_stats_per_grid:
        stats = load_dataframe_csv(open(pjoin(grid_dir, 'pop_stas.csv')),
                                   index_col=0, sep=',')
        stats = stats.data
        stats.columns = range(len(stats.columns))
        geo_locs = load_dataframe_csv(open(pjoin(grid_dir, 'pop_locs.csv')),
                                   index_col=0, sep=',')
        geo_locs = geo_locs.data
        geo_locs.columns = range(len(geo_locs.columns))

        # we create the interpolation to complete the grid with values
        obs_het = stats.xs('obs_het')
        exp_het = stats.xs('exp_het')
        x, y, obs, exp = [], [], [], []
        for pop, loc in geo_locs.iteritems():
            lon = loc['Longitude']
            x.append(lon)
            lat = loc['Latitude']
            y.append(lat)
            obs.append(obs_het[pop])
            exp.append(exp_het[pop])
        obs = interp2d(x, y, obs)
        exp = interp2d(x, y, exp)

        x_min, x_max = min(x), max(x)
        bin_width = (x_max - x_min) / n_bins
        xs = numpy.arange(x_min, x_max, bin_width)
        y_min, y_max = min(y), max(y)
        ys = numpy.arange(y_min, y_max, bin_width)

        iobs = []
        iexp = []
        for y in ys:
            iobs.append([obs(x, y)[0] for x in xs])
            iexp.append([exp(x, y)[0] for x in xs])

        from matplotlib.backends.backend_agg import FigureCanvasAgg as FCanvas
        from matplotlib.figure import Figure
        fig = Figure(figsize=(15.0, 11.0))
        canvas = FCanvas(fig)
        axes = fig.add_subplot(111)
        axes.imshow(iexp)
        canvas.print_figure(open(pjoin(grid_dir, 'exp_het.png'), 'w'))
        fig = Figure(figsize=(15.0, 11.0))
        canvas = FCanvas(fig)
        axes = fig.add_subplot(111)
        axes.imshow(iobs)
        canvas.print_figure(open(pjoin(grid_dir, 'obs_het.png'), 'w'))

    if adze:
        adze_dir = pjoin(RESULTS_DIR, 'adze')
        if not exists(adze_dir):
            os.mkdir(adze_dir)

        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)
        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)
        class_grp = 'Rarefaction_1'
        classification = acc_info.data.ix[accs_for_analysis][class_grp].copy()

        color_mapper = {'SLC_non_Andean': '#00aad4', 'SLC_Andean': '#000080',
                        'SP': '#55d400', 'SLL_vintage': '#ff0000',
                        'SLL_fresh': '#aa0088', 'SLL_processing': '#ff6600',
                        'SLL_contemporary': '#aa0088'}

        results = do_rarefaction(genotypes, classification)

        plots = ['private_alleles', 'allelic_richness']

        for plot in plots:
            plot_groups = generate_plot_groups_for_rarefaction(results[plot],
                                                     color_mapper=color_mapper)
            fhand = open(pjoin(adze_dir, plot + '.' + class_grp + '.svg'),
                         'w')
            scatter_groups(plot_groups, fhand, plot_lines=True)

    if pairwise_adze:

        p_adze_dir = pjoin(RESULTS_DIR, 'pairwise_adze')
        if not exists(p_adze_dir):
            os.mkdir(p_adze_dir)

        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)
        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)
        class_grp = 'group2'
        classification = acc_info.data.ix[accs_for_analysis][class_grp].copy()

        populations = classification.groupby(classification).groups

        pops = populations.keys()

        grouped_results = {}

        for i in range(len(pops)):
            pop1 = pops[i]
            for pop2 in pops[i + 1:]:
                pair = list(populations[pop1])
                pair.extend(populations[pop2])

                pair_classification = classification.data.ix[pair]
                pair_genotypes = MatrixWithMeta(genotypes.data.ix[pair],
                                                genotypes.meta)

                grouped_classification = Series(data='grouped',
                                               index=classification.data.index)
                grp_name = pop1 + '-' + pop2
                grouped_classification.ix[pair] = grp_name

                pair_results = do_rarefaction(pair_genotypes,
                                              pair_classification)
                group_results = do_rarefaction(genotypes,
                                                 grouped_classification,
                                                 standardized_sample_size=16)

                result = group_results['private_alleles'].ix[-1][grp_name]
                grouped_results[grp_name] = result

                plots = ['private_alleles', 'allelic_richness']

                for plot in plots:
                    plot_groups = generate_plot_groups_for_rarefaction(pair_results[plot])
        #                                                     color_mapper=color_mapper)
                    fhand = open(pjoin(p_adze_dir, plot + '.' + pop1 + '_' +
                                       pop2 + '.svg'), 'w')
                    scatter_groups(plot_groups, fhand, plot_lines=True)

        fhand = pjoin(p_adze_dir, 'private_alleles.csv')
        Series(grouped_results).to_csv(fhand)

    if do_structure:

        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)
        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)

        interpoled_map_fhand = open(INTERPOLED_RECOM_MAP)
        interpolated_gen_map = load_markers_map(interpoled_map_fhand,
                                                molecule_col='chromosome',
                                                location_col='genetic_position',
                                                linkage_group_col='chromosome',
                                       genetic_location_col='genetic_position',
                                                index_col=0, sep='\t')

        markers_in_map = interpolated_gen_map.data.index
        markers_in_both = set(markers_in_map).intersection(genotypes.data.columns)
        markers_in_both = list(markers_in_both)
        genotypes = MatrixWithMeta(genotypes.data[markers_in_both],
                                   genotypes.meta)
        structure_dir = pjoin(RESULTS_DIR, 'structure')
        if not exists(structure_dir):
            os.mkdir(structure_dir)

        fpath = pjoin(structure_dir, 'snp.geno')

        mainparams, extraparams = get_default_params()
        mainparams = create_genotype_file(genotypes, fhand=open(fpath, 'w'),
                                          params=mainparams,
                                          markers_map=interpolated_gen_map)[1]

        main_fhand = open(pjoin(structure_dir, 'mainparams'), 'w')
        mainparams_fhand, extraparams_fhand = create_param_files(mainparams,
                                                                 extraparams,
                                                    mainparams_fhand=main_fhand)
        return
        do_structure_for_ks(2, 3, genotype_fhand=open(fpath),
                            directory=structure_dir,
                            mainparams_fhand=open(mainparams_fhand.name),
                            extraparams_fhand=open(extraparams_fhand.name),
                            in_parallel=True)

    if structure_results:
        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)

        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)
        acc_info = MatrixWithMeta(acc_info.data.ix[accs_for_analysis],
                                  acc_info.meta)
        acc_info.meta[CLASSIFICATION_COL] = 'group2'

        structure_dir = pjoin(RESULTS_DIR, 'structure')

        mainparams, extraparams = get_default_params()
        geno_fhand, mainparams = create_genotype_file(genotypes,
                                                      params=mainparams)

        results = read_structure_results(structure_dir,
                                         genotype_fhand=open(pjoin(structure_dir, 'snp.geno')))
        ln_prob_fpath = os.path.join(structure_dir, 'ln_probs.svg')
        ln_prob_fhand = open(ln_prob_fpath, 'w')
        scatter_series(Series(results['ks']), Series(results['ln_probs']),
                       ln_prob_fhand)
        # classifications
        for k_val in range(2, 30):
            k_index = k_val - 2
            indi_class = results['individual_classifications'][k_index]

            pops = GROUPING_ORDER['group2']
            # ancestries stacked bars
            indi_class = sort_individual_classification(indi_class,
                                                        acc_info,
                                                        pops)

            fhand = open(os.path.join(structure_dir,
                                     'indi_ancestries.' + str(k_val) + '.svg'),
                         'w')
            classification = acc_info.data.ix[indi_class.data.index]['group2']
            stacked_bars(indi_class.data, fhand, classification=classification)

    if dist_distributions:
        dist_dir = pjoin(RESULTS_DIR, 'dist_distributions')
        if not exists(dist_dir):
            os.mkdir(dist_dir)
        fpath = pjoin(dist_dir, 'dists_per_group.csv')
        fhand = open(fpath, 'w')
        calculate_dist_distributions(genotypes, acc_info, fhand)

    comparision = 'SLC'

    if indiv_genetic_tree:

        dist_dir = pjoin(RESULTS_DIR, 'indiv_genetic_tree')
        if not exists(dist_dir):
            os.mkdir(dist_dir)

        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)
        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                    genotypes.meta)
        acc_info = MatrixWithMeta(acc_info.data.ix[accs_for_analysis],
                                  acc_info.meta)

        for grouping in ['species', 'group1', 'group2']:
            fhand = open(pjoin(dist_dir, 'hypertree_colors_' + grouping + '.txt'), 'w')
            for ind in genotypes.data.index:
                group = acc_info.data.loc[ind, grouping]
                color = color_definitions[grouping][group]
                rgb_color = html_to_rgb_color(color)
                fhand.write('%s %s,%s,%s\n' % (ind, rgb_color[0], rgb_color[1],
                                                  rgb_color[2]))
            fhand.close()

        return

        dist_matrix = pairwise_individual_genetic_distance(genotypes)

        fhand = open(pjoin(dist_dir, 'genetic_distances.csv'), 'w')
        dist_matrix.to_csv(fhand)

#         dist_matrix = load_dataframe_csv(open(pjoin(dist_dir, 'genetic_distances.csv'))).data
        tree = calc_nj(dist_matrix)
        tree_fpath = pjoin(dist_dir, 'tree.newick')
        newick = tree.getNewick(with_distances=True, with_boot=False)
        open(tree_fpath, 'w').write(newick)

    if isolation_by_distance:

        dist_dir = pjoin(RESULTS_DIR, 'distances_correlation', 'geographic',
                         comparision)
        if not exists(dist_dir):
            os.makedirs(dist_dir)

        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)
        _genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                    genotypes.meta)
        _acc_info = MatrixWithMeta(acc_info.data.ix[accs_for_analysis],
                                  acc_info.meta)

        individuals = _genotypes.data.index
        cols = ['accession', 'Longitude', 'Latitude']
        acc_info_geo = _acc_info.data.ix[individuals][cols].replace('None', '')
        acc_info_geo['Longitude'] = acc_info_geo['Longitude'].map(lon_to_deg)
        acc_info_geo['Latitude'] = acc_info_geo['Latitude'].map(lat_to_deg)
        acc_info_geo = acc_info_geo.dropna()
        acc_info_geo = MatrixWithMeta(acc_info_geo)

        ind_analysis = acc_info_geo.data.index
        genotypes_with_geo = MatrixWithMeta(_genotypes.data.ix[ind_analysis],
                                            genotypes.meta)

        geo_dist = pairwise_geographic_distance(acc_info_geo,
                                                longitude_col='Longitude',
                                                latitude_col='Latitude')

        gen_dist = pairwise_individual_genetic_distance(genotypes_with_geo)

        fhand = open(pjoin(dist_dir, 'correlation_results.txt'), 'w')
        correlation_result = correlation_distance_matrices(geo_dist, gen_dist,
                                                     mantel_permutations=10000)
        fhand.write('genetic vs. geographic distance: ' + comparision + '\n')
        fhand.write('slope: %s\n' % correlation_result['slope'])
        fhand.write('intercept: %s\n' % correlation_result['intercept'])
        fhand.write('r: %s\n' % correlation_result['r'])
        fhand.write('mantel_p-value: %s\n' % correlation_result['mantel_p-value'])
        fhand.close()

        fhand = open(pjoin(dist_dir, 'genetic_geographic.txt'), 'w')
        write_distances_output(fhand, gen_dist, geo_dist)

    if climatic_correlation:

        clim_dir = pjoin(RESULTS_DIR, 'distances_correlation', 'climatic',
                         comparision)
        if not exists(clim_dir):
            os.makedirs(clim_dir)

        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)
        _genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)
        _acc_info = MatrixWithMeta(acc_info.data.ix[accs_for_analysis],
                                  acc_info.meta)

        fpath = pjoin(RESULTS_DIR, 'distances_correlation', 'climatic',
                      'climatic_data.csv')
        climatic_data = load_dataframe_csv(fpath, index_col=0).data

        accs_keep = set(_acc_info.data['accession']).intersection(climatic_data.index)
        climatic_data = climatic_data.ix[accs_keep]

        samples_keep = [sample for sample in _acc_info.data.index if _acc_info.data.get_value(sample, 'accession') in accs_keep]
        acc_info_data = acc_info.data.ix[samples_keep]
        climatic_data = climatic_data.reindex(acc_info_data['accession'])

        pca = do_pca(MatrixWithMeta(climatic_data.iloc[:, 0:20]))

        clim_dist = pdist(pca['projections'].data.values)
        clim_dist = squareform(clim_dist, force='tomatrix')
        clim_dist = DataFrame(clim_dist, index=climatic_data.index,
                              columns=climatic_data.index)

        clim_genotypes = MatrixWithMeta(_genotypes.data.ix[acc_info_data.index],
                                        genotypes.meta)
        gen_dist = pairwise_individual_genetic_distance(clim_genotypes)
        gen_dist.index = clim_dist.index
        gen_dist.columns = clim_dist.columns

        fhand = open(pjoin(clim_dir, 'correlation_results.txt'), 'w')
        correlation_result = correlation_distance_matrices(clim_dist, gen_dist,
                                                     mantel_permutations=10000)
        fhand.write('genetic vs. climatic distance: ' + comparision + '\n')
        fhand.write('slope: %s\n' % correlation_result['slope'])
        fhand.write('intercept: %s\n' % correlation_result['intercept'])
        fhand.write('r: %s\n' % correlation_result['r'])
        fhand.write('mantel_p-value: %s\n' % correlation_result['mantel_p-value'])
        fhand.close()

        fhand.close()

        fhand = open(pjoin(clim_dir, 'genetic_climatic.txt'), 'w')
        write_distances_output(fhand, gen_dist, clim_dist)

    if climatic_geographic_correlation:

        dist_dir = pjoin(RESULTS_DIR, 'distances_correlation',
                         'climatic_geographic', comparision)
        if not exists(dist_dir):
            os.makedirs(dist_dir)

        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)

        acc_info = MatrixWithMeta(acc_info.data.ix[accs_for_analysis],
                                  acc_info.meta)

        fpath = pjoin(RESULTS_DIR, 'distances_correlation', 'climatic',
                      'climatic_data.csv')
        climatic_data = load_dataframe_csv(fpath, index_col=0).data

        accs_keep = set(acc_info.data['accession']).intersection(climatic_data.index)
        climatic_data = climatic_data.ix[accs_keep]
        samples_keep = [sample for sample in acc_info.data.index if acc_info.data.get_value(sample, 'accession') in accs_keep]

        cols = ['accession', 'Longitude', 'Latitude']
        acc_info_geo = acc_info.data.ix[samples_keep][cols].replace('None', '')
        acc_info_geo['Longitude'] = acc_info_geo['Longitude'].map(lon_to_deg)
        acc_info_geo['Latitude'] = acc_info_geo['Latitude'].map(lat_to_deg)
        acc_info_geo = acc_info_geo.dropna()
        acc_info_geo = acc_info_geo.reindex(acc_info_geo['accession'])
        acc_info_geo = MatrixWithMeta(acc_info_geo)

        geo_dist = pairwise_geographic_distance(acc_info_geo,
                                                longitude_col='Longitude',
                                                latitude_col='Latitude')

        acc_info_data = acc_info.data.ix[samples_keep]
        climatic_data = climatic_data.reindex(acc_info_data['accession'])

        pca = do_pca(MatrixWithMeta(climatic_data.iloc[:, 0:20]))

        clim_dist = pdist(pca['projections'].data.values)
        clim_dist = squareform(clim_dist, force='tomatrix')
        clim_dist = DataFrame(clim_dist, index=climatic_data.index,
                              columns=climatic_data.index)

        fhand = open(pjoin(dist_dir, 'correlation_results.txt'), 'w')
        correlation_result = correlation_distance_matrices(clim_dist, geo_dist,
                                                     mantel_permutations=10000)
        fhand.write('genetic vs. geographic distance: ' + comparision + '\n')
        fhand.write('slope: %s\n' % correlation_result['slope'])
        fhand.write('intercept: %s\n' % correlation_result['intercept'])
        fhand.write('r: %s\n' % correlation_result['r'])
        fhand.write('mantel_p-value: %s\n' % correlation_result['mantel_p-value'])
        fhand.close()

        fhand.close()

        fhand = open(pjoin(dist_dir, 'geographic_climatic.txt'), 'w')
        write_distances_output(fhand, geo_dist, clim_dist)

    if arlequin:

        arlequin_dir = pjoin(RESULTS_DIR, 'arlequin')
        if not exists(arlequin_dir):
            os.mkdir(arlequin_dir)

        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)
        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)
        acc_info = MatrixWithMeta(acc_info.data.ix[accs_for_analysis],
                                  acc_info.meta)
        fpath = pjoin(arlequin_dir, 'arlequin_input.arp')
        write_arlequin(fpath, genotypes, acc_info.data['group1'])

        return
        for chrm in range(1, 13):
            chrm_markers = physical_map.data[physical_map.data['chromosome'] == chrm].index
            chrm_markers = set(chrm_markers).intersection(genotypes.data.columns)
            chrm_markers = list(chrm_markers)
            genotypes = MatrixWithMeta(genotypes.data[chrm_markers],
                                       genotypes.meta)
            fpath = pjoin(arlequin_dir, 'arlequin_input.chrom' + str(chrm) + '.arp')
            write_arlequin(fpath, genotypes, acc_info.data['species'])

    if snapp:

        snapp_dir = pjoin(RESULTS_DIR, 'snapp')
        if not exists(snapp_dir):
            os.mkdir(snapp_dir)

        fhand = open(pjoin(snapp_dir, 'accs_for_snapp.csv'))

        samples_for_acc = load_dataframe_csv(fhand, index_col=0, sep='\t')
        samples_for_acc = list(samples_for_acc.data['sample'])

        print set(samples_for_acc).intersection(genotypes.data.index)

        genotypes = MatrixWithMeta(genotypes.data.ix[samples_for_acc],
                                   genotypes.meta)
        fhand = pjoin(snapp_dir, 'snapp.nex')
        write_nexus_file(fhand, genotypes)

    if freq_by_locus:
        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)
        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)
        acc_info = MatrixWithMeta(acc_info.data.ix[accs_for_analysis],
                                  acc_info.meta)

        class_criteria = 'group2'

        freqs = AlleleFreqsDistribs(genotypes, acc_info.data[class_criteria])
        fhand = pjoin(RESULTS_DIR, 'frecuencias_locus.csv')
        freqs_by_locus = freqs.freq_by_locus()

        index = [pop for pop in GROUPING_ORDER[class_criteria] if pop in freqs_by_locus.index]

        freqs_by_locus.ix[index].T.to_csv(fhand, sep='\t')

    if genalex:
        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)
        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)
        acc_info = MatrixWithMeta(acc_info.data.ix[accs_for_analysis],
                                  acc_info.meta)

        fhand = open(pjoin(RESULTS_DIR, 'genalex.xls'), 'w')

        grouping_criteria = 'group2'

        groups = acc_info.data.groupby(grouping_criteria).groups
        codec = get_codec_from_genotypes(genotypes)
        n_markers = len(genotypes.data.columns)
        n_indis = len(genotypes.data.index)
        n_pops = len(groups)
        n_indis_pop = []
        accs_pop = []
        name_pops = []
        for grp, accs in groups.viewitems():
            n_indis_pop.append(len(accs))
            accs_pop.extend(accs)
            name_pops.append(grp)

        line = '%s,%s,%s,' % (n_markers, n_indis, n_pops)
        line += ','.join(map(str, n_indis_pop)) + '\n'
        line += ',,,' + ','.join(name_pops) + '\n'
        line += 'CODE,SITE,' + ',,'.join(genotypes.data.columns) + '\n'
        fhand.write(line)

        decoded_genos = genotypes.data.applymap(codec.decode_to_ints)
        for acc in accs_pop:
            pop = acc_info.data.loc[acc, grouping_criteria]
            line = acc + ',' + pop + ','
            geno = ['%s,%s' % (marker) if marker is not None else '0,0' for marker in decoded_genos.ix[acc]]
            line += ','.join(geno) + '\n'
            fhand.write(line)

        fhand.close()

    if genepop_file:
        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)
        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)
        acc_info = MatrixWithMeta(acc_info.data.ix[accs_for_analysis],
                                  acc_info.meta)

        fhand = open(pjoin(RESULTS_DIR, 'genepop.txt'), 'w')

        grouping_criteria = 'group2'

        groups = acc_info.data.groupby(grouping_criteria).groups
        codec = get_codec_from_genotypes(genotypes)

        line = 'Genepop file\n'
        line += '\n'.join(genotypes.data.columns) + '\n'
        fhand.write(line)

        decoded_genos = genotypes.data.applymap(codec.decode_to_ints)
        for grp, accs in groups.viewitems():
            fhand.write('Pop\n')
            for acc in accs:
                line = grp + ', '
                geno = ['0%s0%s' % (marker) if marker is not None else '0000' for marker in decoded_genos.ix[acc]]
                line += ' '.join(geno) + '\n'
                fhand.write(line)

        fhand.close()

    if bayescan_file:
        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)
        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)
        acc_info = MatrixWithMeta(acc_info.data.ix[accs_for_analysis],
                                  acc_info.meta)

        popcounts = _PopCounts(genotypes=genotypes,
                               classification=acc_info.data['group2'])

        pops = popcounts._genotypic_counts['pop_order']
        n_genotypes = popcounts._genotypic_counts['n_genotypes']
        allele_counts = {allele: popcounts._genotypic_counts['allele_freqs'][allele] * n_genotypes * 2 for allele in range(1,5)}
        reference_allele = [None] * len(genotypes.data.columns)
        for marker_index, ref in enumerate(reference_allele):
            for allele in range(1, 5):
                if any(allele_counts[allele][:, marker_index] != 0):
                    reference_allele[marker_index] = allele
                    break

        fhand = open(pjoin(RESULTS_DIR, 'bayescan.txt'), 'w')

        fhand.write('[loci]=' + '%s' % (len(genotypes.data.columns)) + '\n\n')
        fhand.write('[populations]=' + '%s' % (len(pops)) + '\n\n')

        n_alleles = genotypes.meta[PLOIDY]
        for index_pop, pop in enumerate(pops):
            fhand.write('[pop]=%s\n' % (index_pop + 1))
            for marker_index, ref_allele in enumerate(reference_allele):
                n_marker_genos = int(n_genotypes[index_pop, marker_index])
                line = '%s %s %s ' % (marker_index + 1, n_marker_genos * 2,
                                      n_alleles)
                ref_counts = int(allele_counts[ref_allele][index_pop, marker_index])
                line += '%s %s\n' % (ref_counts, n_marker_genos * n_alleles - ref_counts)
                fhand.write(line)

        fhand.close()

    if tassel_file:
        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)
        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)
        acc_info = MatrixWithMeta(acc_info.data.ix[accs_for_analysis],
                                  acc_info.meta)

        fhand = open(pjoin(RESULTS_DIR, 'tessel_genotypes.txt'), 'w')

        line = 'rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\t'
        line += 'protLSID\tassayLSID\tpanel\tQCcode\t'
        line += '\t'.join(genotypes.data.index) + '\n'
        fhand.write(line)

        codec = get_codec_from_genotypes(genotypes)

        def _decode_to_tassel(genotype):
            decode_geno = codec.decode_to_genotype(genotype)
            if decode_geno is not None:
                decode_geno = '%s%s' % (decode_geno[0], decode_geno[1])
            else:
                decode_geno = 'NN'
            return decode_geno

        decoded_genotypes = genotypes.data.applymap(_decode_to_tassel)

        for marker in genotypes.data.columns:
            chrm = physical_map.data.loc[marker, 'chromosome']
            if chrm == 13:
                chrm = 12
            pos = physical_map.data.loc[marker, 'position']
            line = '%s\tNA\t%.0f\t%.0f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t' % (marker, chrm, pos)
            line += '\t'.join(decoded_genotypes.loc[:, marker]) + '\n'
            fhand.write(line)
        fhand.close()

    if LD:
        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)
        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)
        acc_info = MatrixWithMeta(acc_info.data.ix[accs_for_analysis],
                                  acc_info.meta)

        results = calculate_unphased_LD(genotypes, physical_map,
                                        acc_info.data['Rarefaction_1'],
                                        method='Rogers-Huff')
        fpath = pjoin(RESULTS_DIR, 'LD.csv')
        results.to_csv(fpath)

    if tassel_file_from_phased:

        genotypes_accs = genotypes.data.index
        acc_info_accs = acc_info.data.index
        accs_for_analysis = set(acc_info_accs).intersection(genotypes_accs)
        genotypes = MatrixWithMeta(genotypes.data.ix[accs_for_analysis],
                                   genotypes.meta)
        acc_info = MatrixWithMeta(acc_info.data.ix[accs_for_analysis],
                                  acc_info.meta)

        groups = acc_info.data.groupby(acc_info.data['Rarefaction_1']).groups

        for group, grp_accs in groups.viewitems():

            phase_dir_path = pjoin(RESULTS_DIR, 'PHASE')
            hap_type = 'individual_haplotypes'
            haplotypes = parse_fast_phase_results(phase_dir_path)[hap_type]

            fhand = open(pjoin(RESULTS_DIR, 'tassel_genotypes_' + group + '.txt'), 'w')
            line = 'rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\t'
            line += 'protLSID\tassayLSID\tpanel\tQCcode\t'

            indis = grp_accs

            line += '\t'.join(indis) + '\n'
            fhand.write(line)

            codec = get_codec_from_genotypes(genotypes)

            valid_markers = genotypes.data.columns

            chroms = ['chrom1', 'chrom2', 'chrom3', 'chrom4', 'chrom5',
                      'chrom6', 'chrom7', 'chrom8', 'chrom9', 'chrom10',
                      'chrom11', 'chrom12', 'chrom13']
            for chrom in chroms:
                chrom_haplotypes = haplotypes[chrom].data.ix[indis]
                valid_chrom_markers = [marker for marker in chrom_haplotypes.columns if marker in valid_markers]
                for marker in valid_chrom_markers:
                    m_chrom = chrom
                    marker_genos = chrom_haplotypes[marker]
                    marker_genos = marker_genos.unstack(level=1)
                    decoded_marker_haplo = marker_genos.haplotype1 + marker_genos.haplotype2

                    pos = physical_map.data.loc[marker, 'genetic_position'] * 100000

                    m_chrom = m_chrom.replace('chrom', '')
                    if m_chrom == '13':
                        m_chrom = '12'
                        if pos < 3950000:
                            pos = 5600000
                    line = '%s\tNA\t%s\t%.0f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t' % (marker, m_chrom, pos)
                    line += '\t'.join(decoded_marker_haplo) + '\n'
                    fhand.write(line)
            fhand.close()


    fine_sample_names = []
    fine_accession_names = []

    cols = ['Longitude', 'Latitude']
    acc_info_geo = acc_info.data[cols].replace('None', '')
    acc_info.data['Longitude'] = acc_info_geo['Longitude'].map(lon_to_deg)
    acc_info.data['Latitude'] = acc_info_geo['Latitude'].map(lat_to_deg)

    return

    for sample in acc_info.data.index:
        try:
            fine_name = get_comav_acc_codes(sample)['publication']
            if fine_name == '':
                fine_name = sample
            fine_sample_names.append(fine_name)
            fine_accession_names.append(fine_name)
        except:
            fine_sample_names.append(sample)
            fine_accession_names.append(acc_info.data.ix[sample]['accession'])

    acc_info.data.index = fine_sample_names
    acc_info.data['accession'] = fine_accession_names
    cols = ['accession', 'original_bank', 'source_bank', 'species', 'group1',
            'group2', 'Passport_Species', 'Passport_Classification', 'Country',
            'Geographic_Region', 'Latitude', 'Longitude', 'collection_source',
            'percent_Het']
    acc_info.data[cols].to_csv(pjoin(RESULTS_DIR,'fine_clasification.csv'), sep='\t')

    fine_sample_names = []
    for sample in genotypes.data.index:
        try:
            fine_name = get_comav_acc_codes(sample)['publication']
            if fine_name == '':
                fine_name = sample
            fine_sample_names.append(fine_name)
        except:
            fine_sample_names.append(sample)
    genotypes.data.index = fine_sample_names
    codec = get_codec_from_genotypes(genotypes)
    dgenotypes = genotypes.data.applymap(codec.decode_to_genotype)
    dgenotypes.to_csv(pjoin(RESULTS_DIR,'genotypes.csv'))


if __name__ == '__main__':
    main()
