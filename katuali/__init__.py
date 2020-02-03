import argparse
import collections
import itertools
import logging
import operator
import os
import pathlib
import platform
import re
import shutil
import sys
import yaml

import gpustat
import pkg_resources
import pysam


__version__ = "0.3.0"
__pkg__ = __name__


Unit = collections.namedtuple('Unit', ('symbol', 'unit'))


def check_file_exists(fp, log_level=logging.INFO):
    """Check file exists, recursively following symlinks and forcing NFS cache updates with chown

        ... Note that chown changes owner and group to the present user.

    :param fp: str, filepath.
    :returns: str, real path of file, having followed all symlinks.
    """

    logging.basicConfig(
        format='[%(asctime)s - %(name)s] %(message)s',
        datefmt='%H:%M:%S', level=log_level)
    logger = logging.getLogger('check_file')

    logger.debug('Checking file on {}: {}'.format(platform.node(), fp))

    def _is_link_or_exists(fp):
        return os.path.islink(fp) or os.path.exists(fp)

    if not _is_link_or_exists(fp):
        logger.debug('File not found, forcing NFS cache update for {}'.format(fp))
        # force NFS cache update by changing owner and group to current user
        os.chown(fp, os.getuid(), os.getgid())
        if not _is_link_or_exists(fp):
             raise IOError('File not present even after chown: {}'.format(os.path.abspath(fp)))

    if os.path.islink(fp):  # recursively follow links
        logger.debug('File is symlink, following link. File: {}'.format(fp))
        # support links to absolute and relative paths
        target_path = os.readlink(fp)
        if not os.path.isabs(target_path):
            target_path = os.path.join(os.path.dirname(fp), target_path)
            logger.debug('File is relative symlink: {}'.format(target_path))
        # support links to absolute and relative paths
        return check_file_exists(target_path)
    else:
        logger.debug('File exists! File: {} Size: {}'.format(fp, os.path.getsize(fp)))
        if os.path.basename(fp) in {'basecalls.fasta', 'consensus.fasta'}:
            with pysam.FastxFile(fp) as fx:
                first_rec_name = next(fx).name
                logger.debug('First fastx record: {}'.format(first_rec_name))
        return fp


def _log_level():
    """Parser to set logging level and acquire software version/commit"""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

    #parser.add_argument('--version', action='version', version=get_version())

    modify_log_level = parser.add_mutually_exclusive_group()
    modify_log_level.add_argument('--debug', action='store_const',
        dest='log_level', const=logging.DEBUG, default=logging.INFO,
        help='Verbose logging of debug information.')
    modify_log_level.add_argument('--quiet', action='store_const',
        dest='log_level', const=logging.WARNING, default=logging.INFO,
        help='Minimal logging; warnings only).')

    return parser


def check_files_exist():
    parser = argparse.ArgumentParser(
        description='Check files exist, recursively following symlinks and forcing NFS cache updates with chown.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[_log_level()],
    )
    parser.add_argument('filepaths', nargs='+', help='Filepaths.')
    args=parser.parse_args()

    logging.basicConfig(
        format='[%(asctime)s - %(name)s] %(message)s',
        datefmt='%H:%M:%S', level=args.log_level)
    logger = logging.getLogger('check_files')

    for fp in args.filepaths:
        check_file_exists(fp, args.log_level)
    logger.debug('Finished checking that input files exist.')


def _data_path(filename):
    # entry point to locate datafile
    resource = os.path.join('data', filename)
    filepath = pkg_resources.resource_filename(__pkg__, resource)
    if not os.path.isfile(filepath):
        filepath = None
    return filepath


def print_data_path():
    filename = _data_path(sys.argv[1])
    if filename is not None:
        print(filename)
    else:
        sys.exit(1)


def create_config():
    parser = argparse.ArgumentParser(
        description='Create a template configuration file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('output', help='Output config file.')
    args = parser.parse_args()

    default = _data_path('config.yaml')
    shutil.copyfile(default, args.output)


def process_katuali_config():
    # Helper entry point for katuali shell wrapper
    parser = argparse.ArgumentParser(
        description='Helper entry point for katuali shell wrapper.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('inputs', nargs='+', help='Inputs: base_config out_config additional_config')

    args=parser.parse_args()
    args.base_configfile = args.inputs[0]
    args.out_configfile = args.inputs[1]
    args.additional_config = args.inputs[2] if len(args.inputs) > 2 else []
    if isinstance(args.additional_config, str):
        args.additional_config = [args.additional_config]

    # sys.argv ignores quotes, so would split e.g.
    # MINI_ASSEMBLE_OPTS="-n 10" GUPPY_OPTS="--hp_correct 1"
    # into
    # ['MINI_ASSEMBLE_OPTS=-n', '10', 'GUPPY_OPTS=--hp_correct', '1']

    fixed = []
    for c in args.additional_config:
        if '=' in c:
            fixed.append(c)
        else:
            fixed[-1] += ' ' + c

    d = {}
    for arg in fixed:
        k, v = arg.split('=')
        d[k] = v

    conf = yaml.load(open(args.base_configfile))
    conf.update(d)
    yaml.dump(conf, open(args.out_configfile, 'w'))


def pick_gpu():
    parser = argparse.ArgumentParser(
        description='Choose GPU to use from enviroment variable, falling back '
        'to GPU with lowest current use.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--env_var', default='SGE_HGR_gpu',
        help='Environment variable from which to extract GPU.')
    args=parser.parse_args()


    logging.basicConfig(
        format='[%(asctime)s - %(name)s] %(message)s',
        datefmt='%H:%M:%S', level=logging.INFO)
    logger = logging.getLogger('pick_gpu')

    gpu = os.getenv(args.env_var)
    if gpu is None:
        stats = gpustat.GPUStatCollection.new_query()
        sorter = lambda s: (s.memory_used, s.utilization, s.temperature)
        gpu = sorted(stats.gpus, key=sorter)[0].index
        logger.info('{} was not set, setting GPU to {} based on memory and utilization'.format(args.env_var, gpu))
    else:
        gpu = gpu.replace('cuda', '')
        logger.info('Using gpu {} from {}'.format(gpu, args.env_var))
    print(gpu)


class SafeDict(dict):
     def __missing__(self, key):
        return '{' + key + '}'


def partial_format(s, **kwargs):
    """Partially format a string leaving unused named placeholder unchanged

    >>> partial_format('{a}/{b}', b=1)
    '{a}/1'
    """
    return s.format_map(SafeDict(**kwargs))


def product_dict(**kwargs):
    """Given a dictionary of iterables, yield a dictionary for each combination in the product of all the iterables.

    Strings are not treated as iterables.

    >>> list(product_dict(**{'a': [25, 50], 'b': [1, 2], 'c': 'foo'}))
    [{'a': 25, 'b': 1, 'c': 'foo'},
     {'a': 25, 'b': 2, 'c': 'foo'},
     {'a': 50, 'b': 1, 'c': 'foo'},
     {'a': 50, 'b': 2, 'c': 'foo'}]
    """
    keys, values = zip(*kwargs.items())
    vals = [[v] if isinstance(v, str) or not hasattr(v, '__iter__')
             else v for v in values]
    for value_product in itertools.product(*vals):
        yield dict(zip(keys, value_product))


def int_to_formatted_string(i, fmt='{:.1f}'):
    """Convert int to formatted string with either k, m or g suffix

    >>> int_to_formatted_string(1.5e3)
    '1.5k'
    >>> int_to_formatted_string(4.86e6)
    '4.9m'
    >>> int_to_formatted_string(3e9)
    '3.0g'
    >>> int_to_formatted_string(1e12)
    '1000.0g'
    """
    units = (Unit('k', 1e3), Unit('M', 1e6), Unit('G', 1e9))
    units = sorted(units, key=operator.attrgetter('unit'))
    max_unit = max(u.unit for u in units)
    for u in units:
        scaled = round(i / u.unit, 1)
        if scaled < 100 or u.unit == max_unit:
            return fmt.format(scaled) + u.symbol


def get_region_len(region, ref_fasta):
    """Get region length from samtools region string, using start and end if present, else obtaining contig length from reference fasta.
    """
    if ':' in region:
        if not '-' in region:
            raise ValueError('Regions must be specified just as contig or contig:start-end')
        start, end = map(int, region.split(':')[1].split('-'))
        region_len = end -start
        if region_len < 1:
            raise ValueError('Region length < 1, check your region specification')
    else:
        # we have a full contig
        with pysam.FastaFile(ref_fasta) as fa:
            lengths = dict(zip(fa.references, fa.lengths))
        if region not in fa.references:
            raise KeyError('Region {} is not a contig in the reference'.format(region))
        region_len = lengths[region]

    return region_len


def expand_target_template(template, config):
    """Create target strings by formatting the template with all combinations of referenced config variables.

    Similar to snakemake's expand, but fills with config variables rather than wildcards.
    Also handles variables defined per dataset and the special variable GENOME_SIZE.

    >>> config = {
        'DATA': {'run1': {'REGIONS': ['a1', 'a1'], 'GENOME_SIZE': '4.8m'},
                 'run2': {'REGIONS': ['b1', 'b2'], 'GENOME_SIZE': '10m'}},
        'DEPTH': [25, 50],
        'BASECALLER': 'guppy'
        }
    >>> template = '{DATA}/basecall/guppy/align/{REGIONS}/{DEPTH}X/basecalls.fasta'
    >>> expand_target_template(template, config)
    ['run2/basecall/guppy/align/b2/25X/basecalls.fasta',
     'run2/basecall/guppy/align/b1/25X/basecalls.fasta',
     'run1/basecall/guppy/align/a1/25X/basecalls.fasta',
     'run2/basecall/guppy/align/b2/25X/basecalls.fasta',
     'run2/basecall/guppy/align/b1/25X/basecalls.fasta',
     'run1/basecall/guppy/align/a1/25X/basecalls.fasta',
     'run2/basecall/guppy/align/b2/50X/basecalls.fasta',
     'run2/basecall/guppy/align/b1/50X/basecalls.fasta',
     'run1/basecall/guppy/align/a1/50X/basecalls.fasta',
     'run2/basecall/guppy/align/b2/50X/basecalls.fasta',
     'run2/basecall/guppy/align/b1/50X/basecalls.fasta',
     'run1/basecall/guppy/align/a1/50X/basecalls.fasta']
    """

    to_expand = set(re.findall('\{(.*?)\}', template))
    datasets = tuple(config['DATA'])
    # get parameters which are coupled to the dataset and are in the template
    # e.g. not all genomes/regions are present in all datasets/runs
    dataset_params = set(itertools.chain(*[config['DATA'][v] for v in config['DATA']]))
    dataset_params = dataset_params & to_expand
    # get parameters which are not coupled to the dataset
    non_dataset_params = to_expand - dataset_params - {'GENOME_SIZE'}
    non_dataset_params = {k: config[k] for k in non_dataset_params if k != 'DATA'}

    templates = []
    for dataset in datasets:
        # use partial_format to format only the DATA curly brace in the template
        dataset_tmp = partial_format(template, DATA=dataset)
        # handle genome size if it's in the template
        if 'GENOME_SIZE' in to_expand:
            if 'GENOME_SIZE' in config['DATA'][dataset]:
                # use single genome size defined in config
                genome_sz = config['DATA'][dataset]['GENOME_SIZE']
                dataset_templates = [partial_format(dataset_tmp, GENOME_SIZE=genome_sz)]
            elif 'REFERENCE' in config['DATA'][dataset] and 'REGION' in template:
                # we match any parameter containing 'REGION', e.g.
                # MEDAKA_EVAL_REGIONS, MEDAKA_TRAIN_REGIONS etc
                # get reference lengths per region based either on full samtools
                # region str with start and end, or using reference fasta to
                # pull out contig length
                region_params = [k for k in dataset_params if 'REGION' in k]

                # check we don't have more than one region parameter
                if len(region_params) > 1:
                    raise ValueError('Found more than one region paramter in {}'.format(template))

                region_param = region_params[0]
                if not config['DATA'][dataset][region_param]:
                    # the list of regions is empty, skip this dataset
                    continue

                regions = config['DATA'][dataset][region_param]
                ref = config['DATA'][dataset]['REFERENCE']
                dataset_templates = []
                for region in regions:
                    region_sz = int_to_formatted_string(get_region_len(region, ref))
                    d = {'GENOME_SIZE': region_sz, region_param: region}
                    dataset_templates.append(partial_format(dataset_tmp, **d))
        else:
            dataset_templates = [dataset_tmp]

        # kwargs is dict name: list of values which will be expanded with product_dict
        kwargs = {k: config['DATA'][dataset][k] for k in dataset_params if k in config['DATA'][dataset]}
        if len(kwargs) > 0:
            dataset_templates = list(set([partial_format(t, **k) for t in
                                          dataset_templates for k in product_dict(**kwargs)]))
            templates.extend(dataset_templates)
        else:
            templates.extend(dataset_templates)

    # expand dataset-specific templates with non_dataset_params to create targets
    # this will fail if there are any {} still present
    if len(non_dataset_params) > 0:
        targets = [t.format(**k) for k in product_dict(**non_dataset_params) for t in templates]
    else:
        targets = templates

    return targets


def find_genome_size(target, config):
    """Find an appropriate genome size from a target path."""
    target_parts = pathlib.Path(target).parts
    dataset = target_parts[0]
    found = set()
    if 'REFERENCE' in config['DATA'][dataset]:
        # We have a reference file, need to find a contig name in path.
        with pysam.FastaFile(config['DATA'][dataset]['REFERENCE']) as fh:
            rlengths = dict(zip(fh.references, fh.lengths))
        for folder, ref in itertools.product(target_parts, rlengths.keys()):
            if ref in folder:
                found.add(ref)
    if len(found) == 0:
        # fallback to config value if present
        if 'GENOME_SIZE' in config['DATA'][dataset]:
            return config['DATA'][dataset]['GENOME_SIZE']
        else:
            raise KeyError("Could not find a contig name within target: {}.".format(target))
    elif len(found) > 1:
        raise ValueError("Found multiple contig names within target: {}.".format(target))
    return int_to_formatted_string(rlengths[found.pop()])


def suffix_decorate(func, suffix=''):
    """Add a suffix to a function returning a string

    :param: func which returns string
    :suffix: suffix to add to string before returning.
    """
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs) + suffix

    return wrapper
