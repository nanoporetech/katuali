from glob import glob
import os
import re
from setuptools import setup, find_packages
from setuptools import Distribution
from setuptools.command.install import install
import shutil
import sys

__pkg_name__ = 'katuali'
__author__ = 'mwykes'
__description__ = 'Pipelines for Nanopore sequencing.'
__path__ = os.path.abspath(os.path.dirname(__file__))
__pkg_path__ = os.path.join(os.path.join(__path__, __pkg_name__))


# Use readme as long description and say its github-flavour markdown
kwargs = {'encoding':'utf-8'} if sys.version_info.major == 3 else {}
with open(os.path.join(__path__, 'README.md'), **kwargs) as f:
    __long_description__ = f.read()
__long_description_content_type__ = 'text/markdown'


# Get the version number from __init__.py, and exe_path
version_file = os.path.join(__path__, __pkg_name__, '__init__.py')
verstrline = open(version_file, 'r').read()
vsre = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(vsre, verstrline, re.M)
if mo:
    __version__ = mo.group(1)
else:
    raise RuntimeError('Unable to find version string in "{}".'.format(version_file))


# Scrape requirements from requirements.txt
install_requires = []
with open(os.path.join(__path__, 'requirements.txt')) as fh:
    reqs = (
        r.split('#')[0].strip()
        for r in fh.read().splitlines() if not r.startswith('#')
    )
    for req in reqs:
        if req == '':
            continue
        if req.startswith('git+https'):
            req = req.split('/')[-1].split('@')[0]
        install_requires.append(req)


setup(
    name=__pkg_name__,
    version=__version__,
    author=__author__,
    author_email='{}@nanoporetech.com'.format(__author__),
    description=__description__,
    long_description=__long_description__,
    long_description_content_type=__long_description_content_type__,
    dependency_links=[],
    install_requires=install_requires,
    tests_require=[].extend(install_requires),
    python_requires='>=3.5.2, <3.7',
    packages=find_packages(exclude=['*.test', '*.test.*', 'test.*', 'test']),
    package_data={
        __pkg_name__:['data/*', 'data/test/*']
    },
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'process_katuali_config = {}:process_katuali_config'.format(__pkg_name__),
            'pick_gpu = {}:pick_gpu'.format(__pkg_name__),
            'katuali_datafile = {}:print_data_path'.format(__pkg_name__),
            'katuali_config = {}:create_config'.format(__pkg_name__),
        ]
    },
    scripts=[
        os.path.join('scripts', x) for x in ('katuali',)
    ]
)


