import os
import re
from setuptools import setup

# Allow setup.py to be run from any path.
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))
ROOT = os.path.abspath(os.path.dirname(__file__))

project_name = 'vcf2prs'
# Read the version without importing runtime dependencies during an isolated build.
with open(os.path.join(ROOT, project_name, '__init__.py'), encoding='utf-8') as init_file:
    project_version = re.search(
        r"^__version__\s*=\s*['\"]([^'\"]+)['\"]", init_file.read(), re.MULTILINE
    ).group(1)
project_readme_fname = 'README.rst'
project_author = 'Lorenzo Ficorella; Andrew Lee'


def readRequirements(fname):
    '''
    Read the requirements from the requirements file
    '''
    requirements = []
    if os.path.exists(fname):
        with open(fname) as fp:
            requirements = fp.read().splitlines()
    return [r.replace('==', '>=') for r in requirements]


setup(
    name=project_name,
    version=project_version,
    packages=[project_name],
    entry_points={
        'console_scripts': ['vcf2prs=vcf2prs.cli:main'],
    },
    package_data={},
    author=project_author,
    include_package_data=False,
    zip_safe=False,
    url='http://github.com/CCGE-BOADICEA/SHARE-PRScalculations',
    description=("Calculate a person's PRS (polygenic risk score) from their "
                 "VCF (variant call format) file and/or their ancestry proportions."),
    long_description=open(os.path.join(ROOT, project_readme_fname)).read(),
    install_requires=readRequirements('requirements.txt'),
    classifiers=[
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
    ],
)
