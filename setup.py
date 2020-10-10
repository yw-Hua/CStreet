from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')


setup(
    name='cstreet',  # Required

    version='0.0.3',  # Required

    description='CStreet is a python script (python 3.6 or higher) for cell states trajectory construction by using k-nearest neighbors graph algorithm for time-series single-cell RNA-seq data.',  # Optional

    long_description=long_description,  # Optional

    long_description_content_type='text/markdown',  # Optional (see note above)

    url='https://github.com/yw-Hua/CStreet',  # Optional

    author='Yuwei Hua',  # Optional

    author_email='ywhua@tongji.edu.com',  # Optional

    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish
        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 3.6',
    ],

    package_dir={'': 'src'},
    packages=[''],

    python_requires='>=3.6, <4',

    install_requires=[
        'pandas==0.25.3',
        'scanpy==1.5.0',
        'anndata==0.7.2',
        'networkx==2.4',
        'fa2',
        'retrying'

    ],

)