
# ----------- IMPORTS ---------------------------- ||
from setuptools import setup, find_packages # Always prefer setuptools over distutils
from pathlib import Path

# ----------- GLOBALS ---------------------------- ||
# current location
here = Path(__file__).parent.resolve()
# version
__version__ = ''
with open(here.joinpath('VERSION'), mode='r') as ifd:
	__version__ = ifd.readline().rstrip('\n')
# author & email
__author__ = "Brandon Pickett"
__author_email__ = "pickettbd@byu.edu"
# long description
long_description = ''
with open(here.joinpath('README.md'), encoding='utf-8', mode='r') as f:
    long_description = f.read()


# ----------- CLASSES ---------------------------- ||
# None

# ---------- FUNCTIONS --------------------------- ||
# None

# ------------- MAIN ----------------------------- ||
# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    # This is the name of your project. The first time you publish this
    # package, this name will be registered for you. It will determine how
    # users can install this project, e.g.:
    #
    # $ pip install sampleproject
    #
    # And where it will live on PyPI: https://pypi.org/project/sampleproject/
    #
    # There are some restrictions on what makes a valid project name
    # specification here:
    # https://packaging.python.org/specifications/core-metadata/#name
    name='taxaResiliency',  # Required

    # Versions should comply with PEP 440:
    # https://www.python.org/dev/peps/pep-0440/
    #
    # For a discussion on single-sourcing the version across setup.py and the
    # project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=__version__,   # Required

    # This is a one-line description or tagline of what your project does. This
    # corresponds to the "Summary" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#summary
    description="Calculate the Taxa Resiliency Score for trees",  # Optional

    # This is an optional longer description of your project that represents
    # the body of text which users will see when they visit PyPI.
    #
    # Often, this is the same as your README, so you can just read it in from
    # that file directly (as we have already done above)
    #
    # This field corresponds to the "Description" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#description-optional
    long_description=long_description,  # Optional

    # Denotes that our long_description is in Markdown; valid values are
    # text/plain, text/x-rst, and text/markdown
    #
    # Optional if long_description is written in reStructuredText (rst) but
    # required for plain-text or Markdown; if unspecified, "applications should
    # attempt to render [the long_description] as text/x-rst; charset=UTF-8 and
    # fall back to text/plain if it is not valid rst" (see link below)
    #
    # This field corresponds to the "Description-Content-Type" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#description-content-type-optional
    long_description_content_type='text/markdown',  # Optional (see note above)

    # This should be a valid link to your project's main homepage.
    #
    # This field corresponds to the "Home-Page" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#home-page-optional
    url='https://github.com/pickettbd/taxonResiliency',  # Optional

    # This should be your name or the name of the organization which owns the
    # project.
    author=__author__,  # Optional

    # This should be a valid email address corresponding to the author listed
    # above.
    author_email=__author_email__,  # Optional

    # Classifiers help users find your project by categorizing it.
    #
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[  # Optional
        # How mature is this project?
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
		"Intended Audience :: Science/Research",
		"Topic :: Scientific/Engineering",
		"Topic :: Scientific/Engineering :: Bio-Informatics",

        # Pick your license as you wish
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate you support Python 3. These classifiers are *not*
        # checked by 'pip install'. See instead 'python_requires' below.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3 :: Only',

		# other classifiers
		"Environment :: Console",
		"Natural Language :: English",
		"Operating System :: MacOS",
		"Operating System :: POSIX :: Linux",
    ],

    # This field adds keywords for your project which will appear on the
    # project page. What does your project relate to?
    #
    # Note that this is a list of additional keywords, separated
    # by commas, to be used to assist searching for the distribution in a
    # larger catalog.
    keywords='phylogeny, phylogenetic, taxonomy, taxa, taxon, resilience resiliency, tree, systematics',  # Optional

    # When your source code is in a subdirectory under the project root, e.g.
    # `src/`, it is necessary to specify the `package_dir` argument.
    package_dir={'': 'src'},  # Optional

	# include package data
	include_package_data=True,

	# find packages/modules/etc.
    packages=find_packages(where='src'),  # Required

    # Specify which Python versions you support. In contrast to the
    # 'Programming Language' classifiers above, 'pip install' will check this
    # and refuse to install the project if the version does not match. See
    # https://packaging.python.org/guides/distributing-packages-using-setuptools/#python-requires
    python_requires='>=3.6, <4',

    # This field lists other packages that your project depends on to run.
    # Any package you put here will be installed by pip when your project is
    # installed, so they must be valid existing projects.
    #
    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=["sys", "re", "pkgutil", "argparse", "pathlib"],  # Optional

	# setup_requires
	setup_requires=["pathlib"],  # Optional

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword.
    entry_points={  # Optional
        'console_scripts': [
            'taxaResil=taxaResiliency.calcScore:main'
        ]
    },

    # List additional URLs that are relevant to your project as a dict.
    project_urls={  # Optional
        'Source': 'https://github.com/pickettbd/taxonResiliency',
        'Bug Reports': 'https://github.com/pickettbd/taxonResiliency/issues',
        #'Funding': 'https://donate.pypi.org',
        'Buy me a soda!': 'http://buymeacoff.ee/pickettbd'
    }
)
