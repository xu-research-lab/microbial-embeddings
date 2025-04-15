#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.org') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['click>=8.0', 'tqdm', 'numba',
                'numpy', 'scipy', 'matplotlib',
                'biom-format', 'scikit-bio',
                'scikit-learn', 'torch']

extras_requirements = {'progress_bar': ['tqdm'],
                       'profiler': ['line-profiler'],
                       'doc': ["Sphinx >= 1.4", "sphinx-autodoc-typehints", "nbsphinx"]}

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3', ]

setup(
    author="microbial language model team",
    author_email='zhuzhengnong@xbiome.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="microbial embedding model",
    entry_points={
        'console_scripts': [
            'membed=membed.cli:main',
        ],
    },
    install_requires=requirements,
    extras_requires=extras_requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords=['microbiome', 'language model', 'NLP'],
    name='membed',
    packages=find_packages(include=['membed', 'membed.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='',
    version='0.1.0',
    zip_safe=False,
)
