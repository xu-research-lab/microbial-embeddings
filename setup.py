#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()


requirements = ['click>=8.0', 'tqdm>=4.64.0', 'numba',
                'numpy>=1.22.4', 'scipy>=1.8.1', 'matplotlib',
                'biom-format>=2.1.10', 'scikit-bio>=0.5.6',
                'scikit-learn>=0.24.1', 'pandas>=1.4', 'joblib>=1.2.0',
                'torch>=1.10.0']

extras_requirements = {'progress_bar': ['tqdm'],
                       'profiler': ['line-profiler'],
                       'doc': ["Sphinx >= 1.4", "sphinx-autodoc-typehints", "nbsphinx"]}

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3', ]

setup(
    author="microbial language model team",
    author_email='zhuzhengnong@xbiome.com',
    python_requires='>=3.9',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    description="microbial embedding model",
    entry_points={
        'console_scripts': [
            'membed=membed.cli:main',
        ],
    },
    install_requires=requirements,
    extras_require=extras_requirements,
    license="MIT license",
    long_description=readme, # + '\n\n' + history,
    include_package_data=True,
    keywords=['microbiome', 'language model', 'NLP'],
    name='membed',
    packages=find_packages(include=['membed', 'membed.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='',
    version='1.0.0',
    zip_safe=False,
)
