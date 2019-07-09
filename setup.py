#!/usr/bin/env python

from setuptools import setup

setup(name='pandas_filter',
      version='0.1',
      description='pandas_filter',
      author='Martha Kandziora',
      author_email='martha.kandziora@yahoo.com',
      packages=['pandas_filter'],
      install_requires=[
          'dendropy',
#          'configparser',
          'biopython',
#          'urllib3',
          'peyotl',
      ]
     )
