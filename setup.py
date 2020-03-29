#!/usr/bin/env python

from setuptools import setup

setup(name='PhylUp',
      version='0.1',
      description='Updating of phylogenetic data',
      author='Martha Kandziora',
      author_email='martha.kandziora@yahoo.com',
      packages=['PhylUp'],
      install_requires=[
          'dendropy',
          'biopython'
      ]
     )
