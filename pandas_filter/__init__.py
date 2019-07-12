#!/usr/bin/env python
"""init module"""

#from __future__ import absolute_import

import sys
import os
import json
import contextlib
import csv
from builtins import input
from ete3 import NCBITaxa
from dendropy import Tree, DnaCharacterMatrix, DataSet, datamodel

_DEBUG_MK = 1


def debug(msg):
    """short debugging command
    """
    if _DEBUG_MK == 1:
        print(msg)


@contextlib.contextmanager
def cd(path):
    # print 'initially inside {0}'.format(os.getcwd())
    cwd = os.getcwd()
    os.chdir(path)
    # print 'inside {0}'.format(os.getcwd())
    try:
        yield
    except:
        print('Exception caught: ', sys.exc_info()[0])
    finally:
        # print 'finally inside {0}'.format(os.getcwd())
        os.chdir(cwd)


def get_user_input():
    """Asks for yes or no user input.

    :return: user input
    """
    print("get user input")
    is_valid = 0
    x = None
    while not is_valid:
        try:
            x = input("Please write either 'yes' or 'no': ")
            if x == "yes" or x == "no":
                is_valid = 1  # set it to 1 to validate input and to terminate the while..not loop
        except ValueError as e:
            print("'%s' is not a valid answer." % e.args[0].split(": ")[1])
    return x

def write_msg_logfile(msg, workdir, fn="logfile"):
    """
    Write string to logfile.

    :param msg: string
    :return:
    """
    lfd = "{}/{}".format(workdir, fn)
    with open(lfd, "a") as log:
        log.write(msg)
