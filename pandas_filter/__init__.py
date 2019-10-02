#!/usr/bin/env python
"""init module"""

import sys
import os
import contextlib
from builtins import input
import logging

_DEBUG_MK = 0


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
    :param workdir: working directory as string
    :param fn: filename for logging
    :return:
    """
    lfd = "{}/{}".format(workdir, fn)
    with open(lfd, "a") as log:
        log.write(msg)


def license_print():
    """
    Print the corresponding license
    :return:
    """
    # todo needs to be added to start of program when finished
    sys.stdout.write(
    """
    PhylUp: automatic updating of phylogenetic data
    Copyright (C) 2019  M. Kandziora

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    """)


# global logger
# http://docs.python.org/library/logging.html
# LOG = logging.getLogger("")
# logging.basicConfig(filename='debug.log',level=logging.DEBUG,
#                    format='%(levelname)s [%(asctime)s]: %(message)s')


def setup_logger(name, log_file, level=logging.INFO, writemode="w"):
    """setups as many loggers as you want.
    """
    formatter = logging.basicConfig(filemode="w+", format='%(levelname)s [%(asctime)s]: %(message)s')

    log_fn = "{}".format(log_file)
    handler = logging.FileHandler(log_fn)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

# setup different loggers
log_dict = {}
log_dict2 = {}


def log_info(text, wd):
    """
    setup info logger"""
    if wd not in log_dict:
        dlog = setup_logger("info_log", '{}/info.log'.format(wd), level=logging.INFO, writemode="w+")
        dlog.propagate = True  # logs to file and console
        log_dict[wd] = dlog
    log_dict[wd].info(text)


def log_debug(text, wd):
    """
    setup debug logger"""
    if wd not in log_dict2:
        dlog = setup_logger("debug_log", '{}/debug.log'.format(wd), level=logging.DEBUG, writemode="w+")
        dlog.propagate = False  # logs to file only
        log_dict2[wd] = dlog
    log_dict2[wd].debug(text)
