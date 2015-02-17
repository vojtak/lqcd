#!/bin/bash

CONF_SET_FOLDER=./configurations2
DO_ARG_FOLDER=$CONF_SET_FOLDER

rm in.conf_list

# create list of configurations
find $CONF_SET_FOLDER/ ! -name "*do_arg*" ! -name "*~" -print | tail -n +2 | sort > in.conf_list


