#!/usr/bin/env bash

# SPDX-License-Identifier: MIT
# Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#
#      Author : Luc Grosheintz <forbugrep@zoho.com>
#     Created : 2015-01-02
#
# Description : Push sources to a remote computer
read -r -d '' help_text << EOF
usage: $(basename $0) (push | pull) [--remote <host>]
    Transefer data to or from a remote remote.

    pull    pull data in 'data', optionally specify a stems to match.
    push    copy source to the remote site

    --remote <remote> specify the remote location to copy to/from
EOF

if [[ $# < 1 ]]
then
  echo "$help_text"
  exit -1
fi

if [[ $# -ge 1 ]]
then
  case $1 in
    push|pull)
      mode="$1"
      ;;
    --help)
      echo "$help_text"
      exit 0
      ;;
    *)
      echo "You need to specify a mode, either 'push' or 'pull'."
      exit -1
      ;;
  esac

  shift
fi

# read-in the other for push
while [[ $# -ge 1 ]]
do
  case $1 in
    --remote)
      shift
      remote="$1"
      ;;
    *)
      if [[ $mode == "pull" ]]
      then
        stems="$@"
        break
      else
        echo "$help_text"
        exit -1
      fi
      ;;
  esac

  shift
done

## Configuration section.

# want to only test? (uncomment)
# DRY_RUN=--dry-run

# Set destination for each location
case $remote in
  daint)
    dest='lucg@ela.cscs.ch:~'
    ;;
  ada-*)
    dest="lucg@${remote}:/userdata/lucg/"
    ;;
  euler)
    dest="lucg@euler.ethz.ch:/cluster/scratch/lucg"
    ;;
  leonhard)
    dest="lucg@login.leonhard.ethz.ch:/cluster/scratch/lucg"
    ;;
  *)
    echo "Unknown remote."
    exit -1
    ;;
esac

## Actual copying up and down
case $mode in
  push)
    rsync -urv --delete $DRY_RUN \
          --include-from rsync.include --exclude-from rsync.exclude \
          $(realpath $PWD) $dest
    ;;
  pull)
    rsync -ruv $DRY_RUN "$dest/$(basename $(realpath $PWD))/data/*" data
    ;;
  *)
    echo "Unknown mode."
    exit -1
    ;;
esac
