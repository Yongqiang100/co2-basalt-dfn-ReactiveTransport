#!/bin/bash
# Kept for the commands already in use: the coupled build is build_dirs.sh --variant feedback.
exec bash "$(dirname "$0")/build_dirs.sh" --variant feedback "$@"
