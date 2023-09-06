#!/bin/sh
# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║ NOTE: Execute `pwd` and check whether there is a "'" character            ║
# ║       (single apostrophe) somewhere in the directory path.                ║
# ║       If yes -> Python will not be able to create a virtual environment.  ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# Install `virtualenv` Python module (if not already installed).
  pip3 install virtualenv

# Create virtual environment for Python, if not already done.
  [ -d ../.venv ] || python3 -m virtualenv ../.venv

# Install Python requirements.
  ../.venv/bin/pip install -r ../requirements.txt
