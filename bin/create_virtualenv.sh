#!/bin/sh
# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║ NOTE: Execute `pwd` and check whether there is a """ character            ║
# ║       (single apostrophe) somewhere in the directory path.                ║
# ║       If yes -> Python will not be able to create a virtual environment.  ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

BLUE="\033[0;34m"
NC="\033[0m"

msg_1="Installing \`virtualenv\` package for Python..."
msg_2="Creating \`virtualenv\`..."
msg_3="Installing dependencies specified in \`requirements.txt\`..."

# Install `virtualenv` Python module (if not already installed).
  echo "\n${BLUE}${msg_1}${NC}"
  pip3 install virtualenv

# Create virtual environment for Python, if not already done.
  echo "\n${BLUE}${msg_2}${NC}"
  if [ -d ../.venv ]; then
      echo "Already exists."
  else
      python3 -m virtualenv ../.venv
  fi

# Install Python requirements.
  echo "\n${BLUE}${msg_3}${NC}"
  ../.venv/bin/pip install -r ../requirements.txt
