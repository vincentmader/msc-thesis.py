#!/usr/bin/env python3
import sys
import re
import os
from termcolor import colored

EXEC = "ipython -i --no-confirm-exit"  # todo: choose ipython or python



def get_example_directories():
    examples = os.listdir('.')
    examples = filter(f, examples)
    examples = list(examples)
    examples = sorted(examples)
    examples = {
        g(e): e for e in examples
    }
    return examples


def get_example_id():
    if len(sys.argv) > 1:
        example_id = get_example_id_from_arg()
    else:
        example_id = get_example_id_from_input(examples)
    return example_id


def get_example_id_from_arg():
    args = sys.argv
    example_id = args[1]
    example_id = int(example_id)  # todo: handle error
    return example_id


def get_example_id_from_input(examples):
    for example in examples:
        msg = examples[example]
        msg = msg.split(" - ")
        msg_0 = colored(msg[0], "blue")
        msg = " - ".join([msg_0, msg[1]])
        print(msg)

    example_id = input(">> ")
    example_id = int(example_id)
    return example_id


def get_example_from_id(example_id):
    example = examples[example_id]  # todo: Handle error.
    return example


def run_example(example):
    cmd = f"cd \"{example}\" && {EXEC} main.py"
    os.system(cmd)


def f(x):
    """Filter directory names: Keep valid examples."""
    if re.search(" - ", x) is None:
        return False
    return True


def g(x):
    """Convert example directory name to integer example-id."""
    x = x.split(" - ")[0]
    x = int(x)
    return x


if __name__ == "__main__":
    examples = get_example_directories()
    example_id = get_example_id()
    example = get_example_from_id(example_id)
    run_example(example)
