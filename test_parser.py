### * Description

# Test parser to see if I can use multiple levels of verbs

### * Set up

### ** Import

import argparse

### * Main

parser = argparse.ArgumentParser(prog = "My Program")
subparsers = parser.add_subparsers(help = "My sub-command help")

parser_a = subparsers.add_parser("Command a", help = "Help for a")
parser_a.set_defaults(action = "Run a")
parser_b = subparsers.add_parser("Command b", help = "Help for b")
parser_b.set_defaults(action = "Run b")

args = parser.parse_args()
print(args)
