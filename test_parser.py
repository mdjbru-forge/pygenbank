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
subparsers_a = parser_a.add_subparsers(help = "My subsub-command help")
parser_a_a = subparsers_a.add_parser("a_again", help = "Nested sub command")
parser_a_a.set_defaults(action = "a_a")
parser_a_b = subparsers_a.add_parser("b_again", help = "Another nested sub command")
parser_a_b.set_defaults(action = "action")

parser_b = subparsers.add_parser("Command b", help = "Help for b")
parser_b.set_defaults(action = "Run b")

args = parser.parse_args()
print(args)
