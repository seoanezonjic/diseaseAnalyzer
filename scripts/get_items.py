#!/usr/bin/env python

import json
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-i', '--input_file', dest='input', help='Path to input file')
args = parser.parse_args()

with open(args.input, 'r') as f:
    data = json.load(f)

items = []
for assoc in data['associations']:
    term = assoc['object']['id']
    if term not in items: items.append(term)

for item in items:
    print(item)
