# write a command line tool that takes two parameters as inputs for generation a summary html file

import sys
import os
import re
import argparse

def main():
    output_file = sys.argv[1]
    s1, s2 = sys.argv[2], sys.argv[3]

    with open(output_file, 'w') as f:
        f.write('<html><body><h1>Summary</h1><p>Param 1:{}\nParam 2: {}\n\nHello world!</p></body></html>'.format(s1, s2))

if __name__ == '__main__':
    main()