import sys

def main():
    output_file = sys.argv[1]
    s1, s2 = sys.argv[2], sys.argv[3]

    with open(output_file, 'w') as f:
        f.write(f"Sequence\tCharge\tMass\tUSI\n")
        f.write(f"1\t2\t3\t{s1}\n")
        f.write(f"4\t5\t6\t{s2}\n")

if __name__ == '__main__':
    main()