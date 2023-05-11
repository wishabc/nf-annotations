import sys


def main(fasta_file, window, out_file):
    ref, alt = None, None
    with open(fasta_file) as f, open(out_file) as out:
        for line in f:
            if line.startswith('>'):
                _, _, ref, alt, _ = line.split('@')
                out.write(line.replace('@ref', '@alt'))
            else:
                assert ref == line[window]
                out.write(line[:window + 1] + alt + line[window + 1:], line)

if __name__ == '__main__':
    main(*sys.argv[1:])        