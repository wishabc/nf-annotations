import sys


def main(fasta_file, window, out_file):
    ref, alt = None, None
    with open(fasta_file) as f, open(out_file, 'w') as out:
        for line in f:
            if line.startswith('>'):
                _, _, _, ref, alt, _ = line.split('@')
                out.write(line.replace('@ref', '@alt'))
            else:
                assert ref == line[window]
                out.write(line[:window] + alt + line[window + 1:])

if __name__ == '__main__':
    fasta = sys.argv[1]
    window = int(sys.argv[2])
    out_file = sys.argv[3]
    main(fasta, window, out_file)        