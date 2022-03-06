import sys
import csv
import itertools

def sliding_window(iterable, size):
    """For a string like 'abcde' and size 3 return ['abc', 'bcd', 'cde']."""
    iters = itertools.tee(iterable, size)
    for i in range(1, size):
        for each in iters[i:]:
            next(each, None)
    return zip(*iters)


def all_frames(start_codons, stop_codons, seq):
    """Return a list of three lists, where each of the inner lists contains the
    possible frame shifts of the `seq`, separated by `stop_codons`.

    Information about which stop codon was encountered is discarded."""
    triplets = list(map(lambda t: "".join(t),
                        sliding_window(seq, 3)))
    frames = [[] for _ in range(3)]
    for i in range(3):
        chunk = list(triplets)[i::3]
        putative_frames = []
        current = []
        for t in chunk:
            if t in stop_codons:
                putative_frames.append(current)
                current = []
            else:
                current.append(t)
        for pf in putative_frames:
            r = list(itertools.dropwhile(lambda t: t not in start_codons, pf))
            if r:
                frames[i].append(r)
    return frames

tr = {"A": "T",
      "T": "A",
      "C": "G",
      "G": "C"}


def reverse_complement(seq):
    return reversed("".join([tr[c] for c in seq]))


def translate_frame(code, frame):
    return "".join(map(lambda x: code[x], frame))


def translate(code, intervals, seqs):
    stop_codons = dict([[triple, aa] for triple, aa in code.items() if aa == "*"])
    start_codons = {"ATG"}
    result = []
    for seq_identifier in seqs.keys():
        if seq_identifier in intervals:
            for start, stop, gene_id in intervals[seq_identifier]:
                seq = seqs[seq_identifier][start:stop]
                rcseq = reverse_complement(seq)
                frames = all_frames(start_codons, stop_codons, seq) + all_frames(start_codons, stop_codons, rcseq)
                translations = [[] for _ in range(len(frames))]
                for i in range(len(frames)):
                    translations[i] = map(lambda r: translate_frame(code, r), frames[i])
                result.append([gene_id, max(itertools.chain.from_iterable(translations), key=len)])

    result.sort(key=lambda x: x[0])
    return result


def read_intervals(file_name):
    gene_intervals = {}
    with open(file_name, "r") as file:
        reader = csv.reader(file, delimiter="\t")
        header = next(reader)
        for seq_identifier, start, stop, gene_id in reader:
            if not seq_identifier in gene_intervals:
                gene_intervals[seq_identifier] = [[int(start), int(stop), gene_id]]
            else:
                gene_intervals[seq_identifier].append([int(start), int(stop), gene_id])

    return gene_intervals

def read_standard_code(file_name):
    code = {}
    with open(file_name, "r") as file:
        reader = csv.reader(file, delimiter="\t")
        next(reader)
        for base_triplet, aa in reader:
            code[base_triplet] = aa
    return code

# Note I don't implement any integrity checks on the fasta file.
# If this was supposed to be a robust implementation, this is not sufficient.
# Further note that I just load all the things into memory
# Obviously, with realistically sized genetic data this might be a bad idea.
def read_sequences(file_name):
    seqs = {}
    seq_identifier = None
    current_seq = None
    with open(file_name, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_seq:
                    seqs[seq_identifier] = current_seq
                seq_identifier = line[1:-1]
                current_seq = ""
            else:
                if line[-1] == "\n":
                    current_seq += line[:-1]
                else:
                    current_seq += line

        if current_seq:
            seqs[seq_identifier] = current_seq

    return seqs

def main():
    args = sys.argv[1:]
    if len(args) != 6 or "-code" not in args or "-intervals" not in args or "-sequences" not in args:
        raise ValueError("Please provide arguments like this: \"translate.py -code path/to/code/file -intervals path/to/intervals/file -sequences path/to/fasta/file")

    print("Reading provided information ...")
    args_dict = dict(zip(args[0::2], args[1::2]))
    code = read_standard_code(args_dict["-code"])
    intervals = read_intervals(args_dict["-intervals"])
    sequences = read_sequences(args_dict["-sequences"])
    print("Translating ...")
    translations = translate(code, intervals, sequences)

    with open("genes.fasta", "w") as file:
        for gene_id, aa_seq in translations:
            file.write(">" + gene_id + "\n" + aa_seq + "\n")

    print("Done.")
    print("Output is printed to file 'genes.fasta'.")

if __name__ == "__main__":
    main()
