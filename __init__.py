from decimal import Decimal
import tracemalloc
from time import perf_counter
import pyperclip
import requests
import json
from random import randint
from typing import Tuple, Dict, List


class ProteinMotif():
    def __init__(self, motif_string: str) -> None:
        self.motif_string = motif_string
        self.positions = []

        temp_list: List[str] = []
        bracket_flag = False
        brace_flag = False

        for c in self.motif_string:
            if c == '{':
                brace_flag = True
                continue
            elif c == '[':
                bracket_flag = True
                continue

            if c == '}':
                brace_flag = False
                second_temp_list = []
                for p in protein_alphabet:
                    if p not in temp_list:
                        second_temp_list.append(p)
                self.positions.append(second_temp_list)
                second_temp_list = []
                temp_list = []
                continue
            elif c == ']':
                bracket_flag = False
                self.positions.append(temp_list)
                temp_list = []
                continue

            if bracket_flag:
                temp_list.append(c)

            elif brace_flag:
                temp_list.append(c)

            else:
                self.positions.append([c])

    def __len__(self) -> int:
        return len(self.positions)

    def match(self, substring: str) -> bool:
        for i in range(len(substring)):
            if substring[i] not in self.positions[i]:
                return False
        return True


class Protein():
    def __init__(self, protein: str, id: str = '') -> None:
        self.protein = protein
        self.id = id

    def __len__(self) -> int:
        return len(self.protein)

    def __str__(self) -> str:
        if self.id != '':
            return f'{self.id}: {self.protein}'
        else:
            return f'{self.protein}'

    def find_motifs(self, motif):
        motif_positions = []
        for i in range(len(self.protein) - len(motif) + 1):
            substring_to_check = self.protein[i:i+len(motif)]
            if motif.check(substring_to_check):
                motif_positions.append(i + 1)
        return motif_positions


class Sequence():
    def __init__(self, dna: str, id: str = '',
                 *, is_rna: bool = False,
                 update_immediately=False) -> None:
        if not is_rna:
            self.dna = dna.upper()
            self.rna = self.dna.replace('T', 'U').upper()
        else:
            self.rna = dna
            self.dna = self.rna.replace('U', 'T').upper()

        self.id = id

        if update_immediately:
            self.__update_after_changes__()

    def __len__(self) -> int:
        return len(self.dna)

    def __str__(self) -> str:
        if self.id != '':
            return f'{self.id}: {self.dna}'
        else:
            return f'{self.dna}'

    def __eq__(self, seq2: object) -> bool:
        if not isinstance(seq2, Sequence):
            return NotImplemented
        return self.dna == seq2.dna

    def __form_peptide_chain__(self):
        self.peptide_chain = ''
        rna_copy = [x for x in self.rna]
        rna_copy.reverse()

        for _ in range(len(rna_copy) // 3):
            codon = ''
            for _ in range(3):
                codon += rna_copy.pop()
            self.peptide_chain += codon_table[codon]

    def __find_skew__(self):
        self.skew = [0]
        skew = 0
        for base in self.dna:
            if base == 'C':
                skew -= 1
            elif base == 'G':
                skew += 1
            self.skew.append(skew)

    def __form_reverse_complement_peptide_chain__(self):
        self.reverse_complement_peptide_chain = ''
        rna_copy = [x for x in self.rna_complement]

        for _ in range(len(rna_copy) // 3):
            codon = ''
            for _ in range(3):
                codon += rna_copy.pop()
            self.reverse_complement_peptide_chain += codon_table[codon]

    def __update_after_changes__(self):
        self.dna_base_count = {'G': self.dna.count('G'),
                               'C': self.dna.count('C'),
                               'T': self.dna.count('T'),
                               'A': self.dna.count('A')}
        self.rna_base_count = {'G': self.rna.count('G'),
                               'C': self.rna.count('C'),
                               'U': self.rna.count('U'),
                               'A': self.rna.count('A')}

        try:
            self.dna_complement = ''.join(
                [dna_complements[x] for x in self.dna])
            self.dna_reverse_complement = (self.dna_complement)[::-1]
        except KeyError as e:
            print(f'DNA sequences cannot contain {e}')

        try:
            self.rna_complement = ''.join(
                [rna_complements[x] for x in self.rna])
            self.rna_reverse_complement = (self.rna_complement)[::-1]
        except KeyError as e:
            print(f'RNA sequences cannot contain {e}')

        self.__form_peptide_chain__()
        self.__form_reverse_complement_peptide_chain__()
        self.__find_skew__()

    def gc_content(self, round_places: int = 3) -> Decimal:
        numerator = Decimal(self.dna_base_count['G']
                            + self.dna_base_count['C'])
        denominator = Decimal(len(self.dna))
        return round(numerator / denominator, round_places)

    def find_substrings(self, shortest: int, longest: int = 0,
                        *, window: Tuple[int, int] = (0, 0)) -> Dict[str, int]:
        """Returns a dict containing all the
           unique substrings and their counts"""
        if window == (0, 0):
            window = (0, len(self.dna))
        if longest == 0:
            longest = shortest

        substring_chains: Dict[str, int] = dict()
        for length in range(shortest, longest+1):
            for i in range(window[0], window[1] - length + 1):
                try:
                    substring_chains[self.dna[i:i+length]] += 1
                except KeyError:
                    substring_chains[self.dna[i:i+length]] = 1
        substring_chains = dict(sorted(
            substring_chains.items(), key=lambda x: x[1]))
        return substring_chains

    def find_clumps(self, k: int, clump_size: int, window_size: int) -> set[str]:
        substring_chains: List[str] = []
        for i in range(len(self.dna) - window_size + 1):
            print(f'Searching {i+1}/{len(self.dna) - window_size + 1}')
            substrings = self.find_substrings(k, window = (i, i + window_size))
            for ke, v in substrings.items():
                if v >= clump_size:
                    substring_chains.append(ke)
        substring_chains = set(substring_chains)
        return substring_chains
                    

    def splice(self, intron: str) -> None:
        while intron in self.rna:
            self.rna = self.rna.replace(intron, '')
        self.__update_after_changes__()

    def find_motif(self, dna_chain: str) -> List[int]:
        locations = []
        m_len = len(dna_chain)
        for i in range(0, len(self.dna) - m_len):
            if dna_chain == self.dna[i:i+m_len]:
                locations.append(i)

        return locations

# alanine          = ala = A
# arginine         = arg = R
# asparagine       = asn = N
# aspartic acid    = asp = D
# asparagine       = asx = B **Not used
# cysteine         = cys = C
# glutamic acid    = glu = E
# glutamine        = gln = Q
# glutamine        = glx = Z
# glycine          = gly = G
# histidine        = his = H
# isoleucine       = ile = I
# leucine          = leu = L
# lysine           = lys = K
# methionine       = met = M
# phenylalanine    = phe = F
# proline          = pro = P
# serine           = ser = S
# threonine        = thr = T
# tryptophan       = trp = W
# tyrosine         = tyr = Y
# valine           = val = V


rna_bases = ['U', 'C', 'A', 'G']
dna_bases = ['T', 'C', 'A', 'G']
codons = [a + b + c for a in rna_bases for b in rna_bases for c in rna_bases]
amino_acids_string = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAA \
                      ADDEEGGGG'
codon_table = dict(zip(codons, amino_acids_string))

dna_complements = {'T': 'A', 'A': 'T', 'C': 'G', 'G': 'C'}
rna_complements = {'U': 'A', 'A': 'U', 'C': 'G', 'G': 'C'}


# List of monoisotopic masses of elements
element_mono_mass = {
    'h': Decimal(1.007825035),
    'c': Decimal(12.000000000),
    'n': Decimal(14.003074000),
    'o': Decimal(15.994914630),
    # 'p': Decimal(30.973762000),
    's': Decimal(31.972070700),
    # 'se': Decimal(79.916520)
}

# List of average masses of elements
element_average_mass = {
    'h': Decimal(1.007940000),
    'c': Decimal(12.010700000),
    'n': Decimal(14.006700000),
    'o': Decimal(15.999400000),
    # 'p': Decimal(30.973761000),
    's': Decimal(32.065000000),
    # 'se': Decimal(78.960000000)
}

# Elemental composition of amino acids
amino_acid_composition = {
    'A': {'c': 3, 'h': 5, 'n': 1, 'o': 1, 's': 0},
    'R': {'c': 6, 'h': 12, 'n': 4, 'o': 1, 's': 0},
    'N': {'c': 4, 'h': 6, 'n': 2, 'o': 2, 's': 0},
    'D': {'c': 4, 'h': 5, 'n': 1, 'o': 3, 's': 0},
    'C': {'c': 3, 'h': 5, 'n': 1, 'o': 1, 's': 1},
    'E': {'c': 5, 'h': 7, 'n': 1, 'o': 3, 's': 0},
    'Q': {'c': 5, 'h': 8, 'n': 2, 'o': 2, 's': 0},
    'G': {'c': 2, 'h': 3, 'n': 1, 'o': 1, 's': 0},
    'H': {'c': 6, 'h': 7, 'n': 3, 'o': 1, 's': 0},
    'I': {'c': 6, 'h': 11, 'n': 1, 'o': 1, 's': 0},
    'L': {'c': 6, 'h': 11, 'n': 1, 'o': 1, 's': 0},
    'K': {'c': 6, 'h': 12, 'n': 2, 'o': 1, 's': 0},
    'M': {'c': 5, 'h': 9, 'n': 1, 'o': 1, 's': 1},
    'F': {'c': 9, 'h': 9, 'n': 1, 'o': 1, 's': 0},
    'P': {'c': 5, 'h': 7, 'n': 1, 'o': 1, 's': 0},
    'S': {'c': 3, 'h': 5, 'n': 1, 'o': 2, 's': 0},
    'T': {'c': 4, 'h': 7, 'n': 1, 'o': 2, 's': 0},
    'W': {'c': 11, 'h': 10, 'n': 2, 'o': 1, 's': 0},
    'Y': {'c': 9, 'h': 9, 'n': 1, 'o': 2, 's': 0},
    'V': {'c': 5, 'h': 9, 'n': 1, 'o': 1, 's': 0}
}

# Masses of different amino masses, both average and monoisotopic
amino_acid_mono_mass = dict()
amino_acid_average_mass = dict()

for amino_acid in amino_acid_composition.keys():
    mono_total = 0.0
    average_total = 0.0
    for mono_mass in element_mono_mass:
        mono_total += float((amino_acid_composition[amino_acid][mono_mass]
                             * element_mono_mass[mono_mass]))
    for average_mass in element_average_mass:
        average_total += float((amino_acid_composition
                                [amino_acid][average_mass]
                                * element_average_mass[average_mass]))

    amino_acid_mono_mass[amino_acid] = mono_total
    amino_acid_average_mass[amino_acid] = average_total

protein_alphabet = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L',
                    'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


def concat_sequences(sequences: List[Sequence]) -> Sequence:
    sequence_string = ''
    for i, sequence in enumerate(sequences):
        sequence_string += sequence.dna

    return Sequence(sequence_string)


def generate_random_sequence(gc: int, seq_length: int) -> Sequence:
    temp_seq_string = ''
    for _ in range(seq_length):

        r = randint(1, 100)
        if r <= gc:
            r2 = randint(1, 100)
            if r2 <= 50:
                temp_seq_string += 'G'
            else:
                temp_seq_string += 'C'

        else:
            r2 = randint(1, 100)
            if r2 <= 50:
                temp_seq_string += 'A'
            else:
                temp_seq_string += 'T'

    return Sequence(temp_seq_string)


def generate_random_sequences(num_sequences: int, gc: Tuple[float, float],
                              seq_length: Tuple[int, int],
                              *, equilength: bool = True) -> list[Sequence]:
    sequences = []
    if equilength:
        sequence_length = randint(seq_length[0], seq_length[1])

    for _ in range(num_sequences):
        gc_content = randint(int(gc[0] * 100), int(gc[1] * 100))

        if not equilength:
            sequence_length = randint(seq_length[0], seq_length[1])

        sequences.append(generate_random_sequence(gc_content, sequence_length))

    return sequences


def read(filename: str, prefix: str = 'data/') -> list[str]:
    with open(f'{prefix}{filename}') as f:
        input_data = f.readlines()
        input_data = [x.strip() for x in input_data]
        return input_data


def read_fasta_data(filename: str, prefix: str = 'data/') -> list[Sequence]:
    fasta_data = dict()
    data_array = read(filename, prefix)
    current_id = data_array.pop(0)[1:]
    current_fasta = ''

    for d in data_array:
        if d[0] != '>':
            current_fasta += d
        else:
            fasta_data[current_id] = current_fasta
            current_fasta = ''
            current_id = d[1:]
    fasta_data[current_id] = current_fasta
    sequences = []
    for k, v in fasta_data.items():
        sequences.append(Sequence(v, k))

    return sequences


def read_uniprot(access_ids: list[str]) -> list[Protein]:
    base_url = 'https://rest.uniprot.org/uniprotkb/'
    proteins: list[Protein] = []
    for id in access_ids:
        real_id = id.split('_', 1)[0]
        url = base_url + real_id
        print(f'Fetching {url}')
        downloaded_data = (requests.get(url, allow_redirects=True).text)
        downloaded_data = json.loads(downloaded_data)
        downloaded_data = downloaded_data['sequence']['value']
        proteins.append(Protein(downloaded_data, id))
    return proteins


def hamming_distance(s: Sequence, t: Sequence) -> int:
    hamming = 0
    for base_s, base_t in zip(s.dna, t.dna):
        if base_s != base_t:
            hamming += 1
    return hamming


def find_consensus(consensus_matrix: dict[str, list[int]]) -> Sequence:
    bases = ''
    for i in range(len(consensus_matrix['A'])):
        max_index = 0
        max_base = ''
        for key in consensus_matrix.keys():
            if consensus_matrix[key][i] > max_index:
                max_index = consensus_matrix[key][i]
                max_base = key
        bases += max_base
    return Sequence(bases)


def find_consensus_matrix(sequences: list[Sequence]) -> Dict[str, List[int]]:
    try:
        check_sequence_lengths(sequences)
    except SequenceLengthException as e:
        print(e)
        exit()

    A = [0] * len(sequences[0])
    C = [0] * len(sequences[0])
    G = [0] * len(sequences[0])
    T = [0] * len(sequences[0])
    d = {'A': A, 'C': C, 'G': G, 'T': T}

    for sequence in sequences:
        for i, base in enumerate(sequence.dna):
            d[base][i] += 1

    return d


def check_sequence_lengths(sequences: list[Sequence]) -> None:
    normal_length = len(sequences[0])
    for sequence in sequences:
        if len(sequence) != normal_length:
            raise SequenceLengthException


class SequenceLengthException(Exception):
    def __init__(self,
                 message="The sequences entered are not the same length"):
        self.message = message
        super().__init__(self.message)


n_glycosylation = ProteinMotif("N{P}[ST]{P}")


def time_and_memory_decorator(func):

    def wrapper(*args, **kwargs):
        tracemalloc.start()
        start_time = perf_counter()
        answer = func(*args, **kwargs)
        end_time = perf_counter()
        current, peak = tracemalloc.get_traced_memory()

        peak = f'Peak memory use: {peak / 10**6:0.6f} MB'
        tim = f'Time to execute in seconds: {end_time - start_time:0.6f}'

        tracemalloc.stop()

        try:
            pyperclip.copy(answer)
        except pyperclip.PyperclipException:
            pyperclip.copy("Please manually copy and paste your answer")

        return answer, peak, tim

    return wrapper

def auto_timer(func):

    def wrapper(*args, **kwargs):
        start_time = perf_counter()
        func(*args, **kwargs)
        end_time = perf_counter()
        print(f'{end_time - start_time:0.6f}')
        
