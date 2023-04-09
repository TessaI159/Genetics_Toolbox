from decimal import Decimal
import tracemalloc
from time import perf_counter
import pyperclip
import requests
import json
from random import randint
from typing import Tuple, Dict, List, Union


class ProteinMotif():
    """Class representing a motif do be found in a protein.
    Letters indicate proteins,
    {a, b, c...} indicates any protein except a, b, c...,
    [a, b, c...] indicates any protein a, b, c...

    :param motif_string: A string as described above
    :type motif_string: str
    """
    # TODO (Tess): This constructor might be too long.
    # Try to shorten it, or break it up
    def __init__(self, motif_string: str) -> None:
        """Constructor"""
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
        """Returns the length of the motif. Any protein being compared
        to this motif should be at least this long.

        :return: Length of motif
        :rtype: int
        """
        return len(self.positions)

    def match(self, substring: str) -> bool:
        """Checks if a given string matches the motif pattern

        :param substring: String to check against motif pattern
        :type substring: str
        :return: True if matches, False otherwise
        :rtype: bool
        """
        if len(substring) != len(self.positions):
            return False
        for i in range(len(substring)):
            if substring[i] not in self.positions[i]:
                return False
        return True


class Protein():
    """Class that represents a protein as a string of amino acids

    :param protein: Amino acid string
    :type protein: str
    :param id: Name of protein
    :type id: str
    """
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

    def find_motifs(self, motif: Union[ProteinMotif, str]) -> List[int]:
        """Searches all substrings that match len(motif) to see if they match.

        :param motif: A protein motif you wish to find
        inside the current protein
        :type motif: ProteinMotif
        :return: All starting positions of protein that match motif
        :rtype: List[int]
        """
        if isinstance(motif, str):
            motif = ProteinMotif(motif)
        motif_positions = []
        for i in range(len(self.protein) - len(motif) + 1):
            substring_to_check = self.protein[i:i+len(motif)]
            if motif.match(substring_to_check):
                motif_positions.append(i + 1)
        return motif_positions


class Sequence():
    """Class representing a strand of DNA or RNA

    :param dna: String of A, C, G, T
    :type dna: str
    :param id: Name of DNA strand
    :type id: str
    :param is_rna: Tells if strand is RNA or DNA;
    probably redundant, we could just check for the presence of U
    :type is_rna: bool
    :param update_immediately: Tells if the DNA should
    immediately run __update_after_changes__
    :type update_immediately: bool"""
    def __init__(self, dna: str, id: str = '',
                 *, is_rna: bool = False,
                 update: bool = False) -> None:
        """Constructor"""
        if not is_rna:
            self.dna = dna.upper()
            self.rna = self.dna.replace('T', 'U').upper()
        else:
            self.rna = dna
            self.dna = self.rna.replace('U', 'T').upper()

        self.id = id

        if update:
            self.__update__()

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

    def __update__(self):
        self.__find_reverse_complements__()
        self.__form_peptide_chain__()
        self.__form_reverse_complement_peptide_chain__()
        self.__find_skew__()
        self.__find_suspected_oriC__()
        self.__update_base_counts__()

    def __form_peptide_chain__(self) -> None:
        """Internal method. Defines the peptide chain
        formed by the entire strand"""
        self.peptide_chain = ''
        rna_copy = [x for x in self.rna]
        rna_copy.reverse()

        for _ in range(len(rna_copy) // 3):
            codon = ''
            for _ in range(3):
                codon += rna_copy.pop()
            self.peptide_chain += codon_table[codon]

    def __form_reverse_complement_peptide_chain__(self):
        """Internal method. Defines the reverse complement
        peptide chain formed by the entire DNA strand"""
        self.reverse_complement_peptide_chain = ''
        try:
            self.rna_complement
        except AttributeError:
            self.__find_reverse_complements__()
        rna_copy = [x for x in self.rna_complement]

        for _ in range(len(rna_copy) // 3):
            codon = ''
            for _ in range(3):
                codon += rna_copy.pop()
            self.reverse_complement_peptide_chain += codon_table[codon]

    def __find_skew__(self):
        """Internal method. Defines the skew graph of the DNA"""
        self.skew = [0]
        skew = 0
        for base in self.dna:
            if base == 'C':
                skew -= 1
            elif base == 'G':
                skew += 1
            self.skew.append(skew)

    def __find_suspected_oriC__(self) -> None:
        """Internal method. Defines the index of the suspected oriC"""
        try:
            self.skew
        except AttributeError:
            self.__find_skew__()
        min_skew = min(self.skew)
        self.suspected_oriC = self.skew.index(min_skew)

    def __find_reverse_complements__(self) -> None:
        """Internal method. Finds and assigns the reverse complements
        of the DNA and RNA strands"""
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

    def __update_base_counts__(self) -> None:
        """Internal method. Counts each base in RNA and DNA strands"""
        self.dna_base_count = {'G': self.dna.count('G'),
                               'C': self.dna.count('C'),
                               'T': self.dna.count('T'),
                               'A': self.dna.count('A')}
        self.rna_base_count = {'G': self.rna.count('G'),
                               'C': self.rna.count('C'),
                               'U': self.rna.count('U'),
                               'A': self.rna.count('A')}

    def __find_reading_frames__(self) -> None:
        """Internal method. Finds all reading frames"""
        self.reading_frames = []
        try:
            self.rna_reverse_complement
        except AttributeError:
            self.__find_reverse_complements__()
        for i in range(3):
            beginning = 0 + i
            if i != 0:
                ending = len(self.dna) - (3 - i)
            else:
                ending = len(self.dna)
            self.reading_frames.append(Sequence(self.rna[beginning:ending]))
            self.reading_frames.append(Sequence(self.rna_reverse_complement
                                       [beginning:ending]))

    def __find_open_reading_frames__(self) -> None:
        """Internal method. Finds all open reading frames"""
        try:
            self.reading_frames
        except AttributeError:
            self.__find_reading_frames__()
        self.open_reading_frames = []
        for frame in self.reading_frames:
            start_codons = find_all_indices(frame.rna, 'AUG')
            stop_codons = []
            stop_codons.extend(find_all_indices(frame.rna, 'UAA'))
            stop_codons.extend(find_all_indices(frame.rna, 'UAG'))
            stop_codons.extend(find_all_indices(frame.rna, 'UGA'))

            for start_codon in start_codons:
                for stop_codon in stop_codons:
                    self.open_reading_frames.append(frame.rna
                                                    [start_codon:stop_codon])
        self.open_reading_frames = [Sequence(x, is_rna=True) for x in
                                    self.open_reading_frames]
        for frame in self.open_reading_frames:
            frame.__form_peptide_chain__()
        self.open_reading_frames = [x for x in self.open_reading_frames
                                    if '*' not in x.peptide_chain
                                    and x.peptide_chain != '']

    def find_reverse_palindromes(self, min_length: int,
                                 max_length: int) -> List[Tuple[int, int]]:
        answer = []
        for i in range(min_length, max_length + 1):
            for ii in range(len(self.dna) - i + 1):
                partial_sequence = Sequence(self.dna[ii:ii+i],
                                            update=True)
                if partial_sequence.dna == \
                   partial_sequence.dna_reverse_complement:
                    answer.append((ii + 1, i))
        answer = sorted(answer)
        return answer

    def gc_content(self, round_places: int = 3) -> Decimal:
        """Checks the gc content of the DNA chain

        :param round_places: Number of decimal places to define
        :type round_places: int
        :return: gc_content as a percentage from 0 - 1
        :rtype: Decimal
        """
        numerator = Decimal(self.dna_base_count['G']
                            + self.dna_base_count['C'])
        denominator = Decimal(len(self.dna))
        return round(numerator / denominator, round_places)

    def find_substrings(self, shortest: int, longest: int = 0,
                        *, window: Tuple[int, int] = (0, 0)) -> Dict[str, int]:
        """Finds all unique substrings of any defined length
        inside of a specific window of the DNA strand

        :param shortest: The shortest substrings to find
        :type shorest: int
        :param longest: The longest substrings to find
        :type longest: int
        :param window: The range of the window on the
        DNA strand to search, defaults to the entire DNA strand
        :type window: Tuple[int, int]
        :return: All substrings contained in the DNA strand,
        and the number of times each occurs
        :rtype: Dict[str, int]
        """
        if window == (0, 0):
            window = (0, len(self.dna))
        if longest == 0:
            longest = shortest
        if window[0] > window[1]:
            window = (window[1], window[0])
        if window[1] > len(self.dna):
            print("Window out of bounds. Fixing.")
            window = (window[0], len(self.dna))
        if window[0] < 0:
            print("Window out of bounds. Fixing.")
            window = (0, window[1])

        if longest < shortest:
            shortest, longest = longest, shortest
        if longest > window[1] - window[0]:
            print("Substrings too large. Fixing.")
            longest = window[1] - window[0]

        substring_chains: Dict[str, int] = dict()
        for length in range(shortest, longest+1):
            for i in range(window[0], window[1] - length + 1):
                if self.dna[i:i+length] in substring_chains:
                    substring_chains[self.dna[i:i+length]] += 1
                else:
                    substring_chains[self.dna[i:i+length]] = 1
        substring_chains = dict(sorted(
            substring_chains.items(), key=lambda x: x[1]))
        return substring_chains

    def find_clumps(self, k: int, clump_size: int,
                    window_size: int, search_window:
                    Tuple[int, int] = (0, 0)) -> set[str]:
        """Finds 'clumps' of polymers located inside the DNA strand

        :param k: The length of polymer to search for
        :type k: int
        :param clump_size: The number of identical
        polymers we wish to find in our window
        :type clump_size: int
        :param window_size: The size of the
        SLIDING window. Think, 'I want to find clumps of
        k-mers grouped within window_size base pairs of each other'
        :type window_size: int
        :param search_window: The area of the DNA strand you want to search
        :type search_window: Tuple[int, int]
        :return: All k-mers that appears at
        least clump_size times within a range of window_size base pairs
        :rtype: Set[str]
        """

        if search_window == (0, 0):
            search_window = (0, len(self.dna))

        if search_window[0] > search_window[1]:
            search_window = (search_window[1], search_window[0])
        if search_window[0] < 0:
            print("Search window can't be negative. Fixing.")
            search_window = (0, search_window[1])
        if search_window[1] > len(self.dna):
            print("Can't search outside DNA. Fixing.")
            search_window = (search_window[0], len(self.dna))
        if window_size < search_window[1] - search_window[0]:
            print("Window size can't be less than search window size. Fixing.")
            window_size = search_window[1] - search_window[1]

        substring_chains: List[str] = []
        for i in range(search_window[0], search_window[1] - window_size + 1):
            substrings = self.find_substrings(k, window=(i, i + window_size))
            for ke, v in substrings.items():
                if v >= clump_size:
                    substring_chains.append(ke)
        substring_chains_set = set(substring_chains)
        return substring_chains_set

    def splice(self, intron: str) -> None:
        """Method that removes the specified intron from the RNA strand
        Immediately calls __update_after_changes__ when done

        :param intron: The intron which should be removed
        :type intron: str
        """
        while intron in self.rna:
            self.rna = self.rna.replace(intron, '')
        self.__update__()

    def find_motif(self, dna_chain: str) -> List[int]:
        """Finds all occurences of dna_chain inside the DNA strand

        :param dna_chain Chain to find locations of
        :type dna_chain: str
        :return: A list of all starting indices of dna_chain
        :rtype: List[int]
        """
        locations = []
        m_len = len(dna_chain)
        for i in range(0, len(self.dna) - m_len):
            if dna_chain == self.dna[i:i+m_len]:
                locations.append(i)

        return locations

    def find_approximate_match(self, pattern: str, distance: int) -> List[int]:
        """Same thing as find_motif, but finds similar
        patterns as well as exact matches

        :param pattern: The pattern to find in the DNA strand
        :type pattern: str
        :param distance: The allowable hamming distance
        :type distance: int
        :return: A list of all starting indices of pattern
        :rtype: List[int]
        """
        locations: List[int] = []
        for i in range(len(self.dna) - len(pattern) + 1):
            if hamming_distance(self.dna[i:i+len(pattern)],
                                pattern) <= distance:
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
amino_acids_string = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRR'
amino_acids_string += 'VVVVAAAADDEEGGGG'
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
amino_acid_mono_mass: Dict[str, float] = dict()
amino_acid_average_mass: Dict[str, float] = dict()

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


def ok_adjacency_list(k: int, sequences: List[Sequence]) -> List[str]:
    answer = ''
    for loop1 in sequences:
        for loop2 in sequences:
            if loop1.dna[-k:] == loop2.dna[:k] and loop1.dna != loop2.dna:
                answer += f'{loop1.id} {loop2.id}\n'

    return answer

def find_all_indices(s, ch):
    return [i for i in range(0, len(s), len(ch)) if s[i:i+len(ch)] == ch]

def concat_sequences(sequences: List[Sequence]) -> Sequence:
    """Concatenates multiple sequences into a single sequence

    :param sequences: All sequences to be concatenated
    :type sequences: List[Sequence]
    :return: A single sequence
    :rtype: Sequence
    """
    sequence_string = ''
    for i, sequence in enumerate(sequences):
        sequence_string += sequence.dna

    return Sequence(sequence_string)


def __generate_random_sequence__(gc: int, seq_length: int) -> Sequence:
    """Internal function. Generates a single random sequence

    :param gc: The (approximate) gc content of the
    sequence as a percent, 0 - 100
    :type gc: int
    :param seq_length: The length of the sequences to be generated
    :type seq_length: int
    :return: A random sequence
    :rtype: Sequence
    """
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
                              *, equilength:
                              bool = True) -> Union[list[Sequence], Sequence]:
    """Generates a list of random sequences, or a single random sequence

    :param num_sequences: The number of sequences to generate
    :type num_sequences: int
    :param gc: The (approximate) gc content of the sequences,
    minimum and maximum
    :type gc: Tuple[float, float]
    :param seq_length: The length of of sequences, minimum and maximum
    :type seq_length: Tuple[int, int]
    :param equilength: Whether the sequences should all
    be the same length, default to true
    :type equilength: bool
    :return: A list of randomly generated sequences, or a single sequence
    :rtype: Union[List[Sequence], Sequence]
    """
    sequences = []
    if equilength:
        sequence_length = randint(seq_length[0], seq_length[1])

    for _ in range(num_sequences):
        gc_content = randint(int(gc[0] * 100), int(gc[1] * 100))

        if not equilength:
            sequence_length = randint(seq_length[0], seq_length[1])

        sequences.append(__generate_random_sequence__(
            gc_content, sequence_length))

    if len(sequences) == 1:
        return sequences[0]
    else:
        return sequences


def read(filename: str, prefix: str = 'data/') -> list[str]:
    """Reads a file. Mostly for use with rosalind problems

    :param filename: Filename
    :type filename: str
    :param prefix: Prefix for filename
    :type prefix: str
    :return: Contents of the file, each line as an entry in a list
    :trype: List[str]
    """
    with open(f'{prefix}{filename}') as f:
        input_data = f.readlines()
        input_data = [x.strip() for x in input_data]
        return input_data


def read_fasta_data(filename: str, prefix: str = 'data/',
                    update: bool = False) -> list[Sequence]:
    """Reads a file or fasta sequences. Mostly for use with rosalind problems

    :param filename: Filename
    :type filename: str
    :param prefix: Prefix for filename
    :type prefix: str
    :param update: Update immediately? Defaults to False
    :type update: Bool
    :return: Content of the file
    :trype: List[Sequence]
    """
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
        sequences.append(Sequence(v, k, update=update))

    return sequences


def read_uniprot(access_ids: list[str], quiet=True) -> list[Protein]:
    """Fetches fasta data from uniprot.org

    :param access_ids: Ids to locate in the uniprot database
    :type access_ids: List[str]
    :return: All proteins matching access_ids
    :rtype: List[Protein]
    """
    base_url = 'https://rest.uniprot.org/uniprotkb/'
    proteins: list[Protein] = []
    for id in access_ids:
        real_id = id.split('_', 1)[0]
        url = base_url + real_id
        if not quiet:
            print(f'Fetching {url}')
        downloaded_data = (requests.get(url, allow_redirects=True).text)
        downloaded_data = json.loads(downloaded_data)
        downloaded_data = downloaded_data['sequence']['value']
        proteins.append(Protein(downloaded_data, id))
    return proteins


def hamming_distance(s: Union[Sequence, str], t: Union[Sequence, str]) -> int:
    """Finds the hamming distance between two DNA strands

    :param s: A DNA strand
    :type s: Union[Sequence, str]
    :param t: A DNA strand
    :type t: Union[Sequence, str]
    :return: The hamming distance between s and t
    :rtype: int
    """
    hamming = 0
    if isinstance(s, Sequence):
        s = s.dna
    if isinstance(t, Sequence):
        t = t.dna
    for base_s, base_t in zip(s, t):
        if base_s != base_t:
            hamming += 1
    return hamming


def find_consensus(consensus_matrix: Dict[str, List[int]]) -> Sequence:
    """Finds the 'average' DNA strand from a consensus_matrix

    :param consensus_matrix: A consensus matrix
    created from a list of Sequences.
    The return parameter of find_consensus_matrix below
    :type consensus_matrix: dict[str, List[int]]
    :return: A consensus Sequence
    :rtype: Sequence
    """
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
    """Finds a consensus matrix, as defined at
    https://rosalind.info/problems/cons/
    from a list of Sequences

    :param sequences: The sequences to extract the consensus matrix from
    :type sequences: List[Sequence]
    :return: A consensus matrix
    :rtype: Dict[str, List[int]]
    """
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

# NOTE (Tess): This is very slow. A list of k <= 100 sequences
# each of length <= 1000 takes 7 seconds to compute

def find_longest_common_substring(sequences: List[Sequence]) -> str:
    longest_sequence = Sequence('')
    max_length = 0
    for sequence in sequences:
        if len(sequence) > max_length:
            longest_sequence = sequence
            max_length = len(sequence)
    longest_substring = ''
    substring_length = 0

    substrings = longest_sequence.find_substrings(1, len(longest_sequence))
    for substring in substrings.keys():
        in_all_sequences = True
        for sequence in sequences:
            if substring not in sequence.dna:
                in_all_sequences = False
                break
        if in_all_sequences:
            if len(substring) > substring_length:
                longest_substring = substring
                substring_length = len(substring)

    return longest_substring


def check_sequence_lengths(sequences: list[Sequence]) -> None:
    """Ensures a list of sequences are all the same length

    :param sequences: List of Sequences to check length of
    :type sequences: List[Sequences]
    :raises SequenceLengthException: If all sequences are not equal length
    """
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
