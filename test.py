import unittest
from Genetics_Toolbox import *

sequences = read_fasta_data('testdata.data', '', update=True)
basic_test_sequence = Sequence('AGCTGGCTA')

class TestProteinMotifCreation(unittest.TestCase):
    def runTest(self):
        protein_motif = ProteinMotif('KM[YVW]{YVW}[YVW]')
        self.assertEqual(protein_motif.positions, [['K'], ['M'], ['Y', 'V', 'W'], ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T'], ['Y', 'V', 'W']])
        self.assertTrue(protein_motif.match('KMWLV'))
        self.assertTrue(protein_motif.match('KMYKY'))
        self.assertFalse(protein_motif.match('KMWLVL'))
        self.assertFalse(protein_motif.match('KMYYY'))

class TestProteinCreation(unittest.TestCase):
    def runTest(self):
        protein = read_uniprot(["B5ZC00"])[0]
        self.assertEqual(protein.find_motifs(n_glycosylation),
                    [85, 118, 142, 306, 395],
                    "Protein motifs are incorrect")
        
class TestSequenceUpdates(unittest.TestCase):
    def runTest(self):
        update_test_sequence = generate_random_sequences(1, (0.0,1.0), (10**5, 10**6))
        update_test_sequence.__update__()

class TestPeptideChains(unittest.TestCase):
    def runTest(self):
        test_sequence = Sequence('TTTTTCTTATTGCTTCTCCTACTGATTATCATAATGGTTGTCGTAGTGTCTTCCTCATCGCCTCCCCCACCGACTACCACAACGGCTGCCGCAGCGTATTACTAATAGCATCACCAACAGAATAACAAAAAGGATGACGAAGAGTGTTGCTGATGGCGTCGCCGACGGAGTAGCAGAAGGGGTGGCGGAGGG')
        test_sequence.__update__()
        
        self.assertEqual('FFLLLLLLIIIMVVVVSSSSPPPPTTTTAAAAYY**HHQQNNKKDDEECC*WRRRRSSRRGGGG', test_sequence.peptide_chain)
        self.assertEqual('PSATPSATPSATPSATLFVILFVILLVMLLVIRCGSRCGSRWGRR*GRHYDNHYDNQ*EKQ*EK', test_sequence.reverse_complement_peptide_chain)

class TestSkew(unittest.TestCase):
    pass
        

class TestGCContent(unittest.TestCase):
    def runTest(self):
        for i, sequence in enumerate(sequences):
            if i == 0:
                gc = .476
            elif i == 1:
                gc = .67
            elif i == 2:
                gc = .508
            elif i == 3:
                gc = .622
            elif i == 4:
                gc = .611

            self.assertEqual(float(sequence.gc_content(3)), gc,
                             'GC content is incorrect')

class TestBaseCounts(unittest.TestCase):
    def runTest(self):
        for i, sequence in enumerate(sequences):
            if i == 0:
                dna_base_count = {'G': 1019,
                                  'C': 1050,
                                  'T': 1122,
                                  'A': 1159}
            elif i == 1:
                dna_base_count = {'G': 1494,
                                  'C': 1406,
                                  'T': 734,
                                  'A': 692}
            elif i == 2:
                dna_base_count = {'G': 1261,
                                  'C': 1282,
                                  'T': 1235,
                                  'A': 1231}
            elif i == 3:
                dna_base_count = {'G': 2270,
                                  'C': 2242,
                                  'T': 1358,
                                  'A': 1383}
            elif i == 4:
                dna_base_count = {'G': 2499,
                                  'C': 2499,
                                  'T': 1620,
                                  'A': 1558}
                
            self.assertEqual(dna_base_count, sequence.dna_base_count,
                             'DNA base count is incorrect.')

class TestSubstrings(unittest.TestCase):
    def runTest(self):
        subs = basic_test_sequence.find_substrings(8, 9)
        self.assertEqual(subs, {'AGCTGGCT': 1, 'GCTGGCTA': 1, 'AGCTGGCTA': 1},
                         "Substrings are incorrect")

unittest.main()
