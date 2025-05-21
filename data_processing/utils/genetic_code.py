from collections import Counter
import math
import numpy as np


class GeneticCode:

    syn_short = {
        "GCA": "A",
        "GCC": "A",
        "GCG": "A",
        "GCT": "A",
        "AGA": "R",
        "AGG": "R",
        "CGA": "R",
        "CGC": "R",
        "CGG": "R",
        "CGT": "R",
        "AAC": "N",
        "AAT": "N",
        "GAC": "D",
        "GAT": "D",
        "TGC": "C",
        "TGT": "C",
        "CAA": "Q",
        "CAG": "Q",
        "GAA": "E",
        "GAG": "E",
        "GGA": "G",
        "GGC": "G",
        "GGG": "G",
        "GGT": "G",
        "CAC": "H",
        "CAT": "H",
        "ATA": "I",
        "ATC": "I",
        "ATT": "I",
        "CTA": "L",
        "CTC": "L",
        "CTG": "L",
        "CTT": "L",
        "TTA": "L",
        "TTG": "L",
        "AAA": "K",
        "AAG": "K",
        "ATG": "M",   # non-synonymous
        "TTC": "F",
        "TTT": "F",
        "CCA": "P",
        "CCC": "P",
        "CCG": "P",
        "CCT": "P",
        "AGC": "S",
        "AGT": "S",
        "TCA": "S",
        "TCC": "S",
        "TCG": "S",
        "TCT": "S",
        "ACA": "T",
        "ACC": "T",
        "ACG": "T",
        "ACT": "T",
        "TGG": "W",   # non-synonymous
        "TAC": "Y",
        "TAT": "Y",
        "GTA": "V",
        "GTC": "V",
        "GTG": "V",
        "GTT": "V",
        "TAA": ".",
        "TAG": ".",
        "TGA": "."
    }

    syn = {
        "GCA": "Ala",
        "GCC": "Ala",
        "GCG": "Ala",
        "GCT": "Ala",
        "AGA": "Arg",
        "AGG": "Arg",
        "CGA": "Arg",
        "CGC": "Arg",
        "CGG": "Arg",
        "CGT": "Arg",
        "AAC": "Asn",
        "AAT": "Asn",
        "GAC": "Asp",
        "GAT": "Asp",
        "TGC": "Cys",
        "TGT": "Cys",
        "CAA": "Gln",
        "CAG": "Gln",
        "GAA": "Glu",
        "GAG": "Glu",
        "GGA": "Gly",
        "GGC": "Gly",
        "GGG": "Gly",
        "GGT": "Gly",
        "CAC": "His",
        "CAT": "His",
        "ATA": "Ile",
        "ATC": "Ile",
        "ATT": "Ile",
        "CTA": "Leu",
        "CTC": "Leu",
        "CTG": "Leu",
        "CTT": "Leu",
        "TTA": "Leu",
        "TTG": "Leu",
        "AAA": "Lys",
        "AAG": "Lys",
        "ATG": "Met",  # non-synonymous
        "TTC": "Phe",
        "TTT": "Phe",
        "CCA": "Pro",
        "CCC": "Pro",
        "CCG": "Pro",
        "CCT": "Pro",
        "AGC": "Ser",
        "AGT": "Ser",
        "TCA": "Ser",
        "TCC": "Ser",
        "TCG": "Ser",
        "TCT": "Ser",
        "ACA": "Thr",
        "ACC": "Thr",
        "ACG": "Thr",
        "ACT": "Thr",
        "TGG": "Trp",  # non-synonymous
        "TAC": "Tyr",
        "TAT": "Tyr",
        "GTA": "Val",
        "GTC": "Val",
        "GTG": "Val",
        "GTT": "Val",
        "TAA": "Stop",
        "TAG": "Stop",
        "TGA": "Stop"
    }

    amino_acid_dict = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
    }

    syn_count = {
        "GCA": 4,
        "GCC": 4,
        "GCG": 4,
        "GCT": 4,
        "AGA": 6,
        "AGG": 6,
        "CGA": 6,
        "CGC": 6,
        "CGG": 6,
        "CGT": 6,
        "AAC": 2,
        "AAT": 2,
        "GAC": 2,
        "GAT": 2,
        "TGC": 2,
        "TGT": 2,
        "CAA": 4,
        "CAG": 4,
        "GAA": 4,
        "GAG": 4,
        "GGA": 4,
        "GGC": 4,
        "GGG": 4,
        "GGT": 4,
        "CAC": 2,
        "CAT": 2,
        "ATA": 3,
        "ATC": 3,
        "ATT": 3,
        "CTA": 6,
        "CTC": 6,
        "CTG": 6,
        "CTT": 6,
        "TTA": 6,
        "TTG": 6,
        "AAA": 2,
        "AAG": 2,
        "ATG": 1,  # non-synonymous
        "TTC": 2,
        "TTT": 2,
        "CCA": 4,
        "CCC": 4,
        "CCG": 4,
        "CCT": 4,
        "AGC": 6,
        "AGT": 6,
        "TCA": 6,
        "TCC": 6,
        "TCG": 6,
        "TCT": 6,
        "ACA": 4,
        "ACC": 4,
        "ACG": 4,
        "ACT": 4,
        "TGG": 1,  # non-synonymous
        "TAC": 2,
        "TAT": 2,
        "GTA": 4,
        "GTC": 4,
        "GTG": 4,
        "GTT": 4,
        "TAA": 3,
        "TAG": 3,
        "TGA": 3
    }

    transcription_trans_table = str.maketrans({"A": "T", "T": "A", "G": "C", "C": "G"})

    @staticmethod
    def complementaryRNA(seq):
        return seq.translate(GeneticCode.transcription_trans_table)[::-1]

    @staticmethod
    def translate(codon: str, long_aa_name: bool = True):
        dictionary = GeneticCode.syn if long_aa_name else GeneticCode.syn_short
        return dictionary[codon.replace("U", "T")]
    
    @staticmethod
    def synonymous_of(codon: str):
        am = GeneticCode.translate(codon)
        return [c for c in GeneticCode.syn if GeneticCode.syn[c] == am]
    
    @staticmethod
    def all_aminoacids(long_aa_names: bool = True):
        dictionary = GeneticCode.syn if long_aa_names else GeneticCode.syn_short
        return sorted(set(dictionary.values()))
    
    @staticmethod
    def reverse_translate(aminoacid: str, long_aa_names: bool = True):
        if long_aa_names:
            return sorted([k for k,v in GeneticCode.syn.items() if v == aminoacid])
        else:
            return sorted([k for k,v in GeneticCode.syn_short.items() if v == aminoacid])

    #            start (non syn), another non-syn, 3 stop codons                         
    ignored_codons = ["ATG", "TGG", "TAA", "TAG", "TGA"]

    valid_codons = ["GCA", "GCC", "GCG", "GCT", "AGA", "AGG", "CGA", "CGC", "CGG", "CGT", "AAC", "AAT", "GAC", "GAT", "TGC", "TGT", "CAA", "CAG", "GAA", "GAG", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC", "ATT", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "AAA", "AAG", "TTC", "TTT", "CCA", "CCC", "CCG", "CCT", "AGC", "AGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC", "ACG", "ACT", "TAC", "TAT", "GTA", "GTC", "GTG", "GTT"]

    start_codon = "ATG"

    stop_codons = ["TAA", "TAG", "TGA"]

    valid_characters = {"A", "C", "T", "G"}

    all_codons = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]

    even_rscu = {
        'GCA': 0.25,
        'GCC': 0.25,
        'GCG': 0.25,
        'GCT': 0.25,
        'AGA': 0.16666666666666666,
        'AGG': 0.16666666666666666,
        'CGA': 0.16666666666666666,
        'CGC': 0.16666666666666666,
        'CGG': 0.16666666666666666,
        'CGT': 0.16666666666666666,
        'AAC': 0.5,
        'AAT': 0.5,
        'GAC': 0.5,
        'GAT': 0.5,
        'TGC': 0.5,
        'TGT': 0.5,
        'CAA': 0.5,
        'CAG': 0.5,
        'GAA': 0.5,
        'GAG': 0.5,
        'GGA': 0.25,
        'GGC': 0.25,
        'GGG': 0.25,
        'GGT': 0.25,
        'CAC': 0.5,
        'CAT': 0.5,
        'ATA': 0.3333333333333333,
        'ATC': 0.3333333333333333,
        'ATT': 0.3333333333333333,
        'CTA': 0.16666666666666666,
        'CTC': 0.16666666666666666,
        'CTG': 0.16666666666666666,
        'CTT': 0.16666666666666666,
        'TTA': 0.16666666666666666,
        'TTG': 0.16666666666666666,
        'AAA': 0.5,
        'AAG': 0.5,
        'ATG': 1.0,
        'TTC': 0.5,
        'TTT': 0.5,
        'CCA': 0.25,
        'CCC': 0.25,
        'CCG': 0.25,
        'CCT': 0.25,
        'AGC': 0.16666666666666666,
        'AGT': 0.16666666666666666,
        'TCA': 0.16666666666666666,
        'TCC': 0.16666666666666666,
        'TCG': 0.16666666666666666,
        'TCT': 0.16666666666666666,
        'ACA': 0.25,
        'ACC': 0.25,
        'ACG': 0.25,
        'ACT': 0.25,
        'TGG': 1.0,
        'TAC': 0.5,
        'TAT': 0.5,
        'GTA': 0.25,
        'GTC': 0.25,
        'GTG': 0.25,
        'GTT': 0.25,
        'TAA': 0.3333333333333333,
        'TAG': 0.3333333333333333,
        'TGA': 0.3333333333333333
    }
    
    @staticmethod
    def rscu(cds_list: list):
        codon_counter = Counter()
        aa_counter = Counter()
        frequencies = dict()
        # consider multiple CDS
        for cds in cds_list:
            # extract triplets
            start_p = np.arange(math.floor(len(cds)/3))*3
            stop_p = start_p+3
            # divide in triplets
            codons = [cds[i:e] for i,e in zip(start_p, stop_p)]
            # remove unknown triplets (untranslatable chatacters like n et simila)
            codons = [c for c in codons if set(c) <= GeneticCode.valid_characters] 
            # remove stop codons and non-synonyous codons
            codons = [c for c in codons if c not in GeneticCode.ignored_codons] 
            # update counter codon usage
            codon_counter.update(codons)
            translation = [GeneticCode.translate(c) for c in codons]
            aa_counter.update(translation)
        # number of times the codon is usee / number of times the corresponding amino acid is translated
        frequencies = {c: codon_counter[c] / aa_counter[GeneticCode.translate(c)] for c in codon_counter}
        # rscu 
        rscu = {c: ((frequencies[c] / GeneticCode.even_rscu[c]) if c in frequencies else 1) for c in GeneticCode.valid_codons}
        return rscu

    @staticmethod
    def cds_ranges(sequence):
        orf_start = sequence.find(GeneticCode.start_codon)
        while orf_start >= 0:
            # divide trailing sequence in triplets
            seq_after_start_codon = sequence[orf_start+3:]
            start_p = np.arange(math.floor(len(seq_after_start_codon)/3))*3
            stop_p = start_p+3
            codons = [seq_after_start_codon[i:e] for i,e in zip(start_p, stop_p)]
            # find nearest stop codon
            orf_ends = [] 
            for sc in GeneticCode.stop_codons:
                try:
                    orf_ends.append(codons.index(sc))
                except ValueError:
                    continue    # this stop codon is not present
            try:
                orf_end = min(orf_ends) # this is the triplet number in "codons" list
            except ValueError:
                break   # if there is not any stop codon
            orf_end = orf_start + 3 + orf_end*3 + 3     # this is the right end character index
            yield orf_start, orf_end
            orf_start = sequence.find(GeneticCode.start_codon, orf_end)
    
    @staticmethod
    def extract_cds_list(sequence):
        return [sequence[orf_start:orf_end] for (orf_start,orf_end) in GeneticCode.cds_ranges(sequence)]

    @staticmethod
    def codon_frequency(cds_list):
        codons = []
        for sequence in cds_list:
            start_p = np.arange(math.floor(len(sequence)/3))*3
            stop_p = start_p+3
            codons += [sequence[i:e] for i,e in zip(start_p, stop_p)]
        codon_counter = Counter(codons)
        total_codons = codon_counter.total()
        codon_freq = {k: (codon_counter[k] / total_codons if k in codon_counter else 0.0) for k in GeneticCode.valid_codons}
        return codon_freq

    @staticmethod
    def codon_count(cds_list):
        codons = []
        for sequence in cds_list:
            start_p = np.arange(math.floor(len(sequence) / 3)) * 3
            stop_p = start_p + 3
            codons += [sequence[i:e] for i, e in zip(start_p, stop_p)]
        codon_counter = Counter(codons)
        codon_freq = {k: (codon_counter[k] if k in codon_counter else 0) for k in
                      GeneticCode.valid_codons}
        return codon_freq

    @staticmethod
    def translate_cds(cds_list, long_aa_names: bool = True):
        """
        Transform a list of coding sequences into the corresponding proteins.
        """
        for cds in cds_list:
            # extract triplets
            start_p = np.arange(math.floor(len(cds) / 3)) * 3
            stop_p = start_p + 3
            # divide in triplets
            codons = [cds[i:e] for i, e in zip(start_p, stop_p)]
            # remove unknown triplets (untranslatable characters like n et simila)
            codons = [c for c in codons if set(c) <= GeneticCode.valid_characters]
            # remove stop codons
            if codons[-1] in GeneticCode.stop_codons:
                codons = codons[:-1]
            # remove stop codons and non-synonyous codons
            # codons = [c for c in codons if c not in GeneticCode.ignored_codons]
            translation = "".join([GeneticCode.translate(c, long_aa_names) for c in codons])
            yield translation

    @staticmethod
    def AAnames2letters(sequence_of_aa_names):
        """
        Transform a sequence of aminoacids with three letters encoding into single letter encoding.
        """
        if type(sequence_of_aa_names) is str:
            return "".join(
                [GeneticCode.amino_acid_dict[sequence_of_aa_names[x:x + 3]] for x in range(0, len(sequence_of_aa_names), 3)])
        else:
            return ["".join([GeneticCode.amino_acid_dict[y[x:x + 3]] for x in range(0, len(y), 3)]) for y in sequence_of_aa_names]

    @staticmethod
    def motif_coordinates(motif="RRKR", *sequences):
        """
        Return a list of coordinate ranges (left included, right excluded) where the motif is found in evey given
        sequence. Empty list otherwise.
        """
        sites = []
        for p in sequences:
            idx = p.find(motif)
            while idx != -1:
                sites.append((idx, idx + len(motif)))
                idx = p.find(motif, start=idx + len(motif))
        return sites

    @staticmethod
    def aa_frequency(protein_list, long_aa_names: bool = False):
        if long_aa_names:
            raise NotImplementedError
        aa_counter = Counter("".join(protein_list))
        total_aa = aa_counter.total()
        all_aa = GeneticCode.all_aminoacids(long_aa_names=long_aa_names)
        return {k: ((aa_counter[k] / total_aa) if k in aa_counter else 0) for k in all_aa}

    @staticmethod
    def aa_count(protein_list, long_aa_names: bool = False):
        if long_aa_names:
            raise NotImplementedError
        aa_counter = Counter("".join(protein_list))
        all_aa = GeneticCode.all_aminoacids(long_aa_names=long_aa_names)
        aa_freq = {k: (aa_counter[k] if k in aa_counter else 0) for k in all_aa}
        return aa_freq
    
    @staticmethod
    def extract_triplets(cds):
        # extract triplets
        start_p = np.arange(math.floor(len(cds)/3))*3
        stop_p = start_p+3
        # divide in triplets
        return [cds[i:e] for i,e in zip(start_p, stop_p)]
    
    @staticmethod
    def locate_codons_coordinates(cds: str, codons: list):
        """
        Returns a dictonary with the argument 'codons' as keys. The value of every codon is a list of coordinates corresponding to the positions of the codon in the cds sequence.
        """
        triplets = np.array(GeneticCode.extract_triplets(cds))
        return {cod:np.where(triplets == cod)[0] * 3 for cod in codons}
    
    dinucleotides = ("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
    _dinucleotides_subsymbols = (('AA', 'A', 'A'), ('AC', 'A', 'C'), ('AG', 'A', 'G'), ('AT', 'A', 'T'), ('CA', 'C', 'A'), ('CC', 'C', 'C'), ('CG', 'C', 'G'), ('CT', 'C', 'T'), ('GA', 'G', 'A'), ('GC', 'G', 'C'), ('GG', 'G', 'G'), ('GT', 'G', 'T'), ('TA', 'T', 'A'), ('TC', 'T', 'C'), ('TG', 'T', 'G'), ('TT', 'T', 'T'))
    
    @staticmethod
    def dinucleotides_count(cds_list: list):
        nucleotides_count = Counter()
        _dinucleotides_count = Counter()
        nuc_length = 0
        for cds in cds_list:
            cds = cds.upper().replace("U","T")
            nuc_length += len(cds)
            nucleotides_count.update(cds)
            _dinucleotides_count.update({dn: cds.count(dn) for dn in GeneticCode.dinucleotides})
        dinuc_length = nuc_length - len(cds_list)   # number of existing dinucletides in sequence of length N is N-1
        # at end of loop
        ratio_operands = [(_dinucleotides_count[dn] / dinuc_length, nucleotides_count[n1]*nucleotides_count[n2] / nuc_length**2) for dn,n1,n2 in GeneticCode._dinucleotides_subsymbols]
        _dinucleotides_freq = [op[0] / op[1] if op[1] != 0 else 0 for op in ratio_operands]
        return dict(zip(GeneticCode.dinucleotides,_dinucleotides_freq))

