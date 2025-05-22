""" This is a utility module for calling variants from sequences aligned using Augur align. 
   
    Steps: 
        - Align your reference sequence to one or more target sequences using:
            
            augur align --reference-sequence refseq.fasta --sequences target_seq.fasta --output aligned_seq.fasta 
        
            Two files are expected as a result: 
            - aligned_seq.fasta == a fasta with one aligned sequence for every sequence in target_seq.fasta. The alignment allows reading substitutions and deletions.
            - aligned_seq.insertions.csv == a tabular file containing the insertions of the sequences in target_seq.fasta with respect to the reference.
        
            Documentation (https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/align.html)

        - Create an instance of VariantCaller

            requires:
            - the reference sequence as a string 
            - (optional) the path to the file aligned_seq.insertions.csv (if you want the insertions)

        - Call the variants for each individual sequence

            target_changes = variant_caller.call_variants()

            requires:
            - the aligned sequence as a string (the one you can read from aligned_seq.fasta)
            - (optional) the name of the sequence as reported in the header of each sequence in target_seq.fasta. If None or empty string, insertions are not returned
            - (... other optional arguments) 

            The output object is an instance of AlignmentChanges

    Example: see the __main__ function.
    
    Author Name: Tommaso Alfonsi
    Version 1.0.0

"""

import numpy as np 
import pandas as pd
import miniFasta    # optional dependency used in __main__ for easily reading fasta
from typing import Optional, List, Tuple
from tqdm import tqdm


class AlignmentChanges:
    def __init__(self) -> None:
        self.subs: Optional[list] = None
        self.subs_pos: Optional[np.array] = None

        self.dels: Optional[list] = None
        self.dels_pos: Optional[np.array] = None

        self.ins: Optional[list] = None
        self.ins_pos: Optional[np.array] = None

    def __str__(self) -> str:
        out = "dels: " + str(self.dels) + " " + str(type(self.dels))+ "\n"
        out += "dels_pos: " + str(self.dels_pos) + " " + str(type(self.dels_pos))+ "\n"
        out += "subs: " + str(self.subs) + " " + str(type(self.subs))+ "\n"
        out += "subs_pos: " + str(self.subs_pos) + " " + str(type(self.subs_pos))+ "\n"
        out += "ins: " + str(self.ins) + " " + str(type(self.ins))+ "\n"
        out += "ins_pos: " + str(self.ins_pos) + " " + str(type(self.ins_pos))+ "\n"
        return out
        

class VariantCaller:

    def __init__(self, refseq: str, insertions_file_path: str = None) -> None:
        self.refseq = np.array(list(refseq))
        self.refseq_len = len(refseq)
        self.insertions_df = VariantCaller._read_insertions_file(insertions_file_path) if insertions_file_path else None

    def call_variants(self, name: str, aligned: str, aggregate_close_changes: bool, one_based: bool, strip_edge_insertions: bool = False) -> AlignmentChanges:
        """
        name: name of the aligned sequence in the insertions file. If provided, returns the insertions, otherwise can be an empty string or None.
        aligned: the aligned sequence as string 
        aggregate_close_changes: if False, every change has length 1. If True, changes in close positions are aggregated. Insertions are always aggregated.
        one_based: whether to return changes as 1-based or 0-based.
        strip_edge_insertions: drop insertions at sequence start and end positions (it may be faster to do it in bulk after this method)
        Return a AlignmentChanges.
        """
        aligned = np.array(list(aligned))
        all_changes_pos = np.argwhere(aligned != self.refseq)

        del_changes_pos_mask = aligned[all_changes_pos] == "-"
        del_changes_pos = all_changes_pos[del_changes_pos_mask]
        del_changes_refs = self.refseq[del_changes_pos]

        sub_changes_pos = all_changes_pos[~del_changes_pos_mask]
        sub_changes_refs = self.refseq[sub_changes_pos]
        sub_changes_alts = aligned[sub_changes_pos]

        if aggregate_close_changes:
            del_changes_pos, del_changes_refs = VariantCaller._aggregate_subs(del_changes_pos, del_changes_refs)
            sub_changes_pos, sub_changes_refs, sub_changes_alts = VariantCaller._aggregate_subs(sub_changes_pos, sub_changes_refs, sub_changes_alts)
            # insertions come already aggregated

        if name and self.insertions_df is not None:
            try:
                seq_ins_series = self.insertions_df.loc[name]
            except KeyError:
                # declare no insertions
                ins_changes_pos, ins_changes_alt = np.empty(0, dtype=int), np.empty(0, dtype=str)
            else:
                ins_changes_idx = seq_ins_series.values.nonzero()[0]
                ins_changes_pos = self.insertions_df.columns.values[ins_changes_idx].astype(int)
                if strip_edge_insertions:
                    ins_changes_pos = ins_changes_pos[(ins_changes_pos > 0) & (ins_changes_pos < self.refseq_len)]  # drop insertions at edges
                    try:
                        ins_changes_alt = seq_ins_series.loc[ins_changes_pos].values    
                    except:
                        print("Exception with sequence ", name)
                        raise
                else:
                    ins_changes_alt = seq_ins_series.iloc[ins_changes_idx].values
        else:
            # declare no insertions
            ins_changes_pos, ins_changes_alt = np.empty(0, dtype=int), np.empty(0, dtype=str)

        if one_based:
            del_changes_pos += 1
            sub_changes_pos += 1
            ins_changes_pos += 1

        ac =  AlignmentChanges()
        ac.dels_pos = del_changes_pos
        ac.dels = VariantCaller._make_dels_str_change(del_changes_pos, del_changes_refs)

        ac.subs_pos = sub_changes_pos
        ac.subs = VariantCaller._make_subs_str_change(sub_changes_pos, sub_changes_refs, sub_changes_alts)

        ac.ins_pos = ins_changes_pos
        ac.ins = VariantCaller._make_ins_str_change(ins_changes_pos, ins_changes_alt)
        return ac
    
    @staticmethod
    def _make_dels_str_change(del_changes_pos, del_changes_refs) -> List[str]:
            assert len(del_changes_pos) == len(del_changes_refs)
            dels = np.column_stack((
                del_changes_pos.astype(str), 
                np.repeat("_", len(del_changes_pos)), 
                del_changes_refs, 
                np.repeat("|_", len(del_changes_pos))))
            dels_str = [''.join(dels[x,:]) for x in np.arange(dels.shape[0])]
            return dels_str

    @staticmethod
    def _make_subs_str_change(sub_changes_pos, sub_changes_refs, sub_changes_alt) -> List[str]:
            assert len(sub_changes_pos) == len(sub_changes_refs) == len(sub_changes_alt)
            subs = np.column_stack((
                sub_changes_pos.astype(str), 
                np.repeat("_", len(sub_changes_pos)), 
                sub_changes_refs, 
                np.repeat("|", len(sub_changes_alt)), 
                sub_changes_alt))
            subs_str = [''.join(subs[x,:]) for x in np.arange(subs.shape[0])]
            return subs_str
    
    @staticmethod
    def _aggregate_subs(changes_pos, changes_refs, changes_alts = None) -> Tuple[np.array]:
        """
        changes_pos: coordinates of changes
        changes_refs: ref. letter of single changes
        changes_alts: alt. letter of single changes. If missing, it is assumed that changes are deletions and have no alt. alleles
        """
        assert len(changes_pos) == len(changes_refs)
        if changes_alts is not None: 
            assert len(changes_alts) == len(changes_refs)
        changes_start_pos = []                          # coordinates of aggregated changes
        changes_aggr_refs = []                          # ref. letters of aggregated changes 
        changes_aggr_alts = []                          # alt. letters of aggregated changes
        aggregated_start_stop_idx = []                  # start (included) and stop (excluded) index of aggregated changes
        idx, num_changes = 0, len(changes_pos)
        while(idx < num_changes):
            pos = changes_pos[idx]                      # coordinate of one deletion
            changes_start_pos.append(pos)               
            
            next_idx = idx + 1                              
            if next_idx < num_changes:                          # check for contiguous dels only if there's another del after this 
                next_pos = changes_pos[next_idx]
                while(next_idx < num_changes and (pos + 1) == next_pos):   # loop until next del coord. == this del coord +1
                    pos = next_pos
                    next_idx += 1
                    try:
                        next_pos = changes_pos[next_idx]   # fails if inned_idx is >= num_del
                    except IndexError:
                        break
            aggregated_start_stop_idx.append((idx, next_idx))
            idx = next_idx
        changes_aggr_refs = ["".join(changes_refs[a:b]) for a,b in aggregated_start_stop_idx]
        if changes_alts is None:
            return np.array(changes_start_pos), np.array(changes_aggr_refs)
        else:
            changes_aggr_alts = ["".join(changes_alts[a:b]) for a,b in aggregated_start_stop_idx]
            return np.array(changes_start_pos), np.array(changes_aggr_refs), np.array(changes_aggr_alts)

    @staticmethod
    def _make_ins_str_change(ins_changes_pos, ins_changes_alts) -> List[str]:
        assert len(ins_changes_pos) == len(ins_changes_alts)
        ins = np.column_stack((
            ins_changes_pos.astype(str), 
            np.repeat("_|", len(ins_changes_pos)), 
            ins_changes_alts))
        ins_str = [''.join(ins[x,:]) for x in np.arange(ins.shape[0])]
        return ins_str

    @staticmethod
    def _read_insertions_file(file_path) -> pd.DataFrame:
        """
        Return a table with coordinates (0-based) on the columns and sequence names as the row index. Every cell contains either an empty string or the insertion string.
        """
        ins_coordinates = None
        records = []
        with open(file_path, "r") as ins_file:
            # read_header
            header = ins_file.readline()
            ins_coordinates = [x.strip() for x in header.split(",")[1:]]    # skip first column ("strain")
            ins_coordinates = [x.rsplit(" ", maxsplit=1)[1] for x in ins_coordinates]
            ins_coordinates = np.array(ins_coordinates).astype(int)
            # read rest of file
            for line in ins_file:
                split_columns = line.rstrip("\n").rsplit(",", maxsplit=len(ins_coordinates))   # such weird logic is necessary because "strain" name can contain ","
                records.append(split_columns)
        # make a DataFrame
        column_map = {0: "index"}          
        column_map.update({x+1: ins_coordinates[x] for x in range(len(ins_coordinates))})
        df = pd.DataFrame.from_records(records).rename(columns=column_map).set_index("index")
        return df


def df2fasta(df: pd.DataFrame, out_file_path):
    """
    df: a DataFrame with one row per sequence, where the sequence is stored in a column named "sequence" and 
    the index is used as sequence name. 
    """
    def get_fasta_objects(file_input):
        for idx, row in file_input.iterrows():
            yield miniFasta.fasta_object(head=idx, body=row.sequence)

    miniFasta.write(fasta_pairs=list(get_fasta_objects(df)), file_path=out_file_path)


def extract_annotated_seq_after_augur(md, reference_file, aligned_file, insertions_file, annotation_list, verbose=False):
    """
    md: a DataFrame with a unique index (identifier of the sequence and a column 'sequence')
    reference_file: path to the reference fasta file to read
    aligned_file: path to the aligned fasta file to read
    insertions_file: path to the insertions file to read
    annoation_list: a list with start and stop (excluded) pos of the annotation to extract from the target sequences

    Returns two dictionaries: the first containing the changes of all the sequences, the second dictionary contains the annotated sequences.
    """

    # read refseq
    refseq = next(miniFasta.read(reference_file)).body
    if verbose:
        print("Refseq len", len(refseq))

    # call variants
    all_changes = dict()
    all_annotated_seq = dict()
    caller = VariantCaller(refseq=refseq, 
                        insertions_file_path=insertions_file)
    for fo in tqdm(miniFasta.read(aligned_file), total=md.shape[0]):
        target_name, target = fo.head.lstrip(">"), fo.body
        if not target_name in md.index:
            continue
        changes: AlignmentChanges = caller.call_variants(name=target_name, 
                                                        aligned=target, 
                                                        aggregate_close_changes=True, 
                                                        one_based=False,
                                                        strip_edge_insertions=True)   # strip edge insertions is more convenient in bulk afterwards   
        all_changes[target_name] = changes
        ann_le = annotation_list[0]
        ann_ri = annotation_list[1]
        
        try:

            ltarget = list(target)    
                    
            # # add substitutions (does not change the list length) -- NOT REQUIRED: TARGET ALREADY CONTAINS SUBSTITUTIONS
            # for c_pos, c in zip(changes.subs_pos, changes.subs):
            #     alt = c[c.index("|")+1:]
            #     ltarget[c_pos:c_pos + len(alt)] = list(alt)

            # add insertions (each insertion introduces one more list element -> annotation may be impacted)
            for c_pos, c in zip(changes.ins_pos, changes.ins):
                alt = c[c.index("|")+1:]
                ltarget.insert(c_pos, alt)      # insertion causes the list to acquire 1 element (the elem can contain a long string of inserted chars)
                if c_pos < ann_ri:  # move annotation end in rightward direction
                    ann_ri += 1
                if c_pos < ann_le:  # move annotation start in rightward direction
                    ann_le += 1

            
            # print("AFTER INSERTIONS: List length:", len(ltarget), "string length:", len("".join(ltarget)))
            # print("AFTER INSERTIONS: CDS List length:", len(ltarget[ann_le:ann_ri]), "CDS string length:", len("".join(ltarget[ann_le:ann_ri])))
                                        
            annotated_seq = ltarget[ann_le:ann_ri] 
            annotated_seq = "".join(annotated_seq)
        
            # if len(annotated_seq) % 3 != 0:
            #     raise AssertionError("Become invalid after subs and insertions and transforming list to string")
            
            annotated_seq = annotated_seq.replace("-", "")

            # if len(annotated_seq) % 3 != 0:
            #     raise AssertionError("Become invalid after deletions??")
            
            all_annotated_seq[target_name] = annotated_seq

            

        except AssertionError:
            print(target_name, "has CDS not multiple of 3")
            print("ANNOTATION RANGE", ann_le, ann_ri)
            print("insertions:", changes.ins)
            print("deletions:", changes.dels)
            print("substitutions", changes.subs)
            raise 

    return all_changes, all_annotated_seq


if __name__ == '__main__':
    # read refseq
    for fo in miniFasta.read("reference_h5nx_ha.fasta"):
        refseq = fo.body
    print("Refseq len", len(refseq))

    # read EPI_ISL_219874_A/duck/Taiwan/A3400/2015
    target_name, target = None, None
    for fo in miniFasta.read("aligned_data.fasta"):
        target_name = fo.head.lstrip(">")
        if target_name != "EPI_ISL_219874_A/duck/Taiwan/A3400/2015":
            continue
        else:
            target = fo.body 
            break
    print(f"one_sequence: {target_name}")
    print("one sequence len", len(target))

    # call variants
    # substitutions and deletions are called from the aligned sequence
    # insertions come from a dedicated file (...insertions.csv) returned by augur align. Variant Caller just reads the file and parse them
    # in VariantCaller insertions file_path and name can be omitted if you are not interested in insertions
    caller = VariantCaller(refseq=refseq, 
                           insertions_file_path="data/HA_MP_2.3.4.4bVSother_2000on/3.2_aligned_cleaned_fasta/aligned_data_HA_cleaned.fasta.insertions.csv")
    changes: AlignmentChanges = caller.call_variants(name=target_name, 
                                                     aligned=target, 
                                                     aggregate_close_changes=True, 
                                                     one_based=True,
                                                     strip_edge_insertions=False)   # is more convenient to strip edge insertions in bulk afterwards   

    # print variants
    print(changes.dels, type(changes.dels))
    print(changes.dels_pos, type(changes.dels_pos))
    print(changes.subs, type(changes.subs))
    print(changes.subs_pos, type(changes.subs_pos))
    print(changes.ins, type(changes.ins))
    print(changes.ins_pos, type(changes.ins_pos))