from typing import TypedDict, Any
import pandas as pd
import numpy as np
from collections import Counter
import numpy as np
import pandas as pd
pd.set_option('future.no_silent_downcasting', True)
from cleverdict import CleverDict
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams.update(mpl.rcParamsDefault)     # reset matplotlib style
plt.style.use('fast')
import requests
from sktime.detection.stray import STRAY
from sklearn.preprocessing import MinMaxScaler
from utils.genetic_code import GeneticCode
from utils.data_cleaning import DataCleaning, DateRange
from job_orchestra import Step
import shlex
import subprocess
from utils.data_exploration import Explorer
from math import floor, ceil
from utils.variant_calling import extract_annotated_seq_after_augur
from tqdm import tqdm
import miniFasta
import os

class BaseResultType(TypedDict, total=False):   
    # the following are sugegsted keys. Their presence is not mandatory.
    # data 
    output: Any
    feature_names: list | np.ndarray
    data_feature_names: list | np.ndarray
    metadata_feature_names: list | np.ndarray
    # clustering
    cluster_centroids: np.ndarray
    cluster_names: list | np.ndarray
    cluster_size: list | np.ndarray
    cluster_score: float
    # ml and statistics
    model: dict
    # statistics
    mu: int | list | np.ndarray
    sigma: int | list | np.ndarray 

class ResultType(BaseResultType):
    pass    # subclassing BaseResultType allow declaring a dictionary of type ResultType with additional keys not listed in BaseResultType without triggering a type warning

# DATA
class InputH1N1:

    class AttachBMC_GenomicsCluster(Step):
        def run(self, data: ResultType):
            bmc_data = pd.read_parquet("data/h1n1/H1N1_bmc_cluster_label").rename(columns={'cluster':'bmc_cluster_label'})[['bmc_cluster_label']]              
            output = pd.merge(data.output, bmc_data, left_index=True, right_index=True, how="left")
            return {
                'output': output,
                'data_feature_names':data.data_feature_names,
                'metadata_feature_names': data.metadata_feature_names
            }

class InputH5N1:

    class AttachBMC_GenomicsCluster(Step):
        def run(self, data: ResultType):
            repo_dom_bird = pd.read_parquet("data/h5n1/h5n1_dom_bird_zenodo_repo.parquet")
            repo_dom_bird = repo_dom_bird.rename(columns={'cluster':'bmc_cluster_label'})[['bmc_cluster_label']]
            repo_dom_bird = repo_dom_bird.replace({np.nan: pd.NA})
            repo_dom_bird['bmc_cluster_label'] = repo_dom_bird['bmc_cluster_label'].astype(str).apply(lambda x: "dom_" +x)
            print(f"Collected {repo_dom_bird.shape[0]} cluster labels for H5N1 domestic birds")
            
            repo_wild_bird = pd.read_parquet("data/h5n1/h5n1_wild_bird_zenodo_repo.parquet")
            repo_wild_bird = repo_wild_bird.rename(columns={'cluster':'bmc_cluster_label'})[['bmc_cluster_label']]
            repo_wild_bird = repo_wild_bird.replace({np.nan: pd.NA})
            repo_wild_bird['bmc_cluster_label'] = repo_wild_bird['bmc_cluster_label'].astype(str).apply(lambda x: "wild_" +x)
            print(f"Collected {repo_wild_bird.shape[0]} cluster labels for H5N1 wild birds")

            assert len(set(repo_wild_bird.index) & set(repo_dom_bird.index)) == 0, "Multiple cluster labels assigned for the same sequence in wild + domestic bird"
            bmc_cluster_label = pd.concat([repo_dom_bird, repo_wild_bird])['bmc_cluster_label']
            output = pd.merge(data.output, bmc_cluster_label, left_index=True, right_index=True, how="left")
            return {
                'output': output,
                'data_feature_names':data.data_feature_names,
                'metadata_feature_names': data.metadata_feature_names + ['bmc_cluster_label']
            }

class Input:
    class ReadMetadata(Step):
        def __init__(self, *metadata_file_paths, verbose=False, **kwargs):
            super().__init__(**kwargs)
            self.verbose = verbose
            self.metadata_file_paths = metadata_file_paths
        def run(self):
            meta = Explorer.read_meta(*self.metadata_file_paths, verbose=self.verbose)
            return {
                'output': meta,
                'metadata_feature_names': list(meta.columns)
            }

    class ReadData(Step):
        def __init__(self, *data_file_paths, verbose=False, **kwargs):
            super().__init__(**kwargs)
            self.verbose = verbose
            self.data_file_paths = data_file_paths
        def run(self):
            data = Explorer.read_data(*self.data_file_paths, verbose=self.verbose)
            return {
                'output': data,
                'data_feature_names': list(data.columns),
                'metadata_feature_names': []
            }
    
    class FilterSeqLen(Step):
        def __init__(self, len_range: tuple, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.len_range = len_range
        def run(self, data_input: ResultType):
            return {
                'output': data_input.output[data_input.output.seq_len.between(*self.len_range)].copy(),
                'data_feature_names': data_input['data_feature_names'],
                'metadata_feature_names': data_input['metadata_feature_names'] if 'metadata_feature_names' in data_input.keys() else ['']
            }

    class MetadataCleaning(Step):
        def __init__(self, allow_incomplete_collection_dates, allow_null_submission_dates, *args, additional_mappings=None, **kwargs):
            super().__init__(*args, **kwargs)
            self.allow_incomplete_collection_dates = allow_incomplete_collection_dates
            self.allow_null_submission_dates = allow_null_submission_dates
            self.additional_mappings=additional_mappings
        def run(self, metadata):
            meta = metadata.output.copy()
            # collection date (may be incomplete or null)
            if self.allow_incomplete_collection_dates:
                meta.Collection_Date = pd.to_datetime(meta.Collection_Date, format='ISO8601')
            else:
                meta.Collection_Date =  pd.to_datetime(meta.Collection_Date, format="%Y-%m-%d", errors='coerce')
                meta = meta[~pd.isna(meta.Collection_Date)]
            # submission date is either completely specified or null
            meta.Submission_Date =  pd.to_datetime(meta.Submission_Date, format="%Y-%m-%d")
            if not self.allow_null_submission_dates:
                meta = meta[~pd.isna(meta.Submission_Date)]
            # host and host_type
            meta = DataCleaning.merge_fix_host_species(meta)
            meta = DataCleaning.assign_host_type(meta, additional_mappings=self.additional_mappings, raise_on_error=True)
            meta = meta[meta.Host_Type != "other"].copy()
            # clade
            meta.Clade = meta.Clade.replace({"unassigned": np.nan})
            # flu-season
            meta = DataCleaning.assign_flu_season(meta)
            # genotype
            if 'Genotype' in meta.columns:
                meta['Genotype'] = meta['Genotype'].str.replace(" (<i style='font-size:11px'>GenoFLU</i>)", "", regex=False)
                meta['Genotype'] = meta['Genotype'].str.replace("Not assigned", "Notassigned", regex=False)
                meta['Genotype'] = meta['Genotype'].str.replace("\n", "", regex=False)
                meta['Genotype'] = meta['Genotype'].replace({"Notassigned": np.nan})
                meta['Genotype'] = meta['Genotype'].fillna(np.nan)
            return {
                'output': meta, 
                'metadata_feature_names': list(meta.columns),
                'data_feature_names': []
            }
        
    class MetaDataMerge(Step):
        def run(self, metadata_input, data_input):
            meta = metadata_input.output
            data = data_input.output   
            # merge
            md = pd.merge(meta, data, left_index=True, right_index=True, how="inner")
            # remove duplicate indexes
            md = DataCleaning.remove_duplicate_index(md, verbose=False)
            return {
                'output': md,
                'data_feature_names': data_input['data_feature_names'],
                'metadata_feature_names': metadata_input['metadata_feature_names']
            }
        
    class FilterSegm(Step):
        def __init__(self, segm_name, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.segm_name = segm_name
        def run(self, input_data: ResultType):
            return {
                'output': input_data.output[input_data.output.segm == self.segm_name],
                'metadata_feature_names': input_data.metadata_feature_names,
                'data_feature_names': input_data.data_feature_names
            }
        
    class SegmCleaning(Step):
        def __init__(self, max_n_ratio: float, mode_length_perc_deviation: int, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.max_n_ratio = max_n_ratio
            self.mode_length_perc_deviation = mode_length_perc_deviation

        def run(self, md_input):
            md = md_input.output
            assert {'seq_len', 'n_ratio', 'sequence'}.issubset(md.columns), "Input data frame misses one or more required columns: 'segm', 'seq_len', 'sequence'"
            # remove seqeunces with too many Ns
            md = md[md.n_ratio <= self.max_n_ratio]
            # remove sequences with wrong length
            mode_length = md.seq_len.mode().item()
            three_perc_deviation = mode_length * self.mode_length_perc_deviation / 100
            allowed_length_range = (mode_length - floor(three_perc_deviation), mode_length + ceil(three_perc_deviation))
            md = md[md.seq_len.between(*allowed_length_range)]
            return {
                'output': md, 
                'data_feature_names':  md_input['data_feature_names'],
                'metadata_feature_names': md_input['metadata_feature_names']
            }
        
    class SegmAlignment(Step):
        """
        As bad sequences can influence the global output of augur, it's best to filter out bad sequences (high N% or weird lengths before doing the alignment)
        """
        def __init__(self, target_file_path, reference_file_path: str, nthreads: str|int = 'auto', preview=False,
                    *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.target_file = target_file_path
            self.reference_file = reference_file_path
            self.nthreads = nthreads
            # paths to store future alignments
            self.aligned_file = self.target_file.replace(".fasta", "_aligned.fasta")
            self.insertions_file = f"{self.aligned_file}.insertions.csv"
            self.preview = preview

        def run(self, md_input):
            md = md_input.output
            assert 'sequence' in list(md.columns), "Input data frame misses one or more required columns: 'segm', 'seq_len', 'sequence'"
            # align
            df2fasta(md, self.target_file)
            curr_interpreter = sys.executable
            command = f"{curr_interpreter} -m augur align --reference-sequence {self.reference_file} --sequences {self.target_file} --output {self.aligned_file} --nthreads {self.nthreads}"
            if self.preview:
                print("Will execute:", command)
            else:
                print("Running:", command)
                # check augur is available
                if importlib.util.find_spec('augur') is None:
                    raise ModuleNotFoundError("augur is not available in the current interpreter. Install augur.")
                try:
                    completed_process = subprocess.run(shlex.split(command), check=True, capture_output=True, text=True)   # if completed_process.returncode != 0 raise subprocess.CalledProcessError
                except subprocess.CalledProcessError as cpe:
                    print(f"Error while invoking augur with command {command}")
                    raise cpe
                # alternative nethod to call augur
                    # completed_process = subprocess.Popen(
                    #     shlex.split(command), stdout=subprocess.PIPE
                    # )
                    # output, error = completed_process.communicate()
            return {
                'aligned_file': self.aligned_file,
                'insertions_file': self.insertions_file
            }

    class SegmCDS(Step):
        """
        As bad sequences can influence the global output of augur, it's best to filter out bad sequences (high N% or weird lengths before doing the alignment)
        """
        def __init__(self, reference_file_path: str, annotation_range_1_based: tuple, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.reference_file = reference_file_path
            a = annotation_range_1_based
            self.annotation_range = (a[0]-1, a[1]) # convert to 0-based. The 0-based conversion of the right edge is -1 (because we change the numebring convention) + 1 (to transform the range from [] to [).

        def run(self, md_input, segm_alignment):
            md = md_input.output.copy()
            aligned_file = segm_alignment.aligned_file
            insertions_file = segm_alignment.insertions_file
            changes_dict, annotated_seq_dict = extract_annotated_seq_after_augur(md, self.reference_file, aligned_file, insertions_file, self.annotation_range, verbose=False)
            md["CDS"] = annotated_seq_dict
            return {
                'output': md, 
                'data_feature_names':  md_input['data_feature_names'] + ['CDS'],
                'metadata_feature_names': md_input['metadata_feature_names']
            }
    
    class FilterValidCDS(Step):
        def run(self, input_data: ResultType):
            md = input_data.output
            data_size = md.shape[0]
            # consider only CDS divisible by 3
            md = md.astype({'CDS': "str"})
            md = md[md.CDS.apply(len) % 3 == 0]
            filtered_data_size = md.shape[0]
            print(f"Extracted valid CDS sequence for {filtered_data_size} sequences ({(filtered_data_size / data_size * 100):.2f} % of aligned sequences)")
            return {
                'output': md, 
                'data_feature_names':  input_data['data_feature_names'],
                'metadata_feature_names': input_data['metadata_feature_names']
            }
        
    class AssertCDSLength(Step):
        def __init__(self, reference_length: int, perc_deviation: int, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.reference_length = reference_length
            self.perc_deviation = perc_deviation
        def run(self, data_input: ResultType):
            deviation = self.reference_length * self.perc_deviation / 100
            min_length, max_length = self.reference_length - floor(deviation), self.reference_length + ceil(deviation)
            if not data_input.output.CDS.str.len().between(min_length, max_length).all():
                invalid_cds = data_input.output[~data_input.output.CDS.str.len().between(min_length, max_length)]
                if "HA" in self.name_alias:
                    print(f"Skipping {invalid_cds.shape[0]} sequence having CDS length != {self.perc_deviation}% of the reference CDS length: {invalid_cds.index.tolist()}")
                    return {
                        'output': data_input.output[data_input.output.CDS.str.len().between(min_length, max_length)],
                        'data_feature_names': data_input.data_feature_names,
                        'metadata_feature_names': data_input.metadata_feature_names
                    }
                else:
                    raise AssertionError(f"{invalid_cds.shape[0]} sequences have CDS length != {self.perc_deviation}% of the reference CDS length: {invalid_cds.index.tolist()}")
            return data_input
        
        # EXAMPLE USAGE
        # S_input_meta = ReadMetadata(files)
        # S_meta_cleaning = MetadataCleaning(depends_on=S_input_meta, allow_incomplete_collection_dates=False, allow_null_submission_dates=True)
        # # select partition based on metadata (host, location, time)
        # # consider filtering by segment in ReadData, to avoid huge joins in MetaDataMerge
        # S_input_data = ReadData(other_files)
        # S_metadata_merge = MetaDataMerge(depends_on=[S_meta_cleaning, S_input_data])
        # S_segm_data = FilterSegm(segm_name, depends_on=S_metadata_merge)
        # S_segm_cleaning = SegmCleaning(0.02, 7, depends_on=S_segm_data)
        # S_seq_alignment = SegmAlignment(target_file_path, reference_file_path, nthreads = 'auto', depends_on=S_segm_cleaning)
        # S_seq_cds_raw = SegmCDS(reference_file_path, annotation_range_1_based, depends_on=S_seq_alignment)
        # S_seq_cds = FilterValidCDS(depends_on=S_seq_cds_raw)
            
    class ReadAlignment(Step):
        def __init__(self, reference_file_path, alignment_file_path, indels_file_path, annoation_range, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.ref = reference_file_path
            self.afp = alignment_file_path
            self.infp = indels_file_path
            self.ann_range = annoation_range

        def run(self, data: ResultType):
            md = data.output
            changes_dict, annotated_seq_dict = extract_annotated_seq_after_augur(md, self.ref, self.afp, self.infp, self.ann_range, verbose=False)
            md["CDS"] = annotated_seq_dict
            # consider only CDS divisible by 3
            md = md[md.CDS.apply(len) % 3 == 0]
            return {
                'output': md, 
                'data_feature_names':  data['data_feature_names'] + ['CDS'],
                'metadata_feature_names': data['metadata_feature_names']
            }
    
    class ReindexGISAIDFastaAlignment(Step):
        def __init__(self, alignment_file_path, indels_file_path, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.afp = alignment_file_path
            self.infp = indels_file_path
            afp_dir, _, afp_name = self.afp.rpartition(os.path.sep)
            self.output_afp = afp_dir + os.path.sep + afp_name.replace(".fasta", "_") + "rewritten.fasta"
            infp_dir, _, infp_name = self.infp.rpartition(os.path.sep)
            self.output_infp = infp_dir + os.path.sep + infp_name.replace(".fasta.insertions.csv", "_") + "rewritten.fasta.insertions.csv"
        def run(self):
            # rewrite the input alignements with the correct index
            output_fo = []
            ids = set() # avoid duplicates
            for si, fo in tqdm(enumerate(miniFasta.read(self.afp))):
                foheader = fo.getHead()[1:]     # [1:] removes ">"

                # identify sequence
                try:
                    _,isolate_id,virus_name,type_subtype,lin,clade,segm = foheader.split("|")
                except Exception as e:
                    print(f"Fasta header unpack unusccesful on sequence n. {si} (0-based) of file {self.afp}. Fasta header is: {foheader}.\nSequence ignored.")
                id = (isolate_id + "_" + virus_name).strip()
                if not id in ids:
                    ids.add(id)
                    sequence = fo.getSeq()
                    output_fo.append(miniFasta.fasta_object(id,sequence))
            miniFasta.write(output_fo, self.output_afp)
            # rewrite the isnertions file with the correct index
            ids = set() # avoid duplicates
            with open(self.infp, "r") as inp_ins:
                header = inp_ins.readline()  # skip header
                with open(self.output_infp, "w") as out_ins:
                    out_ins.write(header)
                    for l in inp_ins.readlines():
                        # rewrite with correct index
                        old_id, rest = l.split(",", maxsplit=1)
                        _,isolate_id,virus_name,type_subtype,lin,clade,segm = old_id.split("|")
                        new_id = (isolate_id + "_" + virus_name).strip()
                        if not new_id in ids:
                            ids.add(new_id)
                            out_ins.write(f"{new_id},{rest}")
            return {'aligned_file': self.output_afp, 'insertions_file': self.output_infp}

    class Sequence2Fasta(Step):
        """
        Expects a dataset with at least a column 'sequence'.
        """
        def __init__(self, output_path, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output_path = output_path
        def run(self, input_data: ResultType):
            data = input_data.output
            miniFasta.write([miniFasta.fasta_object(idx, row.sequence) for idx,row in data.iterrows()], self.output_path)               

    class RSCU(Step):
        def run(self, input_data: ResultType):
            data = input_data.output
            df_rscu = data.apply(lambda x: pd.Series(GeneticCode.rscu([x.CDS])), axis=1)
            new_data = pd.merge(data, df_rscu, left_index=True, right_index=True, validate="1:1")
            return {
                'output': new_data,
                'data_feature_names': input_data.data_feature_names + GeneticCode.valid_codons,
                'metadata_feature_names': input_data.metadata_feature_names
            }
        
    class ConcatShuffleDatasets(Step):
        def __init__(self, ignore_index: bool, *args, random_seed=None, **kwargs):
            super().__init__(*args, **kwargs)
            self.ignore_index = ignore_index
            self.random_seed = random_seed
        def run(self, *data_inputs: list[ResultType]):
            metadata_feature_names = data_inputs[0].metadata_feature_names
            data_feature_names = data_inputs[0].data_feature_names
            output_df_colums = data_inputs[0].output.columns.tolist()
            for di in data_inputs:
                assert di.metadata_feature_names == metadata_feature_names, "Datasets have different values of metadata_feature_names"
                assert di.data_feature_names == data_feature_names, "Datasets have different values of data_feature_names"
                assert di.output.columns.tolist() == output_df_colums, "Datasets have different columns of different order of the columns"
            out = pd.concat([d.output for d in data_inputs], ignore_index=self.ignore_index).sample(frac=1, random_state=self.random_seed)
            return {
                'output': out,
                'data_feature_names': data_feature_names,
                'metadata_feature_names': metadata_feature_names
            }

    class DropMetadata(Step):
        def run(self, data_input: ResultType):
            return {
                'output': data_input.output[data_input.data_feature_names],
                'metadata_feature_names': [],
                'data_feature_names': data_input.data_feature_names # ['segm', 'seq_len', 'n_ratio', 'sequence', 'CDS']
            }

    class MoveDataFeaturesToMetadata(Step):
        def run(self, data_input: ResultType):
            return {
                'output': data_input.output, 'metadata_feature_names': data_input.metadata_feature_names + data_input.data_feature_names, 'data_feature_names': []
            }
        
    class NPartitions(Step):
        def __init__(self, n, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.n = n
        def run(self, input_data: ResultType):
            data = input_data.output
            data_feature_names = input_data.data_feature_names
            metadata_feature_names = input_data.metadata_feature_names
            available_data_size = input_data.output.shape[0]
            size_partitions = available_data_size // self.n
            partitions_start_indices = np.array([i*size_partitions for i in range(self.n)])
            partitions_end_indices = partitions_start_indices + size_partitions
            partitions_end_indices[-1] = available_data_size
            partitions = [data.iloc[s:e] for s,e in zip(partitions_start_indices, partitions_end_indices)]
            return [CleverDict({'output': v, 'data_feature_names': data_feature_names, 'metadata_feature_names': metadata_feature_names}) for v in partitions]
         

class AttachWeek(Step):
    def run(self, input_data: ResultType):
        output = input_data.output.copy()
        output['Week'] = output['Collection_Date'].dt.strftime("%G-%V")
        return {
            'output': output,
            'metadata_feature_names': input_data.metadata_feature_names + ['Week'],
            'data_feature_names': input_data.data_feature_names
        }
    
class AttachCountryOrState(Step):
    def __init__(self, *args, exact_replacements: dict = {}, substring_replacements: dict = {}, **kwargs):
        super().__init__(*args, **kwargs)
        self.exact_replacements = exact_replacements
        self.substring_replacements = substring_replacements
    def run(self, input_data: ResultType):
        output = input_data.output.copy()
        country_or_state = (
            output.Location
            .str.replace("nan", "")
            .str.replace(" / EPI_ISL_19541397 refers - with a high likelihood - to the same patients from submissions EPI_ISL_19512045, EPI_ISL_19512046, EPI", "")
            .str.replace("North America / ", "")
            .str.replace("United States / ", "USA, ")
            .str.split(" / ").str[0]
            .replace("North America", np.nan)
            .str.replace("United States", "USA")
            .str.replace("Unknown", "")
            .str.rstrip(" /")
            .str.replace("Conneticut", "Connecticut")
            .replace("USA,", "USA")
            .str.replace("Deleware", "Delaware")
            .str.replace("Massachussetts", "Massachusetts")
            .str.replace("USA, Puerto Rico", "Puerto Rico")
            .str.replace("CaliforniSWL7/2009", "California")
            .replace("", np.nan)
        )
        for k,v in self.exact_replacements:
            country_or_state = country_or_state.replace(k,v)
        for k,v in self.substring_replacements:
            country_or_state = country_or_state.str.replace(k,v)
        output['country_or_state'] = country_or_state
        return {'output': output, 'metadata_feature_names': input_data.metadata_feature_names + ['country_or_state'], 'data_feature_names': input_data.data_feature_names}

class WindowSelector(Step):
    def __init__(self,  date_range, start, end, *args, inclusive="both", **kwargs):
        super().__init__(*args, **kwargs)
        self.date_range = date_range
        self.start = start
        self.end = end
        self.inclusive = inclusive

    def run(self, global_data: ResultType):
        return {
            'start':self.start, 'end':self.end,
            'output': global_data.output[global_data.output.Collection_Date.between(self.date_range.start, self.date_range.end, self.inclusive)].copy(),
            'data_feature_names': global_data.data_feature_names,
            'metadata_feature_names': global_data.metadata_feature_names
        }

class IndexBasedFilter(Step):
    def run(self, input_data_to_filter: ResultType, input_data_filtered: ResultType):
        return_val = input_data_to_filter
        return_val.output = return_val.output.loc[input_data_filtered.output.index].copy()
        return return_val

class AbsoluteSortByCollectionDateAndIndex(Step):
    def run(self, input_data: ResultType):
        data = input_data.output.copy()
        data['Sort_Key'] = data.Collection_Date.dt.strftime("%Y-%m-%d") + "_" + data.index.to_series().astype(str)
        n_duplicates = data.shape[0] - len(set(data['Sort_Key'].tolist()))
        assert n_duplicates == 0, f"Absolute sort of dataset impossible: {n_duplicates} duplicate value(s) in column 'Sort_Key' found."
        data = data.sort_values(by=['Sort_Key'])
        return {
            'output': data,
            'data_feature_names': input_data.data_feature_names,
            'metadata_feature_names': input_data.metadata_feature_names + ['Sort_Key']
        }

class MovingWindowFixedSize(Step):
    """
    Train and test in this context are interpreted as follows:
    Fixed size partitions may need data points that are not included in the given input_data (i.e., data points preceeding the first one). Such partitions are ignored (the output range_names etc. do not contain such partitions).
    """
    def __init__(self, evaluation_date_range: DateRange, train_size=99, test_size=1, window_shift=1, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.evaluation_date_range = evaluation_date_range
        self.train_size = train_size
        self.test_size = test_size
        self.window_shift = window_shift
    def run(self, input_data: ResultType) -> ResultType:
        data = input_data.output
        data_len = data.shape[0]
        assert not data.empty, "Input DataFrame is empty. Cannot extract any partition."

        evaluation_date_start = self.evaluation_date_range.start
        evaluation_date_end = self.evaluation_date_range.end
        try:
            evaluation_data_start_index = data[data.Collection_Date >= evaluation_date_start].iloc[0].name
        except IndexError:
            raise IndexError("Can't produce any partition because the evaluation_date_range correponds to an empty data frame for the given input_data")
        try:
            evaluation_data_end_index = data[data.Collection_Date >= evaluation_date_end].iloc[0].name
        except IndexError:
            evaluation_data_end_index = data.iloc[-1].name

        evaluation_data_positional_start, evaluation_data_positional_end = data.index.get_indexer([evaluation_data_start_index, evaluation_data_end_index])

        # return [data.loc[evaluation_data_start_index].Collection_Date.strftime("%G-%V-%u"), 
        #         data.loc[evaluation_data_end_index].Collection_Date.strftime("%G-%V-%u"),
        #         evaluation_data_positional_start, evaluation_data_positional_end]

        test_ranges = np.arange(evaluation_data_positional_start, evaluation_data_positional_end, self.window_shift)
        test_ranges = list(zip(test_ranges.tolist(), (test_ranges+self.test_size).tolist()))
        
        training_ranges = np.arange(evaluation_data_positional_start-self.train_size, evaluation_data_positional_end-self.train_size, self.window_shift)
        training_ranges = list(zip(training_ranges.tolist(), (training_ranges+self.train_size).tolist()))

        train_plus_test_ranges = [(tr[0],te[1]) for tr,te in zip(training_ranges, test_ranges)]

        # remove invalid partitions
        valid_partitions_k = [k for k,(li,re) in enumerate(train_plus_test_ranges) if li >= 0 and re < data_len]
        if len(valid_partitions_k) > 0:
            valid_partitions_k_first = valid_partitions_k[0]
            valid_partitions_k_last = valid_partitions_k[-1] + 1
            
            train_plus_test_ranges = train_plus_test_ranges[valid_partitions_k_first:valid_partitions_k_last]
            training_ranges = training_ranges[valid_partitions_k_first:valid_partitions_k_last]
            test_ranges = test_ranges[valid_partitions_k_first:valid_partitions_k_last]
        else:
            train_plus_test_ranges = []
            training_ranges = []
            test_ranges = []

        train_range_iterator = iter(training_ranges)
        test_range_iterator = iter(test_ranges)
        train_plus_test_iterator = iter(train_plus_test_ranges)
        
            # return [data.loc[evaluation_data_start_index].Collection_Date.strftime("%G-%V-%u"), 
            #         data.loc[evaluation_data_end_index].Collection_Date.strftime("%G-%V-%u"), "\n", 
            #         evaluation_data_positional_start, evaluation_data_positional_end, "\n",
            #         list(zip(training_ranges, test_ranges))[:5], list(zip(training_ranges, test_ranges))[-5:], "\n",
            #         train_plus_test_ranges[:5], train_plus_test_ranges[-5:]]
            
        def next_test_partition():
            return next(test_range_iterator)
        
        def next_training_partition():
            return next(train_range_iterator)

        def next_train_plus_test_partition():
            return next(train_plus_test_iterator)
        
        return {
            '.next_training_partition': next_training_partition,
            '.next_test_partition': next_test_partition,
            '.next_train_plus_test_partition': next_train_plus_test_partition,
            # USE ISO 8601 for week numbering (%V) and year numbering (%G)
            'train_range_names': dict(enumerate(training_ranges)),
            'test_range_names': dict(enumerate(test_ranges)),
            'train_plus_test_range_names': dict(enumerate(train_plus_test_ranges)),
            'data': data,
            'data_feature_names': input_data['data_feature_names'],
            'metadata_feature_names': input_data['metadata_feature_names']
        }
    
class NumericIndexBasedPartition(Step):
    def __init__(self, partitioner_function_name: str, range_start: str, range_end: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.partitioner_function = partitioner_function_name
        self.range_start = range_start
        self.range_end = range_end

    def run(self, partitioner: ResultType) -> ResultType:
        numeric_index_range = partitioner[self.partitioner_function]() # function call
        
        md = partitioner.data
        md_out = md.iloc[numeric_index_range[0]:numeric_index_range[1],:]
        return {
            'output': md_out,
            'data_feature_names': partitioner['data_feature_names'],
            'metadata_feature_names': partitioner['metadata_feature_names'],
            'range_start': self.range_start,
            'range_end': self.range_end
        }
    
class AbsoluteSortByCollectionDateAndEpisetID(Step):
    def run(self, input_data: ResultType):
        data = input_data.output.copy()
        metadata_feature_names = input_data.metadata_feature_names
        # sort by Collection Date + Isolate Num (EPI_ISL_#####)
        if pd.api.types.is_string_dtype(data.Isolate_Id):
            data['Isolate_Num'] = data['Isolate_Id'].str.lstrip("EPI_ISL_").astype(int)
            data = data.sort_values(by=['Collection_Date', 'Isolate_Num'])
            metadata_feature_names += ['Isolate_Num']
        else:
            assert pd.api.types.is_numeric_dtype(data.Isolate_Id)
            data = data.sort_values(by=['Collection_Date', 'Isolate_Id'])
        data['Sort_Key'] = np.arange(data.shape[0])
        metadata_feature_names += ['Sort_Key']
        return {
            'output': data,
            'metadata_feature_names': metadata_feature_names,
            'data_feature_names': input_data.data_feature_names
        }

# input tranformation
class DataTransform:
                
    class LogBySynCount(Step):
        """
        Transform the RSCUs into log(RSCUs) where the logarithm has a base ranging from 2 to 6, according to the number of synonym codons available for the specific amino acid. 
        """    
        def run(self, data_input: ResultType):
            output = data_input['output']
            data_feature_names = data_input['data_feature_names']
            metadata_feature_names = data_input['metadata_feature_names']
            data = output[data_feature_names]
            metadata = output[metadata_feature_names]
            output = np.emath.logn([GeneticCode.syn_count[c] for c in data_feature_names], data.values)
            output = pd.DataFrame(output, columns=data_feature_names, index=data.index)
            output = pd.merge(output, metadata, left_index=True, right_index=True, how="inner", validate="1:1")
            return {
                'output': output, 
                'data_feature_names': data_feature_names, 
                'metadata_feature_names': data_input['metadata_feature_names']
            }
        
    class PlainAndLogDinucleotide(Step):
        """
        Computes a list of 77 features. The first 59 are the RSCUs of the codons. The next 16 are the dinucleotide frequency ratios.
        """
        new_data_feature_names = GeneticCode.valid_codons + list(GeneticCode.dinucleotides)

        def run(self, input: ResultType):
            md = input.output
            data_feature_names = input['data_feature_names']
            meta_feature_names = input['metadata_feature_names']
            assert 'CDS' in meta_feature_names, "The feature named 'CDS' is not available in the input. Can't compute dinucleotides without the sequence."
            # transform input
            rscu = md[data_feature_names]
            dinucleotides = md[['CDS']].apply(GeneticCode.dinucleotides_count, axis=1, result_type='expand').map(np.log2)
            data_out = pd.merge(left=rscu, right=dinucleotides, left_index=True, right_index=True, how="inner", validate="1:1")
            output = pd.merge(data_out, md[meta_feature_names], left_index=True, right_index=True, how="inner", validate="1:1")
            return {
                'output': output,
                'data_feature_names': DataTransform.PlainAndLogDinucleotide.new_data_feature_names,
                'metadata_feature_names': meta_feature_names
            }


class TrainTestWarningsCollector(Step):
    def run(self, *warnings: ResultType):
        return {
            'train_ranges': [f"{k['train_start']}_{k['train_end']}" for k in warnings],
            'test_ranges': [f"{k['test_start']}_{k['test_end']}" for k in warnings],
            
            'train_total': [k['train_total'] for k in warnings],
            'train_na': [k['train_na'] for k in warnings],
            'train_na_perc': [k['train_na_perc'] for k in warnings],

            'test_total': [k['test_total'] for k in warnings],
            'test_na': [k['test_na'] for k in warnings],
            'test_na_perc': [k['test_na_perc'] for k in warnings]
        }
    
class Stray:

    class MovingWindowFixedSize(Step):
        """
        Train and test in this context are interpreted as follows:
        Fixed size partitions may need data points that are not included in the given input_data (i.e., data points preceeding the first one). Such partitions are ignored (the output range_names etc. do not contain such partitions).
        """
        def __init__(self, evaluation_date_range: DateRange, train_size=99, test_size=1, window_shift=1, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.evaluation_date_range = evaluation_date_range
            self.train_size = train_size
            self.test_size = test_size
            self.window_shift = window_shift

        def run(self, input_data: ResultType) -> ResultType:
            data = input_data.output
            data_len = data.shape[0]
            assert not data.empty, "Input DataFrame is empty. Cannot extract any partition."

            evaluation_date_start = self.evaluation_date_range.start
            evaluation_date_end = self.evaluation_date_range.end
            try:
                evaluation_data_start_index = data[data.Collection_Date >= evaluation_date_start].iloc[0].name
            except IndexError:
                raise IndexError("Can't produce any partition because the evaluation_date_range correponds to an empty data frame for the given input_data")
            try:
                evaluation_data_end_index = data[data.Collection_Date >= evaluation_date_end].iloc[0].name
            except IndexError:
                evaluation_data_end_index = data.iloc[-1].name

            evaluation_data_positional_start, evaluation_data_positional_end = data.index.get_indexer([evaluation_data_start_index, evaluation_data_end_index])
            evaluation_data_positional_end += self.train_size   # makes possible to evaluate the last test sequences as many times as the other ones

            # return [data.loc[evaluation_data_start_index].Collection_Date.strftime("%G-%V-%u"), 
            #         data.loc[evaluation_data_end_index].Collection_Date.strftime("%G-%V-%u"),
            #         evaluation_data_positional_start, evaluation_data_positional_end]

            test_ranges = np.arange(evaluation_data_positional_start, evaluation_data_positional_end, self.window_shift)
            test_ranges = list(zip(test_ranges.tolist(), (test_ranges+self.test_size).tolist()))
            
            training_ranges = np.arange(evaluation_data_positional_start-self.train_size, evaluation_data_positional_end-self.train_size, self.window_shift)
            training_ranges = list(zip(training_ranges.tolist(), (training_ranges+self.train_size).tolist()))

            train_plus_test_ranges = [(tr[0],te[1]) for tr,te in zip(training_ranges, test_ranges)]

            # remove invalid partitions
            valid_partitions_k = [k for k,(li,re) in enumerate(train_plus_test_ranges) if li >= 0 and re < data_len]
            invalid_partitions_k = [k for k,(li,re) in enumerate(train_plus_test_ranges) if li < 0 or re >= data_len]
            if len(valid_partitions_k) > 0:
                valid_partitions_k_first = valid_partitions_k[0]
                valid_partitions_k_last = valid_partitions_k[-1] + 1
                
                train_plus_test_ranges = train_plus_test_ranges[valid_partitions_k_first:valid_partitions_k_last]
                training_ranges = training_ranges[valid_partitions_k_first:valid_partitions_k_last]
                test_ranges = test_ranges[valid_partitions_k_first:valid_partitions_k_last]
                if invalid_partitions_k:
                    first_seq_of_first_valid_partition = data.iloc[train_plus_test_ranges[0][0] + self.train_size + self.test_size]
                    start_date = first_seq_of_first_valid_partition.Collection_Date.strftime("%Y-%m-%d")
                    last_seq_of_last_valid_partition = data.iloc[train_plus_test_ranges[-1][-1] - self.train_size + self.test_size]
                    end_date = last_seq_of_last_valid_partition.Collection_Date.strftime("%Y-%m-%d")
                    print(f"Collection_Date of valid partitions is between {start_date} and {end_date}")
            else:
                train_plus_test_ranges = []
                training_ranges = []
                test_ranges = []

            train_range_iterator = iter(training_ranges)
            test_range_iterator = iter(test_ranges)
            train_plus_test_iterator = iter(train_plus_test_ranges)
            
                # return [data.loc[evaluation_data_start_index].Collection_Date.strftime("%G-%V-%u"), 
                #         data.loc[evaluation_data_end_index].Collection_Date.strftime("%G-%V-%u"), "\n", 
                #         evaluation_data_positional_start, evaluation_data_positional_end, "\n",
                #         list(zip(training_ranges, test_ranges))[:5], list(zip(training_ranges, test_ranges))[-5:], "\n",
                #         train_plus_test_ranges[:5], train_plus_test_ranges[-5:]]
            
            def next_test_partition():
                return next(test_range_iterator)
            
            def next_training_partition():
                return next(train_range_iterator)

            def next_train_plus_test_partition():
                return next(train_plus_test_iterator)
            
            return {
                '.next_training_partition': next_training_partition,
                '.next_test_partition': next_test_partition,
                '.next_train_plus_test_partition': next_train_plus_test_partition,
                # USE ISO 8601 for week numbering (%V) and year numbering (%G)
                'train_range_names': dict(enumerate(training_ranges)),
                'test_range_names': dict(enumerate(test_ranges)),
                'train_plus_test_range_names': dict(enumerate(train_plus_test_ranges)),
                'data': data,
                'data_feature_names': input_data['data_feature_names'],
                'metadata_feature_names': input_data['metadata_feature_names']
            }

    class OutlierDetection(Step):
        def __init__(self, alpha=0.01, k=10, knn_algortithm="brute", p = 0.5, size_threshold=50, outlier_tail="max", *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.alpha = alpha
            self.k = k
            self.knn_algorithm = knn_algortithm
            self.p = p
            self.size_threshold = size_threshold
            self.outlier_tail = outlier_tail

        def run(self, data_input: ResultType):
            X = data_input['output'][data_input['data_feature_names']]
            if X.empty:
                return {'output': []}
            else:
                X = MinMaxScaler().fit_transform(X)
                model = STRAY(self.alpha, self.k, self.knn_algorithm, self.p, self.size_threshold, self.outlier_tail)
                y = model.fit_transform(X)        # for SKTIME.annotation.STRAY
                return {
                    'output': y
                }                

    class MapOutliersToOriginalData(Step):
        def run(self, data_input: ResultType, stray_output: ResultType["Stray.OutlierDetection"]):
            md = data_input.output.copy()
            try:
                md['outlier'] = stray_output.output
            except:
                raise
            if 'scores' in stray_output:    # addition for Andrea's eaon.anomaly_detection.STRAY
                md['stray_scores'] = stray_output.scores
            output = {
                'output': md,
                'data_feature_names': data_input.data_feature_names,
                'metadata_feature_names': data_input.metadata_feature_names + ['outlier']
            }
            if 'bound_val' in stray_output: # addition for Andrea's eaon.anomaly_detection.STRAY
                output['stray_outliers_idx'] = stray_output.outliers_idx
                output['stray_bound_val'] = stray_output.bound_val
                output['stray_bound_idx'] = stray_output.bound_idx
            return output
        
    class FilterOutliers(Step):
        def run(self, data_input: ResultType["Stray.MapOutliersToOriginalData"]):
            return {
                'output': data_input.output[data_input.output.outlier == True],
                'data_feature_names': data_input.data_feature_names,
                'metadata_feature_names': data_input.metadata_feature_names
            }
    
    class CollectOutlierIdMultipleWindows(Step):
        def run(self, *outlier_annot_datasets: ResultType):
            ids = [x for step_output in outlier_annot_datasets for x in step_output.output.index]
            return {'output': ids}
    
    class AnnotateOutliersinOriginalData(Step):
        def run(self, global_data: ResultType, outliers_ids: ResultType["Stray.CollectOutlierIdMultipleWindows"]):
            md = global_data.output.copy()
            md['outlier'] = md.index.isin(outliers_ids.output)
            return {
                'output': md,
                'data_feature_names': global_data.data_feature_names, 
                'metadata_feature_names': global_data.metadata_feature_names + ['outlier']
            }
        
    class CollectOutlierCountMultipleWindows(Step):
        def run(self, *outlier_annot_datasets: ResultType):
            c = Counter()
            for ds in outlier_annot_datasets:
                outliers_index = ds.output[ds.output.outlier].index.tolist()
                c.update(outliers_index)
            return {'count': c}
        
    class AnnotateOutliersCountInOriginalData(Step):
        def run(self, global_data: ResultType, outliers_count: ResultType["Stray.CollectOutlierCountMultipleWindows"]):
            md = global_data.output.copy()
            oc = {k:v for k,v in outliers_count.count.items() if k in md.index}
            md['warnings_counter'] = oc
            md['warnings_counter'] = md['warnings_counter'].fillna(0).astype(int)
            return {
                'output': md,
                'data_feature_names': global_data.data_feature_names, 
                'metadata_feature_names': global_data.metadata_feature_names + ['warnings_counter']
            }
        

class Geolocation:
        
    class GeoapifyBatchJobRequest(Step):
        def __init__(self, geoapiy_api_key, *args, continent_bounding_box=None, cache=None, **kwargs):
            """
            cache: a DF with columns lat, lon, country_or_state, country, state, location.query.text, indexed by Location (as GISAID's attribute)
            """
            super().__init__(*args, **kwargs)
            self.geoapify_api_key = geoapiy_api_key
            self.continent_bounding_box = continent_bounding_box
            self.cache = cache if cache is not None else Geolocation.GeoapifyBatchJobRequest.get_new_cache()
        
        @staticmethod
        def get_new_cache():
            return (
                pd.DataFrame([], columns=['lat', 'lon', 'country_or_state', 'country', 'state', 'location.query.text', 'Location'])
                .astype({'lat': np.float64, 'lon': np.float64, 'country_or_state': str, 'country': str, 'state': str, 'location.query.text': str, 'Location': str})
                .set_index('Location')
            )

        def run(self, input_data: ResultType):
            """
            Returns a job if an API query is necessary, or an empty dictionary when the requested information is completely satisified by the cache
            """
            # identify country or state of a location
            unique_locations = pd.unique(input_data.output.Location)
            assert all(pd.Series(unique_locations).str.len() > 0), "Empty string Location cannot be geolocated"
            
            if not self.cache.empty:
                cacheable_unique_locations = pd.Series(unique_locations[~pd.isna(unique_locations)])
                # strip cached_results from unique_locations (it's going to be the request argument)
                unique_locations = cacheable_unique_locations[~cacheable_unique_locations.isin(self.cache.index)].tolist()
                if not unique_locations:
                    return {}   
            
            body, unique_locations2body = Geoapify.geocoding_batch_request_body_from_locations(unique_locations)
            # send job
            url = Geoapify.geocoding_batch_request(self.geoapify_api_key, self.continent_bounding_box)
            print("POST", url)
            job_response = requests.post(url, json=body, headers=Geoapify.headers)
            job_response.raise_for_status()
            job_url = job_response.json()['url']    
            return {
                'job_url': job_url,
                'unique_locations2body': unique_locations2body
            }

    class GeoapifyParseBatchRequestOutput(Step):
        def __init__(self, geoapiy_api_key, *args, cache=None, **kwargs):
            super().__init__(*args, **kwargs)
            self.geoapify_api_key = geoapiy_api_key
            self.geo_kb = GeoKB()
            self.cache = cache
                
        def run(self, geoapify_job: ResultType["Geolocation.GeoapifyBatchJobRequest"], input_data: ResultType) -> ResultType:   
            data = input_data.output
            output_metadata_columns = data.columns.tolist() + ['lat', 'lon', 'country_or_state', 'country', 'state', 'location.query.text']

            if self.cache is not None:
                augmented_data_part_1 = (pd.merge(data.reset_index(), # .reset_index is to preserve the original index (otherwise deleted on merge)
                                                self.cache, left_on="Location", right_index=True, how="inner", validate="m:1")
                                        .set_index("index"))[output_metadata_columns]
            
            if not geoapify_job:
                augmented_data = augmented_data_part_1
            else:
                output = Geoapify.post_batch_request_with_csv_output(None, None, preexisting_job_url=geoapify_job.job_url)
                country_or_usa_state = Geoapify.parse_country_or_usa_state_from_query_output(output)
                unique_locations_cooordinates = pd.concat([
                    pd.merge(country_or_usa_state[country_or_usa_state.country == 'United States'], self.geo_kb.coordinates_usa_states(), left_on='country_or_state', right_on='location_name', how="left"),
                    pd.merge(country_or_usa_state[country_or_usa_state.country != 'United States'], self.geo_kb.coordinates_countries(), left_on='country_or_state', right_on='location_name', how="left")],
                    axis=0
                )
                # join the input with the new columns (lat, lon, country_or_state, country, state and location.query.text)
                augmented_data_part_2 = data.copy() if self.cache is None else data.loc[~data.index.isin(augmented_data_part_1.index)].copy()

                augmented_data_part_2['location.query.text'] = augmented_data_part_2.Location.map(geoapify_job.unique_locations2body)
                augmented_data_part_2 = (pd.merge(augmented_data_part_2.reset_index(), unique_locations_cooordinates, left_on='location.query.text', right_on='query.text', how="inner", validate="m:1")
                                         .set_index('index')
                                         )[output_metadata_columns]
                augmented_data = augmented_data_part_2 if self.cache is None else pd.concat([augmented_data_part_1, augmented_data_part_2], axis=0, ignore_index=False).reindex(data.index)
                
                # update cache
                cache_update_chunk = (augmented_data_part_2[["lat", "lon", "country_or_state", "country", "state", "location.query.text", "Location"]]
                                    .rename(columns={'query.text':'location.query.text'})
                                    .drop_duplicates(subset='Location')
                                    .set_index("Location"))
                self.cache = cache_update_chunk if self.cache is None else pd.concat([self.cache, cache_update_chunk], axis=0)
                
            return {
                'output': augmented_data,
                'data_feature_names': input_data.data_feature_names, 
                'metadata_feature_names': output_metadata_columns
            }

    class GuessCountryOrStateFromSequenceName(Step):
        def __init__(self, geoapiy_api_key, *args, continent_bounding_box=None, cache=None, **kwargs):
            super().__init__(*args, **kwargs)
            self.geoapify_api_key = geoapiy_api_key
            self.continent_bounding_box = continent_bounding_box
            self.cache = cache

        def run(self, input_data: ResultType):            
            original_data = input_data.output.copy()
            assert 'country_or_state' in input_data.metadata_feature_names, "This Step is a backup method for seuqences where the Location attribute is not sufficiently rich. It is designed to run after a prior identification of the country or state using other methods."

            # select only sequences that needs this processing
            subdata = original_data[pd.isna(original_data.country_or_state)].copy()

            # try to extract the location
            replacement_rule = lambda match_obj: match_obj.group(1)
            identified_location = subdata.country_or_state.index.str.replace(r".*/([A-Z][a-zA-Z_\-,\.]+)/.*", replacement_rule, regex=True)
            failed_identifications = identified_location.str.contains(r"\d", regex=True, case=False)
            identified_location = pd.Series(np.where(failed_identifications, pd.NA, identified_location))

            guessed_location = identified_location.str.replace(",", ", ").str.replace("_", " ").str.replace(".", ". ")
            guessed_location.index = subdata.index
            subdata.Location = guessed_location
        
            geoapify_job_output_columns = ["lat", "lon", "country_or_state", "country", "state", "location.query.text"]
            subdata = subdata.drop(columns=geoapify_job_output_columns)
            subdata4step = CleverDict({'output': subdata, 'data_feature_names': input_data.data_feature_names, 'metadata_feature_names': [x for x in input_data.metadata_feature_names if x not in geoapify_job_output_columns]})

            job_request = Geolocation.GeoapifyBatchJobRequest(self.geoapify_api_key, continent_bounding_box=self.continent_bounding_box, cache=self.cache)
            result_job_request = CleverDict(job_request.run(subdata4step))

            job = Geolocation.GeoapifyParseBatchRequestOutput(self.geoapify_api_key, cache=self.cache)
            job_output = CleverDict(job.run(result_job_request, subdata4step))

            # assemble subdata with the rest of the data
            final_data = pd.concat([job_output.output, original_data[~pd.isna(original_data.country_or_state)]], axis=0).reindex(original_data.index)
            return {
                'output': final_data,
                'data_feature_names': input_data.data_feature_names, 
                'metadata_feature_names': input_data.metadata_feature_names
            }
