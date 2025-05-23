# %%
import os 
print("Current directory:", os.getcwd())

# %%
from job_orchestra import Context, Step
import sys
sys.tracebacklimit = 0
import argparse
import shutil
from utils.pipeline_lib import *

# %% [markdown]
# # Download from GISAID
# Put the downloade isolates metadata and isolates HA genomes in the `input` directory. The software expects one .xls file and one .fasta file. 
# 
# The header of the fasta file must match the pattern: DNA Accession no. | Isolate ID | Isolate name | Type | Lineage | Clade | Segment

# %% [markdown]
# # Copy the Alignments files
# Copy the alignment files into the `alignments` directory. The software expects one .fasta file and one .insertions.csv file. 

# %% [markdown]
# # Config parameters

arg_parser = argparse.ArgumentParser(description="A script that computes warnings for the Flu Warning System (manuscript, data analyses and web application authored by Alfonsi T.; Bernasconi A.; Chiara M.; Ceri S.)")

# Positional argument (required)
arg_parser.add_argument("serotype", type=str, choices=["H1N1", "H5N1"], help="The type of influenza virus. Must be either 'H1N1' or 'H5N1'.")

# Optional flag (boolean)
# arg_parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output.")

args = arg_parser.parse_args()

# %%
def get_alignment_files(alignments_dir):
    print("Checking alignments directory...", end=" ")
    alignments_dir += os.path.sep if not alignments_dir.endswith(os.path.sep) else ""
    alignment_file_path = list(map(lambda x: alignments_dir+x, filter(lambda x: x.endswith(".fasta"), os.listdir(alignments_dir))))
    indels_file_path = list(map(lambda x: alignments_dir+x, filter(lambda x: x.endswith(".insertions.csv"), os.listdir(alignments_dir))))
    # check existence
    alignment_file_exists = len(alignment_file_path) > 0
    indels_file_exists = len(indels_file_path) > 0
    # check ambiguity
    more_alignment_files = len(alignment_file_path) > 1
    more_indels_files = len(indels_file_path) > 1
    if alignment_file_exists and not more_alignment_files and indels_file_exists and not more_indels_files:
        print("OK.")
    else:
        print("")   
        print(f"Was expected to find exactly two files in the alignment directory with extensions:\n(1) .fasta \n(2) .insertions.csv")
        if not alignment_file_exists and indels_file_exists:
            print(f"But (1) was not found.")
        elif alignment_file_exists and not indels_file_exists:
            print(f"But (2) was not found.")
        elif not alignment_file_exists and not indels_file_exists:
            print(f"But (1) and (2) were not found.")
        else:
            print("But more than two files match the criteria for (1) or (2). The input is ambiguous.")
        raise FileNotFoundError(f"Check the existence of alignments directory ({alignments_dir}) and that it contains the expected files.")
    
    # pick first and only file
    alignment_file_path = alignment_file_path[0]
    indels_file_path = indels_file_path[0]
    return alignment_file_path, indels_file_path

def get_gisaid_inputs(inputs_dir):
    print("Checking inputs directory...", end=" ")
    inputs_dir += os.path.sep if not inputs_dir.endswith(os.path.sep) else ""
    data_file_path = list(map(lambda x: inputs_dir+x, filter(lambda x: x.endswith(".fasta"), os.listdir(inputs_dir))))
    metadata_file_path = list(map(lambda x: inputs_dir+x, filter(lambda x: x.endswith(".xls"), os.listdir(inputs_dir))))
    # check existence
    data_file_exists = len(data_file_path) > 0
    metadata_file_exists = len(metadata_file_path) > 0
    # check ambiguity
    more_data_files = len(data_file_path) > 1
    more_metadata_files = len(metadata_file_path) > 1
    if data_file_exists and not more_data_files and metadata_file_exists and not more_metadata_files:
        print("OK.")
    else:
        print("")   
        print(f"Was expected to find exactly two files in the inputs directory with extensions:\n(1) .fasta \n(2) .xls")
        if not data_file_exists and metadata_file_exists:
            print(f"But (1) was not found.")
        elif data_file_exists and not metadata_file_exists:
            print(f"But (2) was not found.")
        elif not data_file_exists and not metadata_file_exists:
            print(f"But (1) and (2) were not found.")
        else:
            print("But more than two files match the criteria for (1) or (2). The input is ambiguous.")
        raise FileNotFoundError(f"Check the existence of inputs directory ({inputs_dir}) and that it contains the expected files.")
    
    # pick first and only file
    data_file_path = data_file_path[0]
    metadata_file_path = metadata_file_path[0]
    return data_file_path, metadata_file_path

def assert_dir_exists(dirname, msg):
    if not os.path.exists(dirname) and os.path.isdir(dirname):
        raise FileNotFoundError(msg)

def remove_dir(dirname):
    try:
        shutil.rmtree(dirname)
    except:
        pass

def read_stray_configurations(file_path):
    try:
        with open(file_path) as f:
            configs = [line.replace(',', ' ').split() for line in f if not line.startswith("#")]
            configs = [tuple(map(int,pair)) for pair in configs]
            assert all([len(pairs) == 2 for pairs in configs])
            assert len(configs) > 0
    except:
        raise ValueError(f"Malformed file {file_path}. Make sure the file contains at least one configuration and that each line is formatted as 'number,number'. Lines starting with # are ignored.")
    return configs

# args
serotype = args.serotype
# output
os.makedirs("output", exist_ok=True)
db_url = f'output{os.path.sep}{serotype}.sqlite'
# cache
cache_dir = ".cache"
os.makedirs(cache_dir, exist_ok=True)
ct_path = f"{cache_dir}{os.path.sep}{serotype}_analysis_ctx"
remove_dir(ct_path)
# cds_range
if serotype == "H1N1":
    cds_range = [(1,1701)]
elif serotype == "H5N1":
    cds_range = [(22,1728)]
else:
    raise ValueError(f"Unkonw serotype {serotype}")
print("CDS range", cds_range)
# prot_name
prot_name = 'HA'
# stray configurations
stray_configurations_file_path = f"inputs{os.path.sep}stray_configurations.csv"
stray_configurations = read_stray_configurations(stray_configurations_file_path)
print("Stray configurations:", *stray_configurations)
# gisaid inputs
inputs_dir = f"inputs{os.path.sep}{serotype}"
assert_dir_exists(inputs_dir, f"This script should run in the same directory containing the {inputs_dir} directory")
data_file_path, metadata_file_path = get_gisaid_inputs(inputs_dir)
# refseq
refseq_path = f"inputs{os.path.sep}references{os.path.sep}{serotype}_ha.fasta"
if not os.path.exists(refseq_path):
    raise FileNotFoundError(f"Reference sequence file not found at the expected location ({refseq_path})")
# alignments
alignments_dir = "alignments"
assert_dir_exists(alignments_dir, f"This script should run in the same directory containing the {alignments_dir} directory")
alignment_file_path = f"alignments/aligned_{serotype}_ha.fasta"
indels_file_path = f"alignments/aligned_{serotype}_ha.fasta.insertions.csv"
if not os.path.exists(alignment_file_path):
    raise FileNotFoundError(f"Alignment file (1-of-2) not found at the expected location ({alignment_file_path})")
if not os.path.exists(indels_file_path):
    raise FileNotFoundError(f"Alignment file (2-of-2) not found at the expected location ({indels_file_path})")

# %% [markdown]
# # Data Preparation

# %%
print("Context path", ct_path)
ct = Context()

# %%
class CustomMetaFilters(Step):
    def run(self, data_input):
        data = data_input.output
        return {
            'output': data[(data.Location.str.startswith("North America"))],
            'metadata_feature_names': data_input.metadata_feature_names,
            'data_feature_names': data_input.data_feature_names,
        }

S_meta_0 = Input.ReadMetadata(metadata_file_path)
S_meta_1 = Input.MetadataCleaning(allow_incomplete_collection_dates=False, allow_null_submission_dates=False, depends_on=S_meta_0)
S_meta_2 = CustomMetaFilters(depends_on=S_meta_1)
S_meta_3 = AttachWeek(depends_on=S_meta_2)
S_meta_4 = AttachCountryOrState(depends_on=S_meta_3)
S_meta = S_meta_4

# %%
# S_segm_alignment = Step.stepify({'aligned_file': alignment_file_path, 'insertions_file': indels_file_path}, name_alias=f'SegmAlignment_{prot_name}')
S_segm_alignment = Input.ReindexGISAIDFastaAlignment(alignment_file_path, indels_file_path, ctx=ct)

# %%
S_data_0 = Input.ReadData(data_file_path)
S_data_1 = Input.FilterSegm(prot_name, depends_on=S_data_0)
S_data_2 = Input.SegmCleaning(0.02, 7, depends_on=S_data_1)

S_prot_cds_0 = Input.SegmCDS(reference_file_path=refseq_path, annotation_range_1_based=cds_range[0], 
                                          depends_on=[S_data_2, S_segm_alignment], 
                                          name_alias=f"ProtCDS_raw_{prot_name}_{cds_range}")
# remove CDS with length not multiple of 3
S_prot_cds_1 = Input.FilterValidCDS(depends_on=S_prot_cds_0, name_alias=f"ProtCDS_1_{prot_name}_{cds_range}")
S_prot_cds_2 = Input.AssertCDSLength(reference_length=sum([b-a+1 for a,b in cds_range]), perc_deviation=7, depends_on=S_prot_cds_1, name_alias=f"ProtCDS_2_{prot_name}_{cds_range}", ctx=ct)
S_data = S_prot_cds_2

# %%
S_md_0 = Input.MetaDataMerge(depends_on=[S_meta, S_data])
S_md_1 = AbsoluteSortByCollectionDateAndEpisetID(depends_on=S_md_0, name_alias=f"MDFinal_{prot_name}_{cds_range}", ctx=ct)
S_md = S_md_1
print("Metadata final input dataset")
S_md.materialize().output.shape

# %%
ct.keys()

# %%
print(f"Save context in {ct_path}")
ct.store(ct_path)

# %% [markdown]
# # Outlier analysis

# %%
ct = Context.load(ct_path)
inputs = [x for x in ct.keys() if x.startswith("MDFinal_")]
print("Outlier analysis inputs", inputs)

# %%
def cds_features(ct_path: str, input_name):   
    ct = Context.load(ct_path)

    # load protein CDS from Context
    S_cds_0 = Step(name_alias=input_name, ctx=ct)
    S_cds_1 = Input.MoveDataFeaturesToMetadata(depends_on=S_cds_0)
    S_features = S_cds_1

    # compute features
    S_transformation_0 = Input.RSCU(depends_on=S_features)
    S_transformation_1 = DataTransform.LogBySynCount(depends_on=S_transformation_0)
    S_transformation_2 = DataTransform.PlainAndLogDinucleotide(depends_on=S_transformation_1, ctx=ct)  # <- its saved in context
    cds_data = S_transformation_2.materialize()

    # Assert the DF is still sorted by Sort_Key
    pd.testing.assert_index_equal(cds_data.output.index, S_cds_0.materialize().output.index, check_order=True)

    return cds_data

def merge_input_output(input_name, stray_outlier_table):
    a = ct[input_name].output
    b = stray_outlier_table.drop(columns=stray_outlier_table.columns.intersection(a.columns))
    pd.testing.assert_index_equal(a.index,b.index, check_order=True)        # check inputs have same rows
    m = a.join(b)
    pd.testing.assert_index_equal(m.index,a.index, check_order=True)        # check output is sorted as inputs
    return m
    
def get_anomlies(n, k, input_data: ResultType):
    # collect data dependencies
    data = input_data.output.copy()
    data_feature_names = input_data.data_feature_names
    assert set(data_feature_names) == set(GeneticCode.valid_codons + list(GeneticCode.dinucleotides))
    metadata_feature_names = input_data.metadata_feature_names

    # windows
    windows = [(i,i+n) for i in range(0, data.shape[0]-n+1)]    # left included, right excluded

    # output 
    outliers_idx = []
    for ws,we in windows:
        window_input_data = {'output': data.iloc[ws:we,:], 'data_feature_names': data_feature_names, 'metadata_feature_names': metadata_feature_names}
        window_anomalies = Stray.OutlierDetection(k=k).run(window_input_data)['output'].tolist()
        if window_anomalies[-1]:
            outliers_idx.append(window_input_data['output'].index[-1])
    data['outlier'] = np.zeros(data.shape[0], dtype=bool)
    data.loc[outliers_idx,'outlier'] = np.True_

    print("Outliers:", len(outliers_idx))
    return data


def loop(stray_configurations, input_data: ResultType):
    for i, (n,k) in enumerate(stray_configurations):
        print(f"Computing configuration {i+1}/{len(stray_configurations)}: (N,K) {n},{k}...")
        return get_anomlies(n,k,input_data)    

print(f"Using data from Context {ct_path}")
for input_name in inputs:
    input_data_features = cds_features(ct_path, input_name)    
    print(f"Running anomaly detection on input: {input_name}")
    for i, (n,k) in enumerate(stray_configurations):
        output_name = f"Outlier_{(n,k)}_" + input_name.replace("MDFinal_", "", 1)
        # avoid recomputing saved results in case the pipeline gets restarted
        if output_name in ct:
            continue
        print(f"Computing configuration {i+1}/{len(stray_configurations)}: (N,K) {n},{k}...", end=" ")
        cds_anomalies = get_anomlies(n,k,input_data_features)    

        # merge result with full CDS input_data
        merged_result = merge_input_output(input_name, cds_anomalies)
        pd.testing.assert_index_equal(merged_result.index, ct[input_name].output.index)                 # check output is sorted as input
        pd.testing.assert_series_equal(merged_result.Sort_Key, merged_result.Sort_Key.sort_values())  # check Sort_Key is sorted as expected
        # save result
        Step.stepify({'output': merged_result, 
                    'metadata_feature_names': ct[input_name].metadata_feature_names, 
                    'data_feature_names': GeneticCode.valid_codons + list(GeneticCode.dinucleotides) + ['outlier']
                    }, name_alias=output_name, ctx=ct).materialize()


# %% [markdown]
# # Write database

# %%
outlier_out = [x for x in ct.keys() if x.startswith("Outlier_")]
print("Outlier analysis outputs", outlier_out)
outlier_out2stray_config = {k: k.split("_")[1].strip("()").split(", ") for k in outlier_out}

# %%
input_data_table_exported_columns = ct[outlier_out[0]].metadata_feature_names
input_data_table = ct[outlier_out[0]].output[input_data_table_exported_columns]
input_data_table_idx = input_data_table.index
for i in range(1, len(outlier_out)):
    assert input_data_table_idx.symmetric_difference(ct[outlier_out[i]].output.index).empty     # check input_data_table contains all isolates

outlier_tables_exported_columns = input_data_table_exported_columns + ['outlier']

# %%
# cnx.dispose()

# %%
from sqlalchemy import create_engine

print("Write database at", db_url)
cnx = create_engine("sqlite:///"+db_url)

# %%
with cnx.connect() as connection:
    input_data_table.to_sql(name='input_data', con=connection, if_exists='replace')

    for outlier_analysis_name in outlier_out:       
        stray_config = outlier_out2stray_config[outlier_analysis_name]
        db_table_name = f"window{stray_config[0]}_k{stray_config[1]}"

        outlier_table = ct[outlier_analysis_name].output[outlier_tables_exported_columns]
        outlier_table = outlier_table[outlier_table.outlier]
        outlier_table.to_sql(name=db_table_name, con=connection, if_exists='replace')


# %%
cnx.dispose()

# %% [markdown]
# # Files cleanup

# %%
# rewritten alignments
for k,v in S_segm_alignment.materialize().items():
    try:
        os.remove(v)
    except:
        pass
# cache
remove_dir(ct_path)
