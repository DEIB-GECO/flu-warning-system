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
from sktime.annotation.stray import STRAY
from sklearn.preprocessing import MinMaxScaler

from utils.location import Geoapify, GeoKB
from utils.genetic_code import GeneticCode
from utils.data_cleaning import DataCleaning, DateRange
from utils.pipeline import Step

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
  
    class AssignHostType(Step):
        def run(input_data):
            meta = input_data.output
            meta = DataCleaning.merge_fix_host_species(meta)
            meta = DataCleaning.assign_host_type(meta, additional_mappings=self.additional_mappings, raise_on_error=True)
            meta = meta[meta.Host_Type != "other"].copy()
            return {
                'output': meta, 
                'data_feature_names': input_data.data_feature_names,
                'metadata_feature_names': input_data.metadata_feature_names + ['Host_Type']
            }

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

class AttachWeek(Step):
    def run(self, input_data: ResultType):
        output = input_data.output.copy()
        output['Week'] = output['Collection_Date'].dt.strftime("%G-%V")
        return {
            'output': output,
            'metadata_feature_names': input_data.metadata_feature_names + ['Week'],
            'data_feature_names': input_data.data_feature_names
        }

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
