import pandas as pd
import plotly.graph_objects as go
import numpy as np


class DataCleaning:

    @staticmethod
    def segment_length_stats(data, histogram_xaxis_range: list = None):
        """
        Prints aggregated stats about the length of each segment in data
        """
        for segm, segm_data in data.groupby(by="segm"):
            print(*list(zip(
                ["SEGMENT", "count", "q2", "q3", "median", "mean", "median-mean", "(1-median/mean)%"],
                [
                    segm,
                    segm_data.shape[0],
                    segm_data.sequence.str.len().quantile(0.5),
                    segm_data.sequence.str.len().quantile(0.75),
                    f"{segm_data.sequence.str.len().median():.2f}",
                    f"{segm_data.sequence.str.len().mean():.2f}",
                    f"{segm_data.sequence.str.len().median() - segm_data.sequence.str.len().mean():.2f}",
                    f"{(1 - segm_data.sequence.str.len().median() / segm_data.sequence.str.len().mean()) * 100:.2f}"
                ])), sep="\n")
        segments = data.segm.unique()
        fig = go.Figure()
        for segm in segments:
            fig.add_trace(go.Histogram(x=data[data.segm == segm].seq_len, name=segm))
        fig.update_layout(title="Histogram of sequence length")
        if histogram_xaxis_range:
            fig.update_xaxes(range=histogram_xaxis_range)
        fig.show()

    @staticmethod
    def simul_filter_segment_length(data, segm, lower_b, upper_b):
        """
        Simulate the final size of the data if only the sequences of segment segm
        having lower_b <= length <= upper_b were selected.
        """
        segm_data = data[data.segm == segm]
        segm_included_sequences_count = segm_data[(segm_data.sequence.str.len() <= upper_b) & (segm_data.sequence.str.len() >= lower_b)].shape[0]
        print(segm, '\t', lower_b, upper_b, segm_included_sequences_count, segm_data.shape[0], segm_included_sequences_count / segm_data.shape[0] * 100)

    @staticmethod
    def remove_duplicate_index(data, verbose=True):
        total_sequences = data.shape[0]
        tmp1 = {segm: segm_data for segm, segm_data in data.groupby(by="segm")}
        segments_affected = []
        for segm, segm_data in tmp1.items():
            if segm_data.index.duplicated().any():
                segments_affected.append(segm)
            tmp1[segm] = segm_data.iloc[np.unique(segm_data.index.values, return_index=True)[1]]
        # reconcile segments
        data = pd.concat(list(tmp1.values()))
        if verbose:
            print("Removed", total_sequences - data.shape[0], "duplicate sequences found in segments", ", ".join(segments_affected))
            print("Data size:", data.shape[0])
        return data

    @staticmethod
    def assign_Y_group(meta):
        try:
            meta = meta.drop(columns=['Y_group'])
        except KeyError:
            pass
        tmp2 = {group_index: group for group_index, group in
                meta.groupby(pd.Grouper(key='Collection_Date', axis=0, freq='Y'))}
        # add group index to df
        for group_index, group in tmp2.items():
            group['Y_group'] = group_index
        # concat df of groups
        meta = pd.concat([group for _, group in tmp2.items()])
        return meta

    @staticmethod
    def assign_M_group(meta):
        try:
            meta = meta.drop(columns=['M_group'])
        except KeyError:
            pass
        tmp2 = {group_index: group for group_index, group in
                meta.groupby(pd.Grouper(key='Collection_Date', axis=0, freq='M'))}
        # add group index to df
        for group_index, group in tmp2.items():
            group['M_group'] = group_index
        # concat df of groups
        meta = pd.concat([group for _, group in tmp2.items()])
        return meta
    
    animal_category_mapping = {
            'accipiter cooperii': 'wild bird',
            'accipiter gentilis': 'wild bird',
            'accipiter nisus': 'wild bird',
            'accipiter trivirgatus': 'wild bird',
            'branta hutchinsii': 'wild bird',
            'aix galericulata': 'wild bird',
            'aix sponsa': 'wild bird',
            'alectoris chukar': 'wild bird',
            'american wigeon': 'wild bird',
            'anas acuta': 'wild bird',
            'anas americana': 'wild bird',
            'anas boschas': 'wild bird',
            'anas carolinensis': 'wild bird',
            'anas clypeata': 'wild bird',
            'anas crecca': 'wild bird',
            'anas cyanoptera': 'wild bird',
            'anas discors': 'wild bird',
            'anas falcata': 'wild bird',
            'anas flavirostris': 'wild bird',
            'anas formosa': 'wild bird',
            'anas georgica': 'wild bird',
            'anas penelope': 'wild bird',
            'anas platyrhynchos': 'wild bird',
            'anas platyrhynchos f. domestica': 'domestic bird',
            'anas platyrhynchos var. domesticus': 'domestic bird',
            'anas poecilorhyncha': 'wild bird',
            'anas querquedula': 'wild bird',
            'anas rubripes': 'wild bird',
            'anas sp.': 'wild bird',
            'anas strepera': 'wild bird',
            'anas undalata': 'wild bird',
            'anas zonorhyncha': 'wild bird',
            'animal': 'other',
            'anser albifrons': 'wild bird',
            'anser anser': 'wild bird',
            'anser anser domesticus': 'domestic bird',
            'anser brachyrhynchus': 'wild bird',
            'anser caerulescens': 'wild bird',
            'anser cygnoides': 'wild bird',
            'anser fabalis': 'wild bird',
            'anser indicus': 'wild bird',
            'anser rossii': 'wild bird',
            'anser sp.': 'wild bird',
            'anseriformes sp.': 'wild bird',
            'arenaria interpres': 'wild bird',
            'avian': 'wild bird',
            'aythya affinis': 'wild bird',
            'aythya collaris': 'wild bird',
            'aythya ferina': 'wild bird',
            'aythya fuligula': 'wild bird',
            'aythya marila': 'wild bird',
            'aythya nyroca': 'wild bird',
            'bean goose': 'wild bird',
            'black-headed gull': 'wild bird',
            'blue-winged teal': 'wild bird',
            'branta bernicla': 'wild bird',
            'branta canadensis': 'wild bird',
            'branta leucopsis': 'wild bird',
            'bucephala clangula': 'wild bird',
            'buteo buteo': 'wild bird',
            'buteo jamaicensis': 'wild bird',
            'buteo japonicus': 'wild bird',
            'buteo lineatus': 'wild bird',
            'cairina moschata': 'domestic bird',
            'calidris alba': 'wild bird',
            'calidris canutus': 'wild bird',
            'calidris minutilla': 'wild bird',
            'canine': 'non-human mammal',
            'chen caerulescens': 'wild bird',
            'chicken': 'domestic bird',
            'chlidonias hybridus': 'wild bird',
            'chroicocephalus cirrocephalus': 'wild bird',
            'chroicocephalus ridibundus': 'wild bird',
            'circus aeruginosus': 'wild bird',
            'common teal': 'wild bird',
            "cooper's hawk": 'wild bird',
            'copsychus saularis': 'wild bird',
            'cormorant': 'wild bird',
            'corvus': 'wild bird',
            'corvus frugilegus': 'wild bird',
            'corvus macrorhynchos': 'wild bird',
            'coturnix': 'wild bird',
            'coturnix sp.': 'wild bird',
            'crow': 'wild bird',
            'curlew': 'wild bird',
            'cygnus atratus': 'wild bird',
            'cygnus columbianus': 'wild bird',
            'cygnus cygnus': 'wild bird',
            'cygnus melancoryphus': 'wild bird',
            'cygnus olor': 'wild bird',
            'dendrocygna autumnalis': 'wild bird',
            'dendrocygna viduata': 'wild bird',
            'domestic goose': 'domestic bird',
            'dove': 'wild bird',
            'duck': 'wild bird',
            'eagle': 'wild bird',
            'egyptian goose': 'domestic bird',
            'emperor goose': 'wild bird',
            'environment': 'other',
            'equine': 'non-human mammal',
            'eurasian curlew': 'wild bird',
            'falco': 'wild bird',
            'falco peregrinus': 'wild bird',
            'falco tinnunculus': 'wild bird',
            'falcon': 'wild bird',
            'falsco rusticolus': 'wild bird',
            'feces': 'other',
            'feline': 'non-human mammal',
            'felis catus': 'non-human mammal',
            'ferret': 'non-human mammal',
            'gallinago gallinago': 'wild bird',
            'gallinula chloropus': 'wild bird',
            'gallus gallus': 'domestic bird',
            'gallus gallus domesticus': 'domestic bird',
            'glaucous gull': 'wild bird',
            'glaucous-winged gull': 'wild bird',
            'goose': 'wild bird',
            'gracula religiosa': 'wild bird',
            'great crested grebe ': 'wild bird',
            'great tit': 'wild bird',
            'green-winged teal': 'wild bird',
            'grey teal': 'wild bird',
            'greylag goose': 'wild bird',
            'guinea fowl': 'domestic bird',
            'guineafowl': 'domestic bird',
            'gull': 'wild bird',
            'haliaeetus leucocephalus': 'wild bird',
            'halietus albicilla': 'wild bird',
            'halietus leucocephalus': 'wild bird',
            'herring gull': 'wild bird',
            'himantopus himantopus melanurus': 'wild bird',
            'hirundo rustica': 'wild bird',
            'host': 'other',
            'human': 'human',
            'insect': 'other',
            'laboratory derived': 'other',
            'larosterna inca': 'wild bird',
            'larus': 'wild bird',
            'larus argentatus': 'wild bird',
            'larus armenicus': 'wild bird',
            'larus brunnicephalus': 'wild bird',
            'larus cachinnans': 'wild bird',
            'larus canus': 'wild bird',
            'larus delawarensis': 'wild bird',
            'larus dominicanus': 'wild bird',
            'larus fuscus': 'wild bird',
            'larus glaucescens': 'wild bird',
            'larus ichthyaetus': 'wild bird',
            'larus marinus': 'wild bird',
            'larus melanocephalus': 'wild bird',
            'larus ridibundus': 'wild bird',
            'larus schistisagus': 'wild bird',
            'larus smithsonianus': 'wild bird',
            'leucophaeus': 'wild bird',
            'leucophaeus atricilla': 'wild bird',
            'lophodytes cucullatus': 'wild bird',
            'lophura nycthemera': 'wild bird',
            'mallard': 'wild bird',
            'mallard duck': 'wild bird',
            'mareca penelope': 'wild bird',
            'meerkat': 'non-human mammal',
            'meleagris gallopavo': 'wild bird',
            'mink': 'non-human mammal',
            'morphnus guianensis': 'wild bird',
            'mulard duck': 'wild bird',
            'necrosyrtes monachus': 'wild bird',
            'nisaetus nipalensis': 'wild bird',
            'northern pintail': 'wild bird',
            'northern shoveler': 'wild bird',
            'numenius arquata': 'wild bird',
            'numida meleagris': 'wild bird',
            'oreortyx': 'wild bird',
            'ostrich': 'wild bird',
            'other avian': 'wild bird',
            'other environment': 'other',
            'other mammals': 'non-human mammal',
            'parabuteo': 'wild bird',
            'parabuteo unicinctus': 'wild bird',
            'partridge': 'wild bird',
            'passerine': 'wild bird',
            'pavo cristatus': 'wild bird',
            'penguin': 'wild bird',
            'phasanius colchicus': 'wild bird',
            'phasanius sp.': 'wild bird',
            'phasianus': 'wild bird',
            'phasianus colchicus': 'wild bird',
            'pheasant': 'wild bird',
            'pica': 'wild bird',
            'pigeon': 'wild bird',
            'podiceps cristatus': 'wild bird',
            'polyplectron bicalcaratum': 'wild bird',
            'poultry': 'domestic bird',
            'primate': 'non-human mammal',
            'pygoscelis antarcticus': 'wild bird',
            'quail': 'wild bird',
            'rails': 'wild bird',
            'ring-necked duck': 'wild bird',
            'rissa tridactyla': 'wild bird',
            'rodent': 'non-human mammal',
            'ruddy shelduck': 'wild bird',
            'ruddy turnstone': 'wild bird',
            'rynchops niger': 'wild bird',
            'sacred ibis': 'wild bird',
            'sandpiper': 'wild bird',
            'scolopax rusticola': 'wild bird',
            'seal': 'non-human mammal',
            'silver teal': 'wild bird',
            'somateria mollissima': 'wild bird',
            'somateria spectabilis': 'wild bird',
            'speckled pigeon': 'wild bird',
            'sterna hirundo': 'wild bird',
            'sterna paradisaea': 'wild bird',
            'sterna sandvicensis': 'wild bird',
            'surface swab': 'other',
            'swan': 'wild bird',
            'swine': 'non-human mammal',
            'tachybaptus ruficollis': 'wild bird',
            'tadorna feruginea': 'wild bird',
            'tadorna tadorna': 'wild bird',
            'teal': 'wild bird',
            'turkey': 'wild bird',
            'us quail': 'wild bird',
            'unknown': 'other',
            'water sample': 'other',
            'air sample': 'other',
            'nan': 'unknown',
            'numenius phaeopus': 'wild bird',
            'white-fronted goose': 'wild bird',
            'whooper swan': 'wild bird',
            'wild bird': 'wild bird',
            'wild birds': 'wild bird',
            'wild waterfowl': 'wild bird',
            'zosterops japonicus': 'wild bird',
            'mammals': 'non-human mammal',
            'american black duck': 'wild bird', 
            'anas platyrhynchos x anas acuta': 'wild bird', 
            'anser canagica': 'wild bird', 
            'aythya valisineria': 'wild bird', 
            'bucephala albeola': 'wild bird', 
            'gelochelidon nilotica': 'wild bird', 
            'pig': 'non-human mammal', 
            'sus scrofa': 'non-human mammal', 
            'sus scrofa domesticus': 'non-human mammal', 
            'sus scrofa scrofa': 'non-human mammal',
            'calidris alpina': 'wild bird',
            'camel': 'non-human mammal', 
            'cinnamon teal': 'wild bird',
            'mouse': 'non-human mammal', 
            'passer montanus': 'wild bird',
            'peafowl': 'wild bird',
            'aythya americana': 'wild bird',
            'numenius phaeopus': 'wild bird',
            'panda': 'non-human mammal',
            'bovine': 'non-human mammal',
            'dairy cow': 'non-human mammal',
            'buteo': 'wild bird',
            'chroicocephalus': 'wild bird',
            'circus': 'wild bird',
            'gyps fulvus': 'wild bird',
            'shearwater': 'wild bird',
            'turnstone': 'wild bird',
            'numida sp.': 'wild bird',
            'dabbling duck': 'wild bird'
        }
    
    @staticmethod
    def merge_fix_host_species(df):
        df.Host = df.Host.str.lower()
        # domestic birds
        df.loc[df.Host.str.startswith("gallus"), 'Host'] = 'chicken'
        df.loc[df.Host.str.startswith("poultry"), 'Host'] = 'chicken'
        df.loc[df.Host.str.startswith("gallus"), 'Host'] = 'chicken'
        df.loc[df.Host.str.startswith("anas platyrhynchos var. domesticus"), 'Host'] = 'anas platyrhynchos f. domestica'
        df.loc[df.Host.str.startswith("guinea"), 'Host'] = 'guineafowl'
        # non-human mammals
        df.loc[df.Host == "other mammals", 'Host'] = 'mammals'
        df.loc[df.Host == "sus scrofa", 'Host'] = 'swine'
        # wild birds
        df.loc[df.Host == "wild bird", 'Host'] = 'wild birds'
        return df
    
    @staticmethod
    def print_known_hosts_classified_by_host_type():
        for t in sorted(set(DataCleaning.animal_category_mapping.values())):
            print("*****", t)
            print(*sorted({x for x,y in DataCleaning.animal_category_mapping.items() if y == t}), sep=",\n")

    @staticmethod
    def assign_host_type(meta, additional_mappings: dict = None, raise_on_error = True):
        """
        Returns hte input dataframe with a new column 'Host_Type' classifying the 'Host' into the categories:
        wild bird, domestic bird, human, non-human mammal, other.
        The mapping is case insensitive
        """
        if isinstance(additional_mappings, dict):
            DataCleaning.animal_category_mapping.update(additional_mappings)
        meta_hosts = set(meta.Host.dropna().str.lower().unique().astype(str))       # take unique values of not null species
        mapped_hosts = set(DataCleaning.animal_category_mapping.keys())
        if raise_on_error and not meta_hosts.issubset(mapped_hosts):
            raise AssertionError("Mapping for the following Hosts is not available:", sorted(meta_hosts.difference(mapped_hosts)))
        meta['Host_Type'] = meta.Host.str.lower().apply(lambda x: DataCleaning.animal_category_mapping.get(x, pd.NA) if pd.notnull(x) else x)
        return meta

    @staticmethod
    def assign_flu_season(meta):
        def season_mapping(date: pd.Timestamp):
            if date.month < 10:
                season = f"{date.year - 1}-{date.year}"
            else:
                season = f"{date.year}-{date.year + 1}"
            return season

        meta['Flu_Season'] = meta.Collection_Date.apply(season_mapping)
        return meta

    @staticmethod
    def transform_series_range(series, new_min, new_max):
        old_max = series.max()
        old_min = series.min()
        old_range = old_max - old_min
        new_range = new_max - new_min
        # shift so that any min value starts from 0, then divide by the old range -> new range is from 0 to 1
        # multiply the range 0-1 by the new range, then shift everything to so that any 0 starts from the new_min value
        out = series.apply(lambda x: ((x - old_min) / old_range) * new_range + new_min)     
        return out

    @staticmethod
    def remove_data_with_invalid_cds(data, min_cds_len, genetic_code_class):
        old_data_size = data.shape[0]
        data = data.copy()
        data['cds_ranges'] = data.apply(lambda x: list(genetic_code_class.cds_ranges(x.sequence)), axis=1)
        data['len_cds_list'] = data['cds_ranges'].map(len)
        data = data[data['len_cds_list'] > 0]
        data['cds_len'] = data.cds_ranges.apply(lambda x: [y[1] - y[0] for y in x])
        data['longest_cds_idx'] = data['cds_len'].apply(lambda z: np.argmax(z))
        data['longest_cds_len'] = data.apply(lambda x: x.cds_len[x.longest_cds_idx], axis=1)
        data['ok_cds_idx'] = data.apply(lambda x: x.longest_cds_idx if x.cds_len[x.longest_cds_idx] > min_cds_len else -1, axis=1)
        data = data[data['ok_cds_idx'] >= 0]
        print(f"Removed {old_data_size - data.shape[0]} sequences with invalid CDS")
        data['longest_cds'] = data.apply(lambda x: x.sequence[x.cds_ranges[x.longest_cds_idx][0]:x.cds_ranges[x.longest_cds_idx][1]], axis=1)
        return data


class DateRange:
    def __init__(self):
        self.start = None
        self.end = None

    def __str__(self):
        return f"{self.start}_{self.end}"

    @staticmethod
    def from_dates(start: str = '2000-01-01', end: str = 'today'):
        """
        The date range is defined as [start,end).
        :param str start: the first date included in the range
        :param str end: the closing date of the range, that should not be included
        """
        obj = DateRange()
        obj.start = pd.to_datetime(start)
        obj.end = pd.to_datetime(end)
        return obj

    @staticmethod
    def from_iso_weeks(start: str, end: str):
        start += "-1" if not start.endswith("-1") else ""
        end += "-1" if not end.endswith("-1") else ""
        obj = DateRange()
        obj.start = pd.to_datetime(start, format="%G-%V-%u") #2023-46-1: Monday of 46th week of 2023 according to ISO8601
        obj.end = pd.to_datetime(end, format="%G-%V-%u")     
        return obj
    
    @staticmethod
    def from_iso_weeks_right_included(start: str, end: str):
        start += "-1" if not start.endswith("-1") else ""
        end += "-7" if not end.endswith("-1") else ""
        obj = DateRange()
        obj.start = pd.to_datetime(start, format="%G-%V-%u") #2023-46-1: Monday of 46th week of 2023 according to ISO8601
        obj.end = pd.to_datetime(end, format="%G-%V-%u") + pd.Timedelta(days=1)
        return obj

    @staticmethod
    def week_step_date_ranges(last_datetime: pd.Timestamp = pd.to_datetime('today'), start_datetime: pd.Timestamp = pd.to_datetime('2000-01-01'), weeks=16):
        """
        Returns pairs of <start,end> pd.Timestamp spanning the given number of weeks. Edges of the range should be considered as left-included, right-excluded. 
        The argument last_datetime is always included in the last range. The argument start_datetime is always included in the first range and may fall well after the 
        left-edge of the first range. 
        :param pd.Timestamp last_datetime: last date (most recent) that must be included in the range. This date is included by the last date range returned and the 
        start date of the range is positioned the argment number of weeks backward.
        :param pd.Timestamp start_datetime: date included in the oldest date range. 
        :param int weeks: the number of weeks defining the date range width.
        """
        end_interval = last_datetime.round(freq="D") + pd.Timedelta(days=1) # is excluded (round(freq="D") sets to 00:00, +1d includes the given date)
        start_interval = end_interval
        output = []
        while start_interval > start_datetime:
            start_interval = end_interval - pd.Timedelta(weeks=weeks)
            output.append((start_interval, end_interval))
            end_interval = start_interval
        return output[::-1]
    
    @staticmethod
    def week_step_date_ranges_split_train_test_moving_window(last_date: pd.Timestamp = pd.to_datetime('today'), start_date: pd.Timestamp = pd.to_datetime('2000-01-01'), train_weeks=16, test_weeks=4, window_weeks=4):
        """
        Returns paired time intervals. The first set of intervals is the training set. The second set is the test set. 
        The i-th interval in the test set always follows the i-th interval in the training set. 
        An interval is expressed as a tuple of <start,end> pd.Timestamp objects.
        The i-th+1 interval starts after the number of weeks specified in window_weeks. If this number is equal to training_weeks+test_weeks, 
        the i-th+1 interval does not overlap with the i-th interval. 
        :param pd.Timestamp last_datetime: last date (most recent). The produced ranges will stop before this date
        :param pd.Timestamp start_datetime: date included in the oldest date range. 
        :param int train_weeks: the number of weeks defining the date range width.
        :param int test_weeks: the number of weeks defining the date range width.
        :param int window_weeks: the number of weeks defining the date range width.
        """
        output_train, output_test = [], []
        window_start = start_date
        while(window_start + pd.Timedelta(weeks=train_weeks+test_weeks)) < last_date:
            train_end = window_start + pd.Timedelta(weeks=train_weeks)
            output_train.append((window_start, train_end))
            test_start = train_end
            test_end = test_start + pd.Timedelta(weeks=test_weeks)
            output_test.append((test_start, test_end))
            window_start = window_start + pd.Timedelta(weeks=window_weeks)
        return output_train, output_test
    
    @staticmethod
    def week_step_date_ranges_split_train_test_moving_window_v2(date_range, train_weeks=16, test_weeks=4, window_shift=4):
        """
        As week_step_date_ranges_split_train_test_moving_window, but last_date and start_date are referred to the test window.
        Receives a date range of type DateRange
        The first window is the week containing the start date, the last window the last possible window that do not iclude the end date of the range
        """
        output_train, output_test = [], []
        start_date, end_date = date_range.start, date_range.end

        # ensure start_date begins with a monday and last_date ends with monday (excluded)
        if start_date.isoweekday() > 1:
            start_date -= pd.Timedelta(days=start_date.isoweekday() - 1)
        if end_date.isoweekday() > 1:
            end_date -= pd.Timedelta(days=end_date.isoweekday() - 1)

        delta_train = pd.Timedelta(weeks=train_weeks)
        delta_test = pd.Timedelta(weeks=test_weeks)
        delta_shift = pd.Timedelta(weeks=window_shift)

        window_start = start_date
        while window_start + delta_test < end_date:
            # training windows
            train_start = window_start - delta_train
            train_end = window_start
            # test windows
            test_start = window_start
            test_end = test_start + delta_test
            # update rule
            output_train.append((train_start, train_end))
            output_test.append((test_start, test_end))
            window_start += delta_shift
        return output_train, output_test
    
    @staticmethod
    def week_step_date_ranges_split_train_test_moving_window_v3(date_range, train_weeks=16, test_weeks=4, window_shift=4):
        """
        As week_step_date_ranges_split_train_test_moving_window, but last_date and start_date are referred to the test window.
        Receives a date range of type DateRange
        The first window is the week containing the start date, the last window the last possible window that do not iclude the end date of the range
        """
        output_train, output_test = [], []
        start_date, end_date = date_range.start, date_range.end

        # ensure start_date begins with a monday and last_date ends with monday (excluded)
        if start_date.isoweekday() > 1:
            start_date -= pd.Timedelta(days=start_date.isoweekday() - 1)
        if end_date.isoweekday() > 1:
            end_date -= pd.Timedelta(days=end_date.isoweekday() - 1)

        delta_train = pd.Timedelta(weeks=train_weeks)
        delta_test = pd.Timedelta(weeks=test_weeks)
        delta_shift = pd.Timedelta(weeks=window_shift)

        window_start = start_date
        while window_start + delta_test <= end_date:
            # training windows
            train_start = window_start - delta_train
            train_end = window_start
            # test windows
            test_start = window_start
            test_end = test_start + delta_test
            # update rule
            output_train.append((train_start, train_end))
            output_test.append((test_start, test_end))
            window_start += delta_shift
        return output_train, output_test
    
