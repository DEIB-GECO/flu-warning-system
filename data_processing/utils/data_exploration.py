import pandas as pd
from tqdm import tqdm
import miniFasta

class Explorer:

    @ staticmethod
    def read_meta(*paths, verbose=True, keep_columns=[]):
        """
        Collection_Date is parsed to datetime ISO8601 (missing values are filled unless allow_incomplete_collection_dates is False)
        Index is Isolate_Id + _ + Isolate_Name (spaces are replaced with _)
        """
        meta = None
        for path in (paths):
            tmp = pd.read_excel(path)
            # select columns (selecting the columns to read in pd.read_excel isn't allowed)
            dropped_columns = [x for x in (
                               'Animal_Vaccin_Product', 'Antigen_Character', 'Update_Date',
                               'Adamantanes_Resistance_geno', 'Oseltamivir_Resistance_geno', 'Zanamivir_Resistance_geno',
                               'Peramivir_Resistance_geno', 'Other_Resistance_geno', 'Adamantanes_Resistance_pheno',
                               'Oseltamivir_Resistance_pheno', 'Zanamivir_Resistance_pheno',
                               'Peramivir_Resistance_pheno', 'Other_Resistance_pheno',
                               'Patient_Status', 'Zip_Code', 'Outbreak', 'Pathogen_Test_Info', 'Animal_Specimen_Source',
                               'Animal_Health_Status', 'Domestic_Status', 'PMID',
                               'PB2 INSDC_Upload', 'PB1 INSDC_Upload', 'PA INSDC_Upload', 'HA INSDC_Upload',
                               'NP INSDC_Upload', 'NA INSDC_Upload', 'MP INSDC_Upload', 'NS INSDC_Upload',
                               'HE INSDC_Upload', 'P3 INSDC_Upload',
                               'Publication', 'Originating_Lab', 'Originating_Sample_Id',
                               'Isolate_Submitter', 'Submitting_Lab', 'Submitting_Sample_Id',
                               'Authors', 'Note', 'Human_Specimen_Source') if x not in keep_columns and x in tmp.columns]
            tmp = tmp.drop(columns=dropped_columns)
            # concat to previous reads
            if meta is None:
                meta = tmp
            else:
                meta = pd.concat([meta, tmp])
        meta = meta.set_index(
            meta[['Isolate_Id', 'Isolate_Name']].apply(lambda row: row.Isolate_Id + "_" + row.Isolate_Name,
                                                        axis=1).str.replace(" ", "_"))
        if verbose:
            print("Explorer.read_meta:")
            print("Meta size:", meta.shape)
            print("Meta Columns:", meta.columns)
        return meta

    @staticmethod
    def read_data(*paths, verbose=True):
        """
        Index is Isolate_Id + _ + Isolate_name (spaces are replaced with _)
        """
        tmp1 = []
        isolates = set()
        for path in (paths):
            for si, fo in tqdm(enumerate(miniFasta.read(path))):
                foheader = fo.getHead()[1:]     # [1:] removes ">"

                # identify sequence
                try:
                    _,isolate_id,virus_name,type_subtype,lin,clade,segm = foheader.split("|")
                except Exception as e:
                    print(f"Fasta header unpack unusccesful on sequence n. {si} (0-based) of file {path}. Fasta header is: {foheader}.\nSequence ignored.")
                id = (isolate_id + "_" + virus_name).strip()
                isolates.add(isolate_id)

                sequence = fo.getSeq()
                seq_len = len(sequence)
                n_ratio = sequence.count("N") / seq_len

                # collect sequence header info for exploration
                tmp1.append([id, segm, seq_len, n_ratio, sequence])
        data = pd.DataFrame(tmp1, columns=['idx', 'segm','seq_len', 'n_ratio', 'sequence']).set_index('idx')
        if verbose:
            print("Data size:", data.shape, "Isolates count:", len(isolates))
            print("Data columns:", data.columns)
        return data
