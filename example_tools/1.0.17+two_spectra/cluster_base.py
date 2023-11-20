import sys
import glob
import os
import bisect
from itertools import groupby
import math

class SpectrumBase(object):

    """Base class for all spectrum-like objects"""

    def __init__(self, scan_number, file_id, precursor_mz, charge):
        self.scan_number = scan_number
        self.file_id = file_id
        self.precursor_mz = precursor_mz
        self.charge = charge

    def __eq__(self, other):
        equality = False
        if issubclass(self.__class__, SpectrumBase) and issubclass(other.__class__, SpectrumBase):
            equality = (
                str(self.scan_number) == str(other.scan_number) and
                self.file_id == other.file_id
            )
        return equality

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return "Base Spectrum (Scan number: {base.scan_number})".format(base=self)

class Cluster(SpectrumBase):

    """Class for individual cluster."""

    def __init__(self, scan_number, file_id, precursor_mz, charge, spectra,
                 purity, consensus_peptide, sqs = 0, ms_cluster_charge = None,
                 spec_prob = None, p_value = None, fdr = None, pep_fdr = None,
                 peaks_at_k = None, top_peaks = None, mix_score = 0):
        self.spectra = spectra
        self.purity = purity
        self.consensus_peptide = consensus_peptide
        self.sqs = sqs
        self.spectrum_weights = self.set_spectrum_weights()
        self.ms_cluster_charge = ms_cluster_charge
        self.spec_prob = spec_prob
        self.p_value = p_value
        self.fdr = fdr
        self.pep_fdr = pep_fdr
        self.of_split = 1
        self.peaks_at_k = peaks_at_k
        self.top_peaks = top_peaks
        self.mix_score = mix_score
        super(Cluster, self).__init__(scan_number, file_id, precursor_mz, charge)
    def avg_spectra_sqs(self):
            output = 0.0
            all_sqs = [
                spectrum.sqs
                for spectrum in self.spectra
            ]
            try:
                 output = sum(all_sqs)/float(len(all_sqs))
            except:
                pass
            return output
    def set_spectrum_weights(self):
        return [
            float(spectrum.precursor_mz)
            for spectrum in self.spectra
        ]
    def guess_charge(self):
        return [
            spectrum.charge
            for spectrum in self.spectra
        ][0]
    def mean_kl(self):
        spec_kl = [
            spectrum.kl_strict
            for spectrum in self.spectra
            if spectrum.kl_strict
        ]
        try:
            return sum(spec_kl)/len(spec_kl)
        except:
            return None
    def max_kl(self):
        spec_kl = [
            spectrum.kl_strict
            for spectrum in self.spectra
            if spectrum.kl_strict
        ]
        try:
            return max(spec_kl)
        except:
            return None


class Spectrum(SpectrumBase):

    """Class for individual spectrum (or cluster solely for SQS)."""

    def __init__(self, scan_number, file_id, precursor_mz, charge, similarity,
                 p_value, sqs, peptide, cluster_number = None, spec_prob = None,
                 fdr = None, pep_fdr = None, kl_strict = None, kl_non_strict = None, kl_delta = None, peaks = [], second_id = None, cluster_member = 'C', cosine = 0,signal_peaks = 0, order = 0):
        self.similarity = similarity
        self.p_value = p_value
        self.sqs = sqs
        self.peptide = peptide
        self.cluster_number = cluster_number
        self.spec_prob = spec_prob
        self.fdr = fdr
        self.pep_fdr = pep_fdr
        self.kl_strict = kl_strict
        self.kl_non_strict = kl_non_strict
        self.kl_delta = kl_delta
        self.peaks = peaks
        self.second_id = second_id
        self.orig_pm = 0
        self.cluster_member = cluster_member
        self.cosine = cosine
        self.signal_peaks = signal_peaks
        self.explained_intensity = 0
        self.explained_peaks = 0
        self.order = 0
        super(Spectrum, self).__init__(scan_number, file_id, precursor_mz, charge)
    def output_pair(self):
        return str(self.file_id) + ":" + str(self.scan_number) + ":" + str(self.sqs) + ":" + str(self.charge) + ":" + str(self.precursor_mz) + ":" + str(self.kl_strict)

class DatabaseSpectrum(SpectrumBase):

    """Class for database results."""

    def __init__(self, scan_number, file_id, precursor_mz, charge, peptide, is_cluster):
        self.peptide = peptide
        self.is_cluster = is_cluster
        super(DatabaseSpectrum, self).__init__(scan_number, file_id, precursor_mz, charge)

class Purity(object):

    """Class for purity measures"""

    def __init__(self, representative_spectrum, purity_incl_undentified,
                 purity_no_undentified, peptide_counts, all_pur, sat):
        self.representative_spectrum = representative_spectrum
        self.purity_incl_undentified = purity_incl_undentified
        self.purity_no_undentified = purity_no_undentified
        self.all = all_pur
        self.sat = sat
        self.peptide_counts = peptide_counts

    def distinct_peptides(self):
        return len(peptide_counts.keys())

    def __str__(self):
        return "%.1f" % self.purity_incl_undentified + ": " + str(self.representative_spectrum)

##Core functions

def countby(func, items):
    sorted_items = sorted(items)
    return dict(
        (key,len(list(value)))
        for key, value in groupby(sorted_items, func)
    )

def safe_ratio(item1, item2):

    """
    Mainly for calculating purity,
    sometimes we want to ignore a ZeroDivisionError
    """

    ratio = 0
    try:
        ratio = (float(item1)/float(item2))
    except ZeroDivisionError:
        pass
    return ratio

def swap(dictionary):
    return dict(
        (dictionary[key],key)
        for key in dictionary
    )

def swap_many(dictionary):
    new_dict = dict(
        (dictionary[key],[])
        for key in dictionary
    )
    for key in dictionary:
        new_dict[dictionary[key]].append(key)
    return new_dict

def parse_xml_file(input_file):
    """From Ming, for params"""

    key_value_pairs = {}
    with open(input_file, 'r') as input_file_lines:
        for line in input_file_lines:
            new_line = line.rstrip().replace("<parameter name=\"", "")
            #new_line = new_line.replace("\">", "=")
            new_line = new_line.replace("</parameter>", "")

            splits = new_line.split("\">")
            #print splits
            if(len(splits) != 2):
                continue


            if(splits[0] in key_value_pairs.keys()):
              key_value_pairs[splits[0]].append(splits[1])
            else:
              key_value_pairs[splits[0]] = []
              key_value_pairs[splits[0]].append(splits[1])
    return key_value_pairs

def upload_file_mapping(params = None, params_xml = None):
    try:
        upload_file_mapping_list = params['upload_file_mapping']
    except:
        try:
            upload_file_mapping_list = parse_xml_file(params_xml)['upload_file_mapping']
        except:
            raise ValueError(str(params) + "Neither xml file or parameters provided")
    mapping = {}
    for item in upload_file_mapping_list:
        file_map = item.split("|")
        mapping[os.path.split(file_map[1])[1]] = file_map[0]
    return mapping

def datalist_to_map(datalist):
    with open(datalist, 'r') as datalist_items:
        return dict(
            (ind, os.path.split(filename)[1].replace("\n",""))
            for ind, filename in enumerate(datalist_items)
        )

#
# def id_sqs(database_spectra, spectrum, bucket_size):
#     output = 0
#     bucket = int(float(spectrum.precursor_mz)) - int(float(spectrum.precursor_mz)) % bucket_size
#     lower = []
#     try:
#         lower = database_spectra[bucket - bucket_size]
#     except:
#         lower = []
#     try:
#         db = [
#             database.sqs
#             for database in database_spectra[bucket] + lower
#             if database.scan_number == spectrum.scan_number.replace("out_0_0.","")
#         ]
#         try:
#             output = db[0]
#         except IndexError:
#             pass
#     except KeyError:
#         pass
#     return output
#
#
# def id_peptide(database_spectra, spectrum, bucket_size):
#     peptide = None
#     bucket = int(float(spectrum.precursor_mz)) - int(float(spectrum.precursor_mz)) % bucket_size
#     lower = []
#     try:
#         lower = database_spectra[bucket - bucket_size]
#     except:
#         lower = []
#     try:
#         db = [
#             database.peptide
#             for database in database_spectra[bucket] + lower
#             if database == spectrum
#         ]
#         try:
#             peptide = db[0]
#         except IndexError:
#             pass
#     except KeyError:
#         pass
#     return peptide
#
def purity_from_spectrum_list(spectra):
    only_id_spectra = [
        spectrum
        for spectrum in spectra
        if spectrum.peptide is not None
    ]
    total_spectra_incl_undentified = len(spectra)
    total_spectra_no_undentified = len(only_id_spectra)
    peptide_counts = countby(lambda x: x.peptide, only_id_spectra)
    grouped_spectra = list(countby(lambda x: x.peptide, only_id_spectra).items())
    representative_spectrum = ""
    number = 0.0
    if total_spectra_no_undentified != 0:
        representative_spectrum, number = (sorted(grouped_spectra, key = lambda x: x[1]))[-1]
    return Purity(
        representative_spectrum,
        safe_ratio(number, total_spectra_incl_undentified),
        0 if number == 0 else base.safe_ratio(number, total_spectra_no_undentified),
        peptide_counts
    )

def purity_from_peptide_list(peptides, all_peptides = [], sat_peptides = []):

    all_count = 0
    sat_count = 0

    id_peptides = [strip_peptide(pep) for pep in peptides if pep != 'PEPTIDE']
    if len(id_peptides) == 0:
        id_peptides = [strip_peptide(pep) for pep in all_peptides if pep != 'PEPTIDE']

    peptides_to_compare = [
        strip_peptide(peptide)
        for peptide in peptides
        if peptide is not None and peptide is not 'PEPTIDE' and peptide is not 'AMG' and peptide is not ''
    ]

    all_peptides_to_compare = [
        strip_peptide(peptide)
        for peptide in all_peptides
        if peptide is not None and peptide is not 'PEPTIDE' and peptide is not 'AMG' and peptide is not ''
    ]

    sat_peptides_to_compare = [
        strip_peptide(peptide)
        for peptide in sat_peptides
        if peptide is not None and peptide is not 'PEPTIDE' and peptide is not 'AMG' and peptide is not ''
    ]

    total_spectra_incl_undentified = len(peptides)
    total_spectra_no_undentified = len(peptides_to_compare)
    all_total_spectra_no_undentified = len(all_peptides_to_compare)

    peptides_to_calc = []
    original_mod = {}


    # peptide_counts = countby(lambda x: x, id_peptides)

    number = 0
    # for mixed_peptide in id_peptides:
    #     no_mod = mixed_peptide[0].split("!")
    #     with_mod = mixed_peptide[1].split("!")
    #     for i in range(0,len(no_mod)):
    #         peptides_to_calc.append((no_mod[i],with_mod[i]))
    #         original_mod[no_mod[i]] = with_mod[i]
    for mixed_peptide in id_peptides:
        peptides_to_calc.append(mixed_peptide)
        original_mod[mixed_peptide[0]] = mixed_peptide[1]
    if len(id_peptides) > 0:
        common_peptide = sorted(
                list(
                    countby(lambda x: x[0], peptides_to_calc).items()
                ),
                key = lambda x: x[1],
                reverse = True
            )
        try:
            most_common = common_peptide[0][0]
            representative_spectrum = original_mod[most_common]
        except:
            most_common = 'PEPTIDE'
            representative_spectrum = 'PEPTIDE'
    else:
        most_common = 'PEPTIDE'
        representative_spectrum = 'PEPTIDE'

    number_no_mod = 0

    cluster_set = set()

    for peptide in peptides_to_compare:
        if "!" in peptide[1]:
            for sub_peptide in peptide[1].split("!"):
                cluster_set.add(sub_peptide)

    if len(cluster_set) > 0:
        most_common = strip_peptide("!".join(list(cluster_set)))[0]
        representative_spectrum = strip_peptide("!".join(list(cluster_set)))[0]

    for peptide in peptides_to_compare:
        found = False
        for sub_peptide in peptide[1].split("!"):
            for sub_common in most_common.split("!"):
                if strip_peptide(sub_peptide)[0] == sub_common:
                    found = True
        if found:
            number_no_mod += 1


    for peptide in all_peptides_to_compare:
        found = False
        for sub_peptide in peptide[1].split("!"):
            for sub_common in most_common.split("!"):
                if strip_peptide(sub_peptide)[0] == sub_common:
                    found = True
        if found:
                all_count += 1

    for peptide in sat_peptides_to_compare:
        found = False
        for sub_peptide in peptide[1].split("!"):
            for sub_common in most_common.split("!"):
                if strip_peptide(sub_peptide)[0] == sub_common:
                    found = True
        if found:
                sat_count += 1

    return Purity(
        representative_spectrum,
        safe_ratio(number_no_mod, total_spectra_incl_undentified),
        0 if number_no_mod == 0 else safe_ratio(number_no_mod, total_spectra_no_undentified),
        dict(),
        0 if all_count == 0 else safe_ratio(all_count, len(all_peptides_to_compare)),
        0 if sat_count == 0 else safe_ratio(sat_count, len(sat_peptides_to_compare))
    )

def genereate_unambig_pep(peptide):
    pep_string_arr = []
    for pep in peptide.split("!"):
        try:
            pep_string_arr.append(''.join([p for p in pep.rsplit('.',1)[0].split('.',1)[1] if p.isalpha()]).replace('L','I'))
        except:
            pep_string_arr.append(''.join([p for p in pep if p.isalpha()]).replace('L','I'))
    return "!".join(pep_string_arr)


def strip_peptide(peptide):

    """Return only the letters in a peptide"""

    return (genereate_unambig_pep(peptide),peptide)

def unpack_peptide(tupled_tuple):
    peptide_search = tupled_tuple[0]
    count = tupled_tuple[1]
    return (peptide_search[1],count)

# def database_by_weight(database_spectra, bucket_size):
#     database = {}
#     for spectrum in database_spectra:
#         bucket = int(float(spectrum.precursor_mz)) - int(float(spectrum.precursor_mz)) % bucket_size
#         try:
#             database[bucket].append(spectrum)
#         except KeyError:
#             database[bucket] = [spectrum]
#     return database

# def read_tsv_database(filename, cluster = False, bucket_size = 1):
#     """Read the database output file -- not the clustering output."""
#     database_spectra = []
#     with open(filename) as tsv:
#         next(tsv)
#         for line in tsv:
#             split_line = line.split("\t")
#             scan = split_line[19]
#             if cluster:
#                 scan = int(split_line[47])
#             database_spectra.append(DatabaseSpectrum(
#                 scan_number=scan,
#                 file_id=split_line[29],
#                 precursor_mz=split_line[36],
#                 charge=split_line[9],
#                 peptide=split_line[0],
#                 is_cluster=cluster
#             ))
#     return database_by_weight(database_spectra, bucket_size)

def median_score(spectra):

    """Until ClusterLike"""

    weights = [
        float(spectrum.precursor_mz)
        for spectrum in spectra
    ]
    return round(statistics.median(weights),3)

def mix_score(spectra):

    """Until ClusterLike"""

    weights = [
        float(spectrum.precursor_mz)
        for spectrum in spectra
    ]
    max_weight = max(weights)
    min_weight = min(weights)
    return round(max_weight - min_weight,3)

def assign_sqs_from_spectra(clusters, cluster_sqs):
    for cluster in clusters:
        try:
            clusters[cluster].sqs = cluster_sqs[clusters[cluster].scan_number].sqs
        except:
            clusters[cluster].sqs = 0
    return clusters

def read_cluster_tsv(cluster_folder, database_spectra, title_prefix,
             database_clusters, data_map = None, bucket_size = 1, only_sqs = False, kl_map = {}):
    """Read the cluster output file -- not the database output."""
    clusters = {}
    for filename in glob.glob(cluster_folder + "/*.clust"):
        with open(filename) as tsv:
            current_cluster = None
            current_spectra = []
            line_number = 0
            for tsv_line in tsv:
                if tsv_line == '\n':
                    current_cluster.purity = 0 #purity_from_spectrum_list(current_spectra)
                    current_cluster.spectra = current_spectra
                    current_cluster.consensus_peptide = "" #id_peptide(database_clusters, current_cluster, bucket_size)
                    clusters.update({str(current_cluster.spectra[0].file_id) + ":" + str(current_cluster.spectra[0].scan_number): current_cluster})
                    current_spectra = []
                else:
                    split_tsv_line = tsv_line.replace('\n', '').split('\t')
                    if len(split_tsv_line) == 4:
                        current_cluster = read_cluster_tsv_line(split_tsv_line, title_prefix, filename)
                    elif len(split_tsv_line) >= 7:
                        spectrum = read_spectrum_tsv_line(split_tsv_line, database_spectra, data_map, only_sqs = only_sqs)
                        # spectrum.kl_strict = kl_map.get(str(spectrum.file_id) + ":" + str(spectrum.scan_number),0.0)
                        current_spectra.append(spectrum)
                    else:
                        raise ValueError(
                            'Line %d is neither a cluster nor a spectrum\n%s' % (line_number,
                                                                                 tsv_line)
                        )
                line_number += 1
    return clusters

def read_cluster_single(filename, database_spectra, title_prefix,
             database_clusters, data_map = None, bucket_size = 1, only_sqs = False, kl_map = {}):
    """Read the cluster output file -- not the database output."""
    clusters = {}
    with open(filename) as tsv:
        current_cluster = None
        current_spectra = []
        line_number = 0
        for tsv_line in tsv:
            if tsv_line == '\n':
                current_cluster.purity = 0 #purity_from_spectrum_list(current_spectra)
                current_cluster.spectra = current_spectra
                current_cluster.consensus_peptide = "" #id_peptide(database_clusters, current_cluster, bucket_size)
                clusters.update({str(current_cluster.spectra[0].file_id) + ":" + str(current_cluster.spectra[0].scan_number): current_cluster})
                current_spectra = []
            else:
                split_tsv_line = tsv_line.replace('\n', '').split('\t')
                if len(split_tsv_line) == 4:
                    current_cluster = read_cluster_tsv_line(split_tsv_line, title_prefix, filename)
                elif len(split_tsv_line) >= 7:
                    spectrum = read_spectrum_tsv_line(split_tsv_line, database_spectra, data_map, only_sqs = only_sqs)
                    # spectrum.kl_strict = kl_map.get(str(spectrum.file_id) + ":" + str(spectrum.scan_number),0.0)
                    current_spectra.append(spectrum)
                else:
                    raise ValueError(
                        'Line %d is neither a cluster nor a spectrum\n%s' % (line_number,
                                                                             tsv_line)
                    )
            line_number += 1
    return clusters

def read_cluster_tsv_line(tsv_line_array, title_prefix, filename):
    return Cluster(
        file_id=filename.replace(".clust",".mgf").replace("clust/",""),
        scan_number=tsv_line_array[0].replace(title_prefix,""),
        precursor_mz=tsv_line_array[2],
        charge=tsv_line_array[3],
        spectra=[],
        purity=None,
        consensus_peptide=None
    )

def read_spectrum_tsv_line(tsv_line_array, database_spectra, data_map, only_sqs, bucket_size = 1):
    try:
        file_id = data_map[int(tsv_line_array[1])]
    except:
        file_id = tsv_line_array[1]
    spectrum = Spectrum(
        scan_number=tsv_line_array[2],
        file_id=0 if only_sqs else file_id,
        precursor_mz=tsv_line_array[3],
        charge=tsv_line_array[6],
        similarity=tsv_line_array[4],
        p_value=tsv_line_array[5],
        sqs=0,
        peptide=None,
        cluster_member=tsv_line_array[0]
    )
    sqs = 0
    try:
        sqs = float(tsv_line_array[7])
    except IndexError:
        pass
    kl = 0
    try:
        kl = float(tsv_line_array[8])
    except IndexError:
        pass
    orig_pm = 0
    try:
        orig_pm = float(tsv_line_array[9])
    except IndexError:
        pass
    cosine = 0
    try:
        cosine = float(tsv_line_array[10])
    except IndexError:
        pass
    signal_peaks = 0
    try:
        signal_peaks = float(tsv_line_array[11])
    except IndexError:
        pass
    order = 0
    try:
        signal_peaks = int(tsv_line_array[12])
    except IndexError:
        pass
    spectrum.sqs = sqs
    spectrum.kl_strict = kl
    spectrum.orig_pm = orig_pm
    spectrum.cosine = cosine
    spectrum.signal_peaks = signal_peaks
    spectrum.peptide = None
    spectrum.order = order
    return spectrum

def as_spectra(cluster_folder, database_spectra, title_prefix,
         database_clusters, bucket_size = 1):

    clusters = read_cluster_tsv(cluster_folder, database_spectra, title_prefix,
             database_clusters, bucket_size, only_sqs = True)



    flat_clusters = [
        clusters[key]
        for key in clusters
    ]

    cluster_spectra = dict(
        (spectrum.scan_number, spectrum)
        for cluster in flat_clusters
        for spectrum in cluster.spectra
    )

    return cluster_spectra

# def init_from_pickle(cluster_pickle):
#     clusters = []
#     with open(cluster_pickle,"rb") as f:
#         clusters = pickle.load(f)
#     print("Successfully loaded " + str(len(clusters)) + " clusters from " + cluster_pickle)
#     return cluster_session.ClusterSession(clusters)

def initialize_unidentified_clusters(cluster_folder, clustered_clusters_folder, title_prefix, data_map, kl_map):
    # clusters_as_spectra = as_spectra(clustered_clusters_folder, {}, title_prefix, {})
    clusters = read_cluster_tsv(cluster_folder, {}, title_prefix, {}, data_map, kl_map)
    return assign_sqs_from_spectra(clusters, clusters_as_spectra)

def find_cluster(clusters, cluster_id):
    found_cluster = None
    for cluster in clusters:
        if int(cluster.scan_number) == cluster_id:
            found_cluster = cluster
            break
    return found_cluster

def initialize_from_proteosafe(clusters_file, spectra_file, parameter_file, kl_file):
    clusters = {}
    file_name_map = swap(upload_file_mapping(params_xml=parameter_file))
    with open(clusters_file, 'r') as cluster_lines:
        cluster_lines.readline()
        for cluster_line in cluster_lines:
            clust_split = cluster_line.split('\t')
            scan_number = clust_split[2]
            charge  = 0
            try:
                charge = int(clust_split[6].replace('+',''))
            except:
                charge = int(clust_split[6].replace('+','').split('_')[1])
            new_cluster = Cluster(
                scan_number = scan_number,
                file_id = clust_split[0],
                precursor_mz = clust_split[4],
                charge = charge,
                spectra = [],
                purity = None,
                consensus_peptide = clust_split[7],
                ms_cluster_charge = clust_split[21],
                spec_prob = clust_split[11],
                p_value = float(clust_split[12]),
                fdr = float(clust_split[13]),
                pep_fdr = float(clust_split[14]),
                sqs = float(clust_split[18])
            )
            clusters[scan_number] = new_cluster
    strict = {}
    non_strict = {}
    for filename in glob.glob(kl_file + "/*.csv"):
        with open(filename) as f:
            for line in f:
                split_line = line.split(",")
                scan_file = os.path.basename(split_line[0]).split(".")[0]
                try:
                    scan_number = split_line[1].replace(" ","")
                except:
                    print(line)
                try:
                    kl = float(split_line[6])
                except:
                    print(line)
                if "_t" in filename:
                    strict[scan_file + ":" + scan_number] = kl
                else:
                    non_strict[scan_file + ":" + scan_number] = kl
    with open(spectra_file, 'r') as spectrum_lines:
        spectrum_lines.readline()
        for spectrum_line in spectrum_lines:
            spec_split = spectrum_line.split('\t')
            cluster_number = spec_split[15]
            scan_number = spec_split[2]
            file_id = file_name_map[spec_split[0]]
            try:
                charge = int(spec_split[6])
            except:
                charge = int(spec_split[6].split('_')[1])
            kl_strict = strict.get(file_id.split(".")[0] + ":" + scan_number)
            kl_non_strict = non_strict.get(file_id.split(".")[0] + ":" + scan_number)
            kl_delta = None
            if kl_strict and kl_non_strict:
                kl_delta = math.log(kl_strict/kl_non_strict);
            new_spectrum = Spectrum(
                scan_number = scan_number,
                file_id = file_id,
                precursor_mz = spec_split[4],
                charge = charge,
                similarity = None,
                p_value = float(spec_split[12]),
                sqs = float(spec_split[16]),
                peptide = spec_split[7],
                cluster_number = cluster_number,
                spec_prob = float(spec_split[11]),
                fdr = float(spec_split[13]),
                pep_fdr = float(spec_split[14]),
                kl_strict = kl_strict,
                kl_non_strict = kl_non_strict,
                kl_delta = kl_delta
                )
            clusters[cluster_number].spectra.append(new_spectrum)
    for cluster in clusters:
        peptides = [
            spectrum.peptide
            for spectrum in clusters[cluster].spectra
        ]
        clusters[cluster].purity = purity_from_peptide_list(peptides)
    # all_peptides = [
    #     clusters[cluster].consensus_peptide
    #     for cluster in clusters
    # ]
    # grouped_and_sorted = {
    #     key: len(list(value))
    #     for key, value in groupby(sorted(all_peptides),key=lambda x: x)
    #     if key != ''
    # }
    # for cluster in clusters:
    #     try:
    #         clusters[cluster].of_split = grouped_and_sorted[clusters[cluster].consensus_peptide]
    #     except:
    #         clusters[cluster].of_split = 0
    return clusters

def initialize_from_proteosafe_file(clusters_file, spectra_file, parameter_file, kl_file):
    clusters = {}
    file_name_map = swap(upload_file_mapping(params_xml=parameter_file))
    with open(clusters_file, 'r') as cluster_lines:
        cluster_lines.readline()
        for cluster_line in cluster_lines:
            clust_split = cluster_line.split('\t')
            scan_number = clust_split[2]
            charge  = 0
            try:
                charge = int(clust_split[6].replace('+',''))
            except:
                charge = int(clust_split[6].replace('+','').split('_')[1])
            new_cluster = Cluster(
                scan_number = scan_number,
                file_id = clust_split[0],
                precursor_mz = clust_split[4],
                charge = charge,
                spectra = [],
                purity = None,
                consensus_peptide = clust_split[7],
                ms_cluster_charge = clust_split[21],
                spec_prob = clust_split[11],
                p_value = float(clust_split[12]),
                fdr = float(clust_split[13]),
                pep_fdr = float(clust_split[14]),
                sqs = float(clust_split[18])
            )
            clusters[scan_number] = new_cluster
    strict = {}
    non_strict = {}
    for filename in glob.glob(kl_file + "/*.csv"):
        with open(filename) as f:
            for line in f:
                split_line = line.split(",")
                scan_file = os.path.basename(split_line[0]).split(".")[0]
                try:
                    scan_number = split_line[1].replace(" ","")
                except:
                    print(line)
                try:
                    kl = float(split_line[6])
                except:
                    print(line)
                if "_t" in filename:
                    strict[scan_file + ":" + scan_number] = kl
                else:
                    non_strict[scan_file + ":" + scan_number] = kl
    with open(spectra_file, 'r') as spectrum_lines:
        spectrum_lines.readline()
        for spectrum_line in spectrum_lines:
            spec_split = spectrum_line.split('\t')
            cluster_number = spec_split[15]
            scan_number = spec_split[2]
            file_id = file_name_map[spec_split[0]]
            try:
                charge = int(spec_split[6])
            except:
                charge = int(spec_split[6].split('_')[1])
            kl_strict = strict.get(file_id.split(".")[0] + ":" + scan_number)
            kl_non_strict = non_strict.get(file_id.split(".")[0] + ":" + scan_number)
            kl_delta = None
            if kl_strict and kl_non_strict:
                kl_delta = math.log(kl_strict/kl_non_strict);
            new_spectrum = Spectrum(
                scan_number = scan_number,
                file_id = file_id,
                precursor_mz = spec_split[4],
                charge = charge,
                similarity = None,
                p_value = float(spec_split[12]),
                sqs = float(spec_split[16]),
                peptide = spec_split[7],
                cluster_number = cluster_number,
                spec_prob = float(spec_split[11]),
                fdr = float(spec_split[13]),
                pep_fdr = float(spec_split[14]),
                kl_strict = kl_strict,
                kl_non_strict = kl_non_strict,
                kl_delta = kl_delta
                )
            clusters[cluster_number].spectra.append(new_spectrum)
    for cluster in clusters:
        peptides = [
            spectrum.peptide
            for spectrum in clusters[cluster].spectra
        ]
        clusters[cluster].purity = purity_from_peptide_list(peptides)
    # all_peptides = [
    #     clusters[cluster].consensus_peptide
    #     for cluster in clusters
    # ]
    # grouped_and_sorted = {
    #     key: len(list(value))
    #     for key, value in groupby(sorted(all_peptides),key=lambda x: x)
    #     if key != ''
    # }
    # for cluster in clusters:
    #     try:
    #         clusters[cluster].of_split = grouped_and_sorted[clusters[cluster].consensus_peptide]
    #     except:
    #         clusters[cluster].of_split = 0
    return clusters
#
# def load_to_pickle(cluster_folder, clusters_clusters_folder, identified_clusters,
#                    identified_spectra, output_pickle, title_prefix):
#     print("Loading clusters...this might take a moment.")
#     database_spectra = database.read_tsv(identified_spectra)
#     database_clusters = database.read_tsv(identified_clusters, cluster = True)
#     clusters_as_spectra = as_spectra(clusters_clusters_folder, database_spectra, title_prefix, database_clusters)
#     clusters = read_tsv(cluster_folder, database_spectra, title_prefix, database_clusters)
#
#     new_clusters = assign_sqs_from_spectra(clusters, clusters_as_spectra)
#
#     with open(output_pickle,"wb") as f:
#         pickle.dump(new_clusters, f)
#
#     print("Saved " + str(len(clusters)) + " clusters to " + output_pickle)

def calculate_precursor_mix(cluster,tolerance):
    mix = 0
    spectra = cluster.spectra
    tolerance_percent = float(tolerance)/1000000
    charge = 0
    pep_rep = cluster.purity.representative_spectrum
    for spectrum in spectra:
        try:
            spec_id = spectrum.file_id + ":" + spectrum.scan_number
            if spectrum.peptide == pep_rep:
                charge = spectrum.charge
        except:
            pass
        if charge != 0:
            break
    monoisotopic_mass = 0
    if pep_rep != 'PEPTIDE' and pep_rep !='':
        rep_isotopes = isotopes(pep_rep,int(charge))
        monoisotopic_mass = rep_isotopes[0]
        similar = 0
        total = len(spectra)
        for spectrum in spectra:
            for isotope in rep_isotopes:
                # print(float(spectrum.precursor_mz))
                # print(isotope)
                if (float(spectrum.precursor_mz) > (isotope - tolerance_percent*isotope) and float(spectrum.precursor_mz) < (isotope + tolerance_percent*isotope)):
                    similar += 1
                    break
        mix = float(similar)/float(total)
    return mix,monoisotopic_mass
