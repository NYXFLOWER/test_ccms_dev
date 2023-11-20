import glob
from csv import DictReader,DictWriter
from collections import Counter, defaultdict, namedtuple
import pyteomics
from pyteomics import mzxml, mzml, mgf
import ming_psm_library as mpl
import ming_spectrum_library as msl
import ming_numerical_utilities as mnu
import sys

INPUT_JOBS = sys.argv[1]
OUTPUT_FOLDER = sys.argv[2]
FILENAME_HEADER = sys.argv[3]
SCAN_HEADER = sys.argv[4]
PEPTIDE1_HEADER = sys.argv[5]
PEPTIDE2_HEADER = sys.argv[6]
CHARGE1_HEADER = sys.argv[7]
CHARGE2_HEADER = sys.argv[8]
TYPE_HEADER = sys.argv[9]
INTENSITY_HEADER = sys.argv[10]
SNR_FILTER = int(sys.argv[11])
LOW_MASS_FILTER = int(sys.argv[12])
FILTER_PRECURSOR = sys.argv[13] == 'Yes'
PEAK_TOLERANCE = float(sys.argv[14])
USE_MGF_INDEX = sys.argv[15] == 'Yes'

def peak_intersection(peaks1, peaks2, tol):
    common_peaks = []
    for mz1,_ in peaks1:
        found = False
        for mz2,_ in peaks2:
            if abs(mz1-mz2) <= tol:
                found = True
                common_peaks.append(mz2)
        if found:
            common_peaks.append(mz1)
    return set(common_peaks)

def normalize_peptide(peptide):
    if peptide and peptide != '--':
        if len(peptide.split(".")) > 1 and len(peptide.split(".")[0]) == 1 and len(peptide.split(".")[-1]) == 1:
            peptide = ".".join(peptide.split(".")[1:-1])
        return peptide.split("}")[-1].replace(",-","-").replace(",","+").replace("[-","-").replace("[","+").replace("(","").replace(")","").replace("]","")
    else:
        return 'PEPTIDE'

def determine_peptides_if_mix(peptide1, charge1, peptide2, charge2):
    split_peptide1 = peptide1.split('!')
    split_charge1 = charge1.split('!')
    if not peptide2 or peptide2 == '--' or '!' in peptide2:
        peptide1 = split_peptide1[0]
        peptide2 = split_peptide1[1]
        charge1 = split_charge1[0]
        charge2 = split_charge1[1]
    else:
        if normalize_peptide(peptide1.split('!')[0]) == normalize_peptide(peptide2):
            peptide1 = split_peptide1[1]
            peptide2 = peptide2
            charge1 = split_charge1[1]
            charge2 = charge2
        else:
            peptide1 = split_peptide1[0]
            peptide2 = peptide2
            charge1 = split_charge1[0]
            charge2 = charge2
    return peptide1, charge1, peptide2, charge2

def determine_peptides(peptide1, charge1, peptide2, charge2):
    if peptide1 and '!' in peptide1:
        peptide1, charge1, peptide2, charge2 = determine_peptides_if_mix(peptide1, charge1, peptide2, charge2)
    elif peptide2 and '!' in peptide2:
        peptide2, charge2, peptide1, charge1 = determine_peptides_if_mix(peptide2, charge2, peptide1, charge1)
    int_charge1 = int(charge1) if charge1 else None
    int_charge2 = int(charge2) if charge2 else None
    return (peptide1, int_charge1), (peptide2, int_charge2)

def create_intensity_output(s, intensities, psms_per_scan, scan, precursor_func):
    mz = s['m/z array']
    intn = s['intensity array']
    precursor = precursor_func(s)
    peaks = [[mz[i],intn[i]] for i in range(len(mz))]
    unfiltered_peaks = peaks
    peaks = filter_peaks(peaks, precursor)
    for spec in psms_per_scan:
        cluster_type = spec[0]
        snr_annotations = []
        snr_signal = []
        annotations = []
        original_peptide = spec[1]
        for peptide in spec[2:]:
            charge = peptide[1]
            peptide = normalize_peptide(peptide[0])
            if peptide and peptide != 'PEPTIDE':
                annotated_peaks, annotated_peaks_at_snr, signal_peaks_at_snr = scan_peptide_to_intensity2(peaks,unfiltered_peaks,charge,peptide)
                snr_annotations.append(annotated_peaks_at_snr)
                snr_signal.append(signal_peaks_at_snr)
                annotations.append(annotated_peaks)

        total_intensity = sum(peak[1] for peak in peaks)

        if len(annotations) > 0 and total_intensity > 0:

            explained_intensity_unique1 = 0
            explained_intensity_unique2 = 0
            explained_intensity_shared = 0


            explained_intensity1 = sum(peak[1] for peak in annotations[0])/total_intensity

            if len(annotations) > 1:
                shared_peaks = peak_intersection(annotations[0], annotations[1], PEAK_TOLERANCE)
                explained_intensity_unique1 = sum(peak[1] for peak in annotations[0] if peak[0] not in shared_peaks)/total_intensity
                explained_intensity_unique2 = sum(peak[1] for peak in annotations[1] if peak[0] not in shared_peaks)/total_intensity
                explained_intensity_shared = sum(peak[1] for peak in annotations[1] if peak[0] in shared_peaks)/total_intensity

            intensities[(str(scan), str(cluster_type), str(original_peptide))] = (explained_intensity1, explained_intensity_unique1, explained_intensity_unique2, explained_intensity_shared, snr_annotations[0], snr_signal[0])

    return intensities

def filter_peaks(peaks, precursor):
    peaks = [p for p in peaks if p[0] >= LOW_MASS_FILTER]
    if FILTER_PRECURSOR:
        peaks = msl.filter_precursor_peaks(peaks,PEAK_TOLERANCE,precursor)
    #apply SNR filter
    peaks = msl.filter_noise_peaks(peaks, SNR_FILTER)
    return peaks

def explained_intensity2(input_file,use_scans,spec_dict):
    intensities = {}
    input_filetype = input_file.split('.')[-1]
    if input_filetype == 'mzXML':
        precursor_func = lambda spectrum: float(spectrum["precursorMz"][0]["precursorMz"])
        with mzxml.read(input_file) as reader:
            for s in reader:
                scan = str(s['id'])
                if scan in spec_dict and s['msLevel'] == '2' or s['msLevel'] == 2:
                    intensities = create_intensity_output(s, intensities, spec_dict[scan], scan, precursor_func)
    elif input_filetype == 'mzML':
        precursor_func = lambda spectrum: float(s['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z'])
        with mzml.read(input_file) as reader:
            for s in reader:
                scan = str(s['id'].split('scan=')[1])
                if scan in spec_dict and s['ms level'] == '2' or s['ms level'] == 2:
                    intensities = create_intensity_output(s, intensities, spec_dict[scan], scan, precursor_func)
    elif input_filetype == 'mgf':
        precursor_func = lambda spectrum: float(s['params']['pepmass'][0])
        with mgf.read(input_file) as reader:
            for i,s in enumerate(reader):
                if USE_MGF_INDEX:
                    scan = str(i+1)
                else:
                    scan = str(s['params']['scans'])
                if scan in spec_dict:
                    intensities = create_intensity_output(s, intensities, spec_dict[scan], scan, precursor_func)

    return intensities

def scan_peptide_to_intensity2(peaks,unfiltered_peaks,charge,peptide):
    # explained_intensity = mpl.calculated_explained_intensity(
    #     spectrum['peaks'],                                                   # MS2 peaks
    #     int(spectrum['precursorMz'][0]['precursorCharge'])-1,                # highest charge ions to look for, eg. 3 means 1+ ions and 2+ ions
    #     peptide,                                                          # Peptide string in InSPECT format N+1, not (N,1)
    #     PEAK_TOLERANCE                                                             # tolerance
    # )
    _,annotated_peaks = mpl.calculated_number_annotated_peaks(
        peaks,                                                   # MS2 peaks
        max(charge-1,1),                # highest charge ions to look for, eg. 3 means 1+ ions and 2+ ions
        peptide,                                                          # Peptide string in InSPECT format N+1, not (N,1)
        PEAK_TOLERANCE                                                             # tolerance
    )
    annotated_peaks_at_snr = [0,0,0,0,0]
    signal_peaks_at_snr = [0,0,0,0,0]
    for snr in range(len(annotated_peaks_at_snr)):
        peaks_at_snr = mnu.return_signal_peaks(unfiltered_peaks,snr+1)
        if peaks_at_snr != 0.0:
            annotated_peaks_at_snr[snr],_ = mpl.calculated_number_annotated_peaks(
                peaks_at_snr,                                                   # MS2 peaks
                max(charge-1,1),                # highest charge ions to look for, eg. 3 means 1+ ions and 2+ ions
                peptide,                                                          # Peptide string in InSPECT format N+1, not (N,1)
                PEAK_TOLERANCE                                                             # tolerance
            )
            signal_peaks_at_snr[snr] = len(peaks_at_snr)

    return annotated_peaks, annotated_peaks_at_snr, signal_peaks_at_snr


for job in open(INPUT_JOBS):
    split_job = job.replace("\n","").split("\t")
    id_file = split_job[0]
    spec_file = split_job[1]
    spec_file2 = split_job[2]
    spec_dict = defaultdict(list)
    header = {}
    lines = []
    with open(id_file) as f:
        header = f.readline().replace("\n","").split("\t")
    header.extend([INTENSITY_HEADER,"% Unique Intensity PSM1", "% Unique Intensity PSM2", "% Shared Intensity"] + [h for i in range(5) for h in ["Annotated Peaks at SNR{}".format(i+1),"Signal Peaks at SNR{}".format(i+1)] ])
    if "InternalFilename" not in header:
        header.append("InternalFilename")
    if "InternalFilename2" not in header:
        header.append("InternalFilename2")
    if "Scan" not in header:
        header.append("Scan")
    if "Peptide1" not in header:
        header.append("Peptide1")
    if "Peptide2" not in header:
        header.append("Peptide2")
    if "Charge1" not in header:
        header.append("Charge1")
    if "Charge2" not in header:
        header.append("Charge2")
    with open(id_file) as f:
        r = DictReader(f, delimiter = "\t")
        for l in r:
            if TYPE_HEADER == "":
                cluster_type = 'C'
            else:
                try:
                    cluster_type = l[TYPE_HEADER]
                except:
                    cluster_type = 'C'

            (peptide1,charge1), (peptide2,charge2) = determine_peptides(l.get(PEPTIDE1_HEADER),l.get(CHARGE1_HEADER),l.get(PEPTIDE2_HEADER),l.get(CHARGE2_HEADER))
            spec_dict[l[SCAN_HEADER]].append([cluster_type,l.get(PEPTIDE1_HEADER),(peptide1,charge1),(peptide2,charge2)])

            lines.append(l)
    intensities = explained_intensity2(spec_file,False,spec_dict)
#    print(spec_dict)
#    print(intensities)

    with open(OUTPUT_FOLDER + "/" + id_file.split("/")[-1], "w") as w:
        wr = DictWriter(w, header, delimiter = "\t")
        wr.writeheader()
        for line in lines:
            if TYPE_HEADER == "":
                cluster_type = 'C'
            else:
                try:
                    cluster_type = line[TYPE_HEADER]
                except:
                    cluster_type = 'C'
            # spec_dict[l[SCAN_HEADER]] = [cluster_type,line[PEPTIDE_HEADER]]
            line["InternalFilename"] = spec_file
            line["InternalFilename2"] = spec_file2
            line["Scan"] = line[SCAN_HEADER]
            (peptide1,charge1), (peptide2,charge2) = determine_peptides(line.get(PEPTIDE1_HEADER),line.get(CHARGE1_HEADER),line.get(PEPTIDE2_HEADER),line.get(CHARGE2_HEADER))
            line["Peptide1"] = normalize_peptide(peptide1)
            line["Peptide2"] = normalize_peptide(peptide2)
            line["Charge1"] = charge1
            line["Charge2"] = charge2
            if (line[SCAN_HEADER],cluster_type,line[PEPTIDE1_HEADER]) in intensities:
                explained_intensity1, explained_intensity_unique1, explained_intensity_unique2, explained_intensity_shared, snr_annotations, snr_signal = intensities[(line[SCAN_HEADER],cluster_type,line[PEPTIDE1_HEADER])]
                combined_intensity = explained_intensity_shared + explained_intensity_unique1 + explained_intensity_unique2
                line[INTENSITY_HEADER] = combined_intensity if combined_intensity > 0 else explained_intensity1
                line["% Unique Intensity PSM1"] = "{0:.2f}".format(explained_intensity_unique1)
                line["% Unique Intensity PSM2"] = "{0:.2f}".format(explained_intensity_unique2)
                line["% Shared Intensity"] = "{0:.2f}".format(explained_intensity_shared)
                for i in range(len(snr_annotations)):
                    line["Annotated Peaks at SNR{}".format(i+1)] = snr_annotations[i]
                    line["Signal Peaks at SNR{}".format(i+1)] = snr_signal[i]
            else:
                line[INTENSITY_HEADER] = 0
                line["% Unique Intensity PSM1"] = 0
                line["% Unique Intensity PSM2"] = 0
                line["% Shared Intensity"] = 0
                for i in range(5):
                    line["Annotated Peaks at SNR{}".format(i+1)] = 0
                    line["Signal Peaks at SNR{}".format(i+1)] = 0
            wr.writerow(line)
        # w.write("\t".join([
        #             "Filename",
        #             "Scan",
        #             "Type",
        #             "Peptide",
        #             "ExplainedIntensity"
        #         ]) + "\n")
        #
        # for row in intensities:
        #     w.write("\t".join([id_file.split("/")[-1].split(".")[0]] + [str(r) for r in row]) + "\n")
