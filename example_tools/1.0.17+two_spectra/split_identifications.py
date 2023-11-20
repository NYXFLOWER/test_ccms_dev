import glob
import sys
from itertools import groupby
from csv import DictReader, DictWriter
import cluster_base
import math
from collections import defaultdict

INPUT_IDS = sys.argv[1]
PARAMS = sys.argv[2]
OUTPUT_FOLDER = sys.argv[3]
SPECTRUM_FOLDER = sys.argv[4]
OUTPUT_JOBS = sys.argv[5]
FILENAME_HEADER = sys.argv[6]
PARALLEL_FACTOR = int(sys.argv[7])
SPECTRUM_FOLDER2 = sys.argv[8]
FILENAME_HEADER2 = sys.argv[9]

all_files = []
all_lines = []
header = []

if PARAMS != 'NONE':
    parameters = cluster_base.parse_xml_file(PARAMS)
    param_dict = cluster_base.upload_file_mapping(params = parameters)
else:
    param_dict = {}

with open(INPUT_IDS) as f:
    header = f.readline().replace("\n","").split("\t")

with open(INPUT_IDS) as f:
    r = DictReader(f, delimiter = "\t")
    for l in r:
        l['partition'] = 0
        all_lines.append(l)

file_count = len(set([(x[FILENAME_HEADER],x['partition']) for x in all_lines]))
grouped_lines = groupby(sorted(all_lines,key=lambda x: (x[FILENAME_HEADER], x[FILENAME_HEADER2], x['partition'])), key=lambda x: (x[FILENAME_HEADER], x[FILENAME_HEADER2], x['partition']))

if file_count < PARALLEL_FACTOR:
    max_jobs_per_file = math.floor(PARALLEL_FACTOR/file_count)
    job_per_file = defaultdict(int)
    for line in all_lines:
        current = job_per_file[line[FILENAME_HEADER]]
        line['partition'] = current
        job_per_file[line[FILENAME_HEADER]] = (current+1)%max_jobs_per_file
    grouped_lines = groupby(sorted(all_lines,key=lambda x: (x[FILENAME_HEADER], x[FILENAME_HEADER2], x['partition'])), key=lambda x: (x[FILENAME_HEADER], x[FILENAME_HEADER2], x['partition']))

print(param_dict)

for ((filename_full,filename_full2,partition),lines) in grouped_lines:
    filename = filename_full.replace("f.","").replace("d.","").split("/")[-1]
    filename2 = filename_full2.replace("f.","").replace("d.","").split("/")[-1]
    print(filename)
    if filename in param_dict or PARAMS == 'NONE':
        partitioned_filename = "{}_{}.tsv".format(filename.split(".")[0], partition)
        print(partitioned_filename)
        all_files.append((partitioned_filename,filename,filename2))
        with open(OUTPUT_FOLDER + "/" + partitioned_filename, "w") as w:
            wr = DictWriter(w, header, delimiter = "\t")
            wr.writeheader()
            for line in lines:
                line.pop('partition')
                wr.writerow(line)

files_per_split = max(int(len(all_files)/PARALLEL_FACTOR),1)
#
# if len(all_files) < PARALLEL_FACTOR:
#     files_per_split

print(files_per_split)

file_number = 1

current_file = open(OUTPUT_JOBS + "/" + str(file_number), 'w')

print("Filename header:", FILENAME_HEADER)
print("Filename header2:", FILENAME_HEADER2)

for (i,(filename, original_filename, original_filename2)) in enumerate(all_files):
        print("current file:", current_file)
        print("line data:", filename, original_filename, original_filename2)
        if i % (files_per_split) == 0 and i != 0:
            current_file.close()
            file_number += 1
            current_file = open(OUTPUT_JOBS + "/" + str(file_number), 'w')
        demanged_filename = param_dict.get(original_filename)
        if not demanged_filename:
            if PARAMS == 'NONE':
                demanged_filename = original_filename
            else:
                raise Exception
        demanged_filename2 = param_dict.get(original_filename2)
        if not demanged_filename2:
            if PARAMS == 'NONE':
                demanged_filename2 = original_filename2
            else:
                raise Exception
        current_file.write(OUTPUT_FOLDER + "/" + filename + "\t" + SPECTRUM_FOLDER + "/" + demanged_filename + "\t" + SPECTRUM_FOLDER2 + "/" + demanged_filename2 + "\n")
