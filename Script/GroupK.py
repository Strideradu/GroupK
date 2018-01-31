# use minimap as filter for the experiment
import subprocess
from Bio import SeqIO
from utils import *
import os
from groupsChain import *
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("input", help="path of input fasta file", type=str)
parser.add_argument("--k1", help="kmer size for filtration", type=int, default=15)
parser.add_argument("--threshold", help="count threshold for shared k1", type=int, default=2)
parser.add_argument("--k2", help="kmer size for group", type=int, default=9)
parser.add_argument("--accuracy", help="accuracy of the reads", type=float, default=0.85)
parser.add_argument("--gap", help="gap rate", type=float, default=0.12)
parser.add_argument("--chain", help="number of kmer in the chain to report", type=int, default=3)
parser.add_argument("--groupc", help="the coefficeint c in paper to control chaining threshold, must be non-zero float",
                    type=float, default=4)
parser.add_argument("--idbase",
                    help="if the number of identical base meet this threshold the output will always be reported",
                    type=int, default=400)
parser.add_argument("--ratio", help="ratio of two overlap region threshold", type=float, default=0.5)
parser.add_argument("--size", help="group size threshold", type=int, default=12)
parser.add_argument("--large", help="release threhold for large overlap", type=float, default=1.5)

try:
    args = parser.parse_args()

except:
    parser.print_help()
    sys.exit(1)

save_path = os.path.dirname(args.input)
fasta_path = args.input
k = args.k1
threshold = args.threshold

filter_output = os.path.join(save_path, "sa_k{}_c{}_filter.out".format(k, threshold))

filter_command = filter_path + " -i {} -k {} -o {}".format(fasta_path, k, filter_output)
filter_process = subprocess.Popen(filter_command, stdout=subprocess.PIPE, shell=True)
filter_process.wait()

# parse minimap output to generate dict to set
filter_dict = {}
with open(filter_output) as f:
    for line in f:
        line = line.rstrip()
        if line != "" and line[0] != "#":
            line_sp = line.split("\t")
            query_id = line_sp[0]
            target_id = line_sp[1]
            count = int(line_sp[2])
            if count >= threshold:
                if filter_dict.get(query_id, False):
                    filter_dict[query_id].add(target_id)
                else:
                    filter_dict[query_id] = set([target_id])

group_output = os.path.join(save_path, "groups_k{}_c{}.out".format(k, threshold))
record_dict = SeqIO.index(fasta_path, "fasta")
with open(group_output, "w") as fout:
    for query_id in filter_dict.keys():
        record = record_dict[query_id]
        query_save = save_path + "query_k{}_c{}.fasta".format(k, threshold)
        record.description = ""
        SeqIO.write(record, query_save, "fasta")

        target_output = []
        target_list = list(filter_dict[query_id])
        for target_id in target_list:
            record = record_dict[target_id]
            target_output.append(record)
        target_save = save_path + "target_k{}_c{}.fasta".format(k, threshold)
        SeqIO.write(target_output, target_save, "fasta")

        temp_output = save_path + "groups_temp_k{}_c{}.out".format(k, threshold)
        yass_command = yass_path + ' -p "{}" -m {} -i {} -o {} {} {}'.format('#' * args.k2, args.accuracy, args.gap,
                                                                             temp_output, query_save, target_save)
        process = subprocess.Popen(yass_command, stdout=subprocess.PIPE, shell=True)
        process.wait()

        with open(temp_output) as f:
            for line in f.readlines():
                print(line, end="", file=fout)


L = statistical_bound_of_waiting_time(args.accuracy, args.k2)

with open(group_output) as f:
    lines = f.readlines()

    for line in lines:
        group_hit = GroupHit(line, args.size)
        group_hit.set_ratio(args.ratio)

        # print group_hit.groups

        group_hit.chain_groups(accuracy=args.accuracy, group_distance=L, rechain_threshold=args.chain,
                               span_coefficient=args.groupc, identity=args.idbase, release= args.large)
        # print group_hit.chain_align
        if group_hit.aligned:
            output_str = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                group_hit.query, group_hit.target, str(group_hit.aligned), group_hit.aligned_base,
                group_hit.overlap_length, group_hit.query_ali_start, group_hit.query_ali_end,
                group_hit.target_ali_start, group_hit.target_ali_end, group_hit.query_overlap_start,
                group_hit.query_overlap_end, group_hit.target_overlap_start, group_hit.target_overlap_end
            )

            print(output_str)
