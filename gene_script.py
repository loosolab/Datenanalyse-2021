import os
import sys

import pandas as pd


def overlaps(a, b):
    #  based on: https://stackoverflow.com/questions/2953967/built-in-function-for-computing-overlap-in-python
    return min(a[1], b[1]) - max(a[0], b[0])


def get_gencode(path):
    table = pd.read_csv(path, sep=" ", names=["Chr", "Start", "End", "Gene"])

    return table.values.tolist()


def get_cluster(path):
    table = pd.read_csv(path, sep="\t", names=["Chr", "Start", "End", "Name", "Score", "Strand", "Signal Value",
                                               "p value", "q value", "Peak"])
    table.drop(["Name", "Score", "Strand", "p value", "q value", "Peak"], axis=1, inplace=True)
    table.sort_values(by="Signal Value", ascending=False, inplace=True)
    table = table.head(100)

    return table.values.tolist()


def annotate_genes(top_peaks, gene_table):
    merged_lines = []
    for line in top_peaks:
        new_line = line
        for gene in gene_table:
            if line[0] == gene[0]:
                if overlaps([line[1], line[2]], [gene[1], gene[2]]) > 0:
                    if len(new_line) == 4:
                        new_line.append(gene[3])
                    elif gene[3] not in new_line[4]:
                        new_line[4] = new_line[4] + ", " + gene[3]

        merged_lines.append(new_line)

    return merged_lines


def write_annotated_cluster(lines, file_name):
    with open(file_name, "w") as cluster_file:
        for line in lines:
            for col in line:
                cluster_file.write(str(col) + "\t")
            cluster_file.write("\n")


def write_clusters():
    cluster_folder = sys.argv[1]
    cluster_files = os.listdir(cluster_folder)
    gene_file = sys.argv[2]
    print("Reading " + gene_file)
    gencode = get_gencode(gene_file)
    for file in cluster_files:
        cluster_peaks = get_cluster(cluster_folder + "/" + file)
        write_annotated_cluster(annotate_genes(cluster_peaks, gencode), file.split("_")[0] + ".peaks.annotated.bed")
        print(file + " annotated")


def main():
    write_clusters()
    exit(0)


if __name__ == '__main__':
    main()
