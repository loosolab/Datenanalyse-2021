import pandas as pd


def overlaps(a, b):
    #  based on: https://stackoverflow.com/questions/2953967/built-in-function-for-computing-overlap-in-python
    return min(a[1], b[1]) - max(a[0], b[0])


def get_genecode():
    path = "/home/mfk/Google Drive/JLU/1. Semester/Datenanalyse/Gene_script/gencode.filtered.bed"
    table = pd.read_csv(path, sep=" ", names=["Chr", "Start", "End", "Gene"])

    return table.values.tolist()


def get_cluster(path):
    table = pd.read_csv(path, sep="\t", names=["Chr", "Start", "End", "Name", "Score", "Strand", "Signal Value",
                                               "p value", "q value", "Peak"])
    table.drop(["Name", "Score", "Strand", "p value", "q value", "Peak"], axis=1, inplace=True)
    table.sort_values(by="Signal Value", ascending=False, inplace=True)
    table = table.head(10)

    return table.values.tolist()


def annotate_genes(top_peaks, gene_table):
    merged_lines = []
    for line in top_peaks:
        new_line = line
        for gene in gene_table:
            if line[0] == gene[0]:
                if overlaps([line[1], line[2]], [gene[1], gene[2]]) > 0:
                    if gene[3] not in new_line:
                        new_line.append(gene[3])
        merged_lines.append(new_line)

    return merged_lines


def main():
    pd.set_option('display.max_columns', None)
    cluster1_peaks = "/home/mfk/Google Drive/JLU/1. Semester/Datenanalyse/Gene_script/right-lobe-of-liver.1_peaks.narrowPeak"

    print(annotate_genes(get_cluster(cluster1_peaks), get_genecode()))


if __name__ == '__main__':
    main()
