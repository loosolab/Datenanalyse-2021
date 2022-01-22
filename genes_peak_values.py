import os


def sum_peaks_per_gene(p_file, a_file):
    annot_peaks = {}
    with open(a_file, "r") as annot_file, open(p_file, "r") as peak_file:
        annot_file.readline()
        annot_lines = annot_file.readlines()
        peak_lines = peak_file.readlines()

    for i, p_line in enumerate(peak_lines):
        p_s = p_line.split()
        a_s = annot_lines[i].split()
        if a_s[21] not in annot_peaks.keys():
            annot_peaks[a_s[21]] = [p_s[6]]
        else:
            annot_peaks[a_s[21]].append(p_s[6])

    for gene in annot_peaks.keys():
        sum = 0
        for peak in annot_peaks[gene]:
            sum += float(peak)
        annot_peaks[gene] = sum

    annot_peaks = dict(sorted(annot_peaks.items(), key=lambda item: item[1], reverse=True))

    with open(p_file + "_sum", "w") as outfile:
        for gene in annot_peaks.keys():
            outfile.write(gene + "\t" + str(round(annot_peaks[gene])) + "\n")


def main():
    files = os.listdir()
    p_files = []
    a_files = []
    for file in files:
        if ".narrowPeak" in file:
            p_files.append(file)
        elif "finalhits.txt" in file:
            a_files.append(file)

    for p_file in p_files:
        p_name = p_file.split("_peaks")[0]
        for a_file in a_files:
            a_name = a_file.split("_peaks")[0]
            if p_name == a_name:
                print("match")
                sum_peaks_per_gene(p_file, a_file)


if __name__ == '__main__':
    main()
