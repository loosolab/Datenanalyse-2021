import math
import os
import statistics


def get_c_dict(clusters):
    """
    Count all genes which appear in any cluster
    :param clusters: dictionary
    :return: dictionary
    """
    c_dict = {}

    for c in clusters.keys():
        for gene in clusters[c].keys():
            if gene in c_dict.keys():
                c_dict[gene] += 1
            else:
                c_dict[gene] = 1

    return c_dict


def get_duplicates(c_dict, th):
    """
    Detect genes which appear in >= th clusters
    :param c_dict: dictionary
    :param th: int
    :return: list of strings
    """
    filtered = []

    for gene in c_dict.keys():
        if c_dict[gene] >= th:
            filtered.append(gene)

    return filtered


def get_panglao(tissue, connect=False):
    """
    Read and parse panglao database file
    :param tissue: string
    :param connect: boolean
    :return: dictionary
    """
    tissues = [tissue]
    path = "CTI/panglao_markers"
    panglao_dict = {}
    panglao_rank_dict = {}

    if connect:
        tissues.append("tissue")

    with open(path, "r") as panglao_file:
        panglao_file.readline()
        for line in panglao_file.readlines():
            species, gene_symb, ct, n_genes, ub_i, organ = line.split("\t")
            us = float(ub_i)
            if us != 0:
                us = round(math.sqrt(1 / us))
            else:
                us = 32
            if any(t in organ.lower() for t in tissues) and "Hs" in species:
                if ct not in panglao_dict.keys():
                    panglao_dict[ct] = []
                genes = [(us, gene_symb)]
                if len(n_genes.split("|")) > 1:
                    for gene in n_genes.split("|"):
                        genes.append((us, gene.upper()))
                elif n_genes != "NA":
                    genes.append((us, n_genes))
                for gene in genes:
                    panglao_dict[ct].append(gene)

    for ct in panglao_dict.keys():
        rank_dict = {}
        for gene in panglao_dict[ct]:
            if gene[1] in rank_dict.keys():
                rank_dict[gene[1]] = gene[0]
            else:
                rank_dict[gene[1]] = gene[0]

        panglao_rank_dict[ct] = rank_dict

    return panglao_rank_dict


def get_cell_types(cpath, tissue, connect=False):
    """
    Prepare database and clusters for upcoming ranking calculations
    :param cpath: string
    :param tissue: string
    :param connect: boolean
    :return: dictionary
    """
    clusters = get_annotated_clusters(path=cpath)
    dupls = get_duplicates(get_c_dict(clusters), math.ceil(len(clusters) / 1.5))
    pang_dict = get_panglao(tissue, connect=connect)

    filtered_clusters = {}
    for cluster_num in clusters.keys():
        filtered_clusters[cluster_num] = {}
        for gene in clusters[cluster_num].keys():
            if gene not in dupls:
                filtered_clusters[cluster_num][gene] = clusters[cluster_num][gene]

    mean_peaks = {}
    for cluster_num in filtered_clusters.keys():
        mean_peaks[cluster_num] = (statistics.mean(clusters[cluster_num].values()))

    def calc_ranks(cm_dict):
        """
        Identify cell types of each cluster by ranking each fitting gene by peak value,
        quantity and panglao ubiquitousness index
        :param cm_dict: dictionary
        :return: dictionary
        """
        ct_dict = {}

        for key in clusters.keys():
            ct_dict[key] = {}

        for celltype in cm_dict.keys():
            gene_count = len(cm_dict[celltype])
            for c in filtered_clusters.keys():
                count = 0
                ranks = []
                ub_scores = []

                for mgene in filtered_clusters[c].keys():
                    if filtered_clusters[c][mgene] > mean_peaks[c]:
                        if mgene in cm_dict[celltype].keys():
                            ranks.append(cm_dict[celltype][mgene])
                            ub_scores.append(cm_dict[celltype][mgene])
                            count += 1

                if count > 4:
                    ranks = ranks[:5]

                if count > 0:
                    ub_score = round(statistics.mean(ub_scores))
                    ct_dict[c][celltype] = [round(sum(ranks) * ub_score * count / gene_count), count, gene_count,
                                            ub_score]

        return ct_dict

    return calc_ranks(pang_dict)


def annot_peaks(apath="CTI/annot/", ppath="CTI/npeaks/", cpath="CTI/cluster/"):
    """
    Annotate parsed narrow peak files with parsed uropa finalhits files and save annotated cluster files to cpath
    :param apath: string
    :param ppath: string
    :param cpath: string
    """
    for afile in os.listdir(apath):
        with open(cpath + afile.split("_peaks")[0], "w") as cfile:
            for pfile in os.listdir(ppath):
                if afile.split("_peaks")[0] == pfile.split("_peaks")[0]:
                    with open(apath + afile) as a_file:
                        alines = a_file.readlines()
                    with open(ppath + pfile) as p_file:
                        plines = p_file.readlines()
                    for i in range(len(alines)):
                        cfile.write(alines[i].split("\t")[3].rstrip() + "\t" + plines[i].split("\t")[3])


def get_annotated_clusters(path="CTI/liver/cluster/"):
    """
    Read cluster files and sum peak values if genes appear more than once per file
    :param path: string
    :return: dictionary
    """
    annotated_clusters = {}
    for file in os.listdir(path):
        cname = file.split(".")[1]
        annotated_dict = {}
        with open(path + file) as cfile:
            annotated_dict[cname] = []
            lines = cfile.readlines()
            for line in lines:
                split = line.split("\t")
                annotated_dict[cname].append([split[0], float(split[1].rstrip())])

        sum_dict = {}
        for gene in annotated_dict[cname]:
            if gene[0] in sum_dict.keys():
                sum_dict[gene[0]] += gene[1]
            else:
                sum_dict[gene[0]] = gene[1]

        annotated_clusters[cname] = sum_dict

    return annotated_clusters


def perform_cell_type_annotation(sample, tissue, connect=False):
    """
    Perform cell type identification and annotation, create cluster folder, generate cell type assignment table
    and create ranks folder with files for further investigation
    :param sample: string
    :param tissue: string
    :param connect: boolean
    """
    apath, ppath, cpath, opath = ["CTI/" + sample + "/annot/", "CTI/" + sample + "/npeaks/",
                                  "CTI/" + sample + "/cluster/", "CTI/" + sample + "/ranks/"]

    if not os.path.exists(cpath):
        os.makedirs(cpath)
        annot_peaks(apath=apath, ppath=ppath, cpath=cpath)

    if not os.path.exists(opath):
        os.makedirs(opath)

    ct_dict = get_cell_types(cpath, tissue, connect=connect)
    with open("CTI/" + sample + "/annotation.txt", "w") as c_file:
        for dic in ct_dict.keys():
            sorted_dict = dict(sorted(ct_dict[dic].items(), key=lambda r: (r[1][0], r[1][1]), reverse=True))
            with open(opath + "cluster_" + dic, "w") as d_file:
                for key in sorted_dict.keys():
                    d_file.write(key)
                    for value in sorted_dict[key]:
                        d_file.write("\t" + str(value))
                    d_file.write("\n")
            if len(sorted_dict.keys()) > 0:
                c_file.write(dic + "\t" + str(next(iter(sorted_dict))) + "\n")


def main():
    perform_cell_type_annotation("colon", "gi tract", connect=True)
    print("done")
    perform_cell_type_annotation("stomach", "gi tract", connect=True)
    print("done")
    perform_cell_type_annotation("liver", "liver", connect=True)
    print("done")
    perform_cell_type_annotation("lung", "lungs", connect=True)
    print("done")
    perform_cell_type_annotation("fat", "adipose tissue", connect=True)
    print("done")
    perform_cell_type_annotation("muscle", "muscle", connect=True)
    print("done")
    perform_cell_type_annotation("muscle2", "muscle", connect=True)
    print("done")
    perform_cell_type_annotation("skin", "skin", connect=True)
    print("done")
    perform_cell_type_annotation("pancreas", "pancreas")
    print("done")


if __name__ == '__main__':
    main()
