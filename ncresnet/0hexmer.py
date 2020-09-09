from .FrameKmer import kmer_freq_file


def coding_nocoding_potential(input_file):
    coding = {}
    noncoding = {}
    for line in open(input_file).readlines():
        fields = line.split()
        if fields[0] == 'hexamer':
            continue
        coding[fields[0]] = float(fields[1])
        noncoding[fields[0]] = float(fields[2])
    return coding, noncoding


def maker_hexmer(nc_path, pc_path, hexmer_path):
    f = open(hexmer_path, "w")
    cod = kmer_freq_file(fastafile=pc_path, word_size=6, step_size=3, frame=0)
    noncod = kmer_freq_file(fastafile=nc_path, word_size=6, step_size=1, frame=0)
    cod_sum = 0.0
    cod_sum += sum(cod.values())
    noncod_sum = 0.0
    noncod_sum += sum(noncod.values())
    print('hexamer' + '\t' + 'coding' + '\t' + 'noncoding', file=f)
    for kmer in cod:
        if 'N' in kmer:
            continue
        print(kmer + '\t' + str(float(cod[kmer] / cod_sum)) + '\t' + str(float(noncod[kmer] / noncod_sum)), file=f)
    f.close()


# nc_path = "/home/ys/work/lncrna_identify_20190829/code/fasta/0Human/nc_human_train.fasta"
# pc_path = "/home/ys/work/lncrna_identify_20190829/code/fasta/0Human/pc_human_train.fasta"
#
# maker_hexmer(nc_path=nc_path, pc_path=pc_path, hexmer_path="hexmer_human.txt")
