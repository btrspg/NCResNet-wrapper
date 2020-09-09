from optparse import OptionParser

import keras
from Bio import SeqIO
from sklearn.preprocessing import StandardScaler

from ncresnet.calc_feature import get_feature
from ncresnet.utils import get_data


def get_features(seq_path):
    features = get_feature(seq_path=seq_path)
    return features


def get_seq_id(seq_path):
    records = list(SeqIO.parse(seq_path, format="fasta"))
    ids = [str(record.id) for record in records]
    return ids


'''
python test.py -i ./demo.fasta -m NCResNet.h5 -o result.tsv
'''

parser = OptionParser()
parser.add_option("-i", "--seq_path", dest="seq_path", help="sequences path with Fasta format")
parser.add_option("-m", "--model_file", dest="model_file", default=get_data() + 'NCResNet.h5',
                  help="NCResNet model file")
parser.add_option("-o", "--output_file", dest="output_file", help="Output file path with tsv format")
option, args = parser.parse_args()

model = keras.models.load_model(option.model_file, compile=False)
ids = get_seq_id(option.seq_path)
features = get_features(option.seq_path)
features = StandardScaler().fit_transform(X=features)
pred_res = model.predict(features).flatten()
pred_label = []
for x in pred_res:
    if x > 0.5:
        pred_label.append("pcRNA")
    else:
        pred_label.append("ncRNA")
f = open(option.output_file, "w")
f.write("id\tlabel\tscore\n")
for i in range(len(features)):
    print("%s\t%s\t%f" % (ids[i], pred_label[i], pred_res[i]), file=f)
f.close()
