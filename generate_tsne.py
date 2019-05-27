import pandas
from sklearn.manifold import TSNE

from sys import argv

pd_file = argv[1] 

csv_data = pandas.read_csv(pd_file, index_col=0)

tsne = TSNE(n_iter = 2000, perplexity = 50)
tsne_results = tsne.fit_transform(csv_data.values)
tsne_df = pandas.DataFrame(data = tsne_results, index = csv_data.index, columns = ["TSNE-1", "TSNE-2"])
tsne_df.to_csv(argv[2])
