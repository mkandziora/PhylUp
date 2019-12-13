import pandas as pd

colnames = ['accession', 'ncbi_txid', 'ncbi_txn']
fn = 'Senecioneae.acc'
data = pd.read_csv(fn, names=colnames, sep=" ", header=None)

print(len(data['accession']))
print(len(data['ncbi_txid']))
#print(data['ncbi_txid'].tolist())
print(len(set(data['ncbi_txid'].tolist())))


