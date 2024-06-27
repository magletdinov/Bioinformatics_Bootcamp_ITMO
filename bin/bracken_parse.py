from pathlib import Path
import pandas as pd
import argparse

def create_path_to_root(taxid, nodes):
    try:
        path_to_root = []
        path_to_root.append(taxid)
        parent = nodes.loc[taxid]["parent tax_id"]
        if parent != 1:
            path_to_root.extend(create_path_to_root(parent, nodes))
        return path_to_root
    except:
        return ["None"]

parser = argparse.ArgumentParser(description='Process some files.')

parser.add_argument("-i", "--input", dest="INPUT", help="The path to the dir with the bracken reports", required=True)
parser.add_argument("-t", "--tax", dest="TAX", help="The path to the nodes.dmp", required=True)
args = parser.parse_args()

INPUT = Path(args.INPUT)
to_df = Path(args.TAX)
#to_df = Path("../ncbi_taxonomy/nodes.dmp")
#to_df = Path("nodes.dmp")
nodes = pd.read_csv(to_df, sep="\t", header=None, usecols=[0, 2, 4], names=["tax_id", "parent tax_id", "tax_name"], index_col="tax_id")

df_dict = {i.stem:pd.read_csv(i, sep="\t") for i in INPUT.glob("*.tsv")}
vir_df_list = []
for sample in df_dict:
    df = df_dict[sample]
    mask = df["taxonomy_id"].apply(lambda x: 10239 == create_path_to_root(x, nodes)[-1])
    df_vir = df[mask][["name", "taxonomy_id", "kraken_assigned_reads"]]
    df_vir["sample"] = ["_".join(sample.split("_")[:3]) for i in range(df_vir.shape[0])]
    vir_df_list.append(df_vir)
    
final_vir_df = pd.concat(vir_df_list).sort_values(by=["sample", "kraken_assigned_reads"], ascending=False)
final_vir_df.to_csv("bracken_vir_concat_mqc.tsv", sep="\t", index=False)