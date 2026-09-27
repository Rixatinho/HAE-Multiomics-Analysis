import pandas as pd, numpy as np, os
from scipy import stats

B = "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学/02_analysis/results"
E = f"{B}/enhancement14_meta_analysis"
OUT = f"{B}/enhancement52_human_subgroup_meta"
os.makedirs(OUT, exist_ok=True)

df = pd.read_csv(f"{E}/meta_analysis_all_genes.csv")
print("total genes:", len(df))

# human cohorts (from effective_sample_summary): HAE_our, GSE124362
lg1,se1 = df["logFC_HAE_our"], df["SE_HAE_our"]
lg2,se2 = df["logFC_GSE124362"], df["SE_GSE124362"]
mask = lg1.notna()&lg2.notna()&se1.notna()&se2.notna()&(se1>0)&(se2>0)
sub = df[mask].copy()
print("genes with both human cohorts:", mask.sum())

w1=1/se1[mask]**2; w2=1/se2[mask]**2
mu=(w1*lg1[mask]+w2*lg2[mask])/(w1+w2)
se=np.sqrt(1/(w1+w2))
z=mu/se; p=2*stats.norm.sf(np.abs(z))
Q=w1*(lg1[mask]-mu)**2+w2*(lg2[mask]-mu)**2
I2=np.maximum(0,(Q-1)/np.maximum(Q,1e-12))*100
agree=(np.sign(lg1[mask])==np.sign(lg2[mask])).astype(int)

sub["human_logFC"]=mu; sub["human_SE"]=se; sub["human_P"]=p
sub["human_I2"]=I2; sub["human_dir_agree"]=agree
sub["human_padj"]=stats.false_discovery_control(p) if hasattr(stats,'false_discovery_control') else np.minimum(1,p*len(p)/ (np.argsort(np.argsort(p))+1))

# key genes
key=["CYP3A4","CYP3A5","CYP3A7","CYP1A2","CYP2C9","CYP2C19","CYP2D6","CYP2E1",
     "HNF1A","HNF4A","FMO3","FMO1","ALB","CPS1","ASS1","HAL","GLUL","CYP2A6","UGT2B7","ABCB11"]
kg=sub[sub["gene"].str.strip('"').isin(key)].copy()
kg["gene"]=kg["gene"].str.strip('"')
cols=["gene","logFC_HAE_our","P_HAE_our","logFC_GSE124362","P_GSE124362","human_logFC","human_SE","human_P","human_I2","human_dir_agree"]
kg=kg[cols].sort_values("human_P")

# summary stats
n_human = mask.sum()
n_sig = (sub["human_P"]<0.05).sum()
n_sig_agree = ((sub["human_P"]<0.05)&(sub["human_dir_agree"]==1)).sum()
med_I2 = sub["human_I2"].median()
print(f"\nhuman-subgroup meta: n_genes={n_human}, P<0.05={n_sig}, P<0.05 & dir-agree={n_sig_agree}, median I2={med_I2:.1f}%")
print("\n=== KEY GENES (human-subgroup fixed-effect meta) ===")
print(kg.to_string(index=False))

sub.to_csv(f"{OUT}/human_subgroup_meta_all.csv", index=False)
kg.to_csv(f"{OUT}/human_subgroup_meta_key_genes.csv", index=False)
with open(f"{OUT}/summary.txt","w") as f:
    f.write(f"human cohorts: HAE_our (12v12) + GSE124362 (6v6); fixed-effect, 2-cohort\n")
    f.write(f"genes with both: {n_human}; human_P<0.05: {n_sig}; sig&dir-agree: {n_sig_agree}; median I2: {med_I2:.1f}%\n")
print("\nWROTE", OUT)
