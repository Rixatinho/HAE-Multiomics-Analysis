import pandas as pd, numpy as np
from scipy import stats
import statsmodels.api as sm

B = "/Users/rishat/Library/Mobile Documents/com~apple~CloudDocs/个人文档/Word/文献写作/2.肝包虫/20260317-肝包虫多组学/02_analysis/results"
zs = pd.read_csv(f"{B}/zonation_collapse/zonation_scores_per_sample.csv")
bp = pd.read_csv(f"{B}/enhancement16_deconvolution/bayesprism_cell_proportions.csv", index_col=0)
bp.index = [str(x).strip('"') for x in bp.index]

rows=[]
for _,r in zs.iterrows():
    p=int(r['patient'])
    rows.append(dict(patient=p,
        hep_n=bp.loc[f'Normal{p}','Hepatocyte'], hep_a=bp.loc[f'Adjacent{p}','Hepatocyte'],
        pp_n=r['periportal_normal'], pp_a=r['periportal_adjacent'],
        pc_n=r['pericentral_normal'], pc_a=r['pericentral_adjacent'],
        zdi_n=r['ZDI_normal'], zdi_a=r['ZDI_adjacent']))
df=pd.DataFrame(rows)

def dz(a,b):
    d=np.asarray(b)-np.asarray(a); return d.mean()/d.std(ddof=1)
def wp(a,b):
    try: return stats.wilcoxon(b,a).pvalue
    except Exception: return np.nan

print("n pairs =", len(df), "patients:", list(df.patient))

print("\n=== 1. Hepatocyte composition (paired, Adjacent vs Normal) ===")
print(f"Normal {df.hep_n.mean():.4f}+-{df.hep_n.std():.4f}  Adjacent {df.hep_a.mean():.4f}+-{df.hep_a.std():.4f}")
print(f"delta={(df.hep_a-df.hep_n).mean():.4f}  Wilcoxon p={wp(df.hep_n,df.hep_a):.4f}  dz={dz(df.hep_n,df.hep_a):.3f}")
df2=df[df.patient!=3]
print(f"[excl P3 outlier] Adjacent {df2.hep_a.mean():.4f}  p={wp(df2.hep_n,df2.hep_a):.4f}  dz={dz(df2.hep_n,df2.hep_a):.3f}")

print("\n=== 2. Zone scores (paired, replication) ===")
for name,n,a in [('Periportal','pp_n','pp_a'),('Pericentral','pc_n','pc_a'),('ZDI gradient','zdi_n','zdi_a')]:
    print(f"{name}: Normal {df[n].mean():.3f} -> Adjacent {df[a].mean():.3f}  p={wp(df[n],df[a]):.4f}  dz={dz(df[n],df[a]):.3f}")

print("\n=== 3. Composition-adjusted (ANCOVA on change scores, intercept=adjusted change) ===")
for name,n,a in [('Periportal','pp_n','pp_a'),('Pericentral','pc_n','pc_a'),('ZDI','zdi_n','zdi_a')]:
    dzone=(df[a]-df[n]).values; dhep=(df['hep_a']-df['hep_n']).values
    X=sm.add_constant(dhep); m=sm.OLS(dzone,X).fit()
    rho,rp=stats.spearmanr(dzone,dhep)
    print(f"{name}: adj_change={m.params[0]:.3f}(p={m.pvalues[0]:.4f})  hep_slope={m.params[1]:.3f}(p={m.pvalues[1]:.4f})  cor rho={rho:.3f}(p={rp:.4f})")

print("\n=== 4. Gradient preservation (periportal - pericentral) ===")
g_n=df.pp_n-df.pc_n; g_a=df.pp_a-df.pc_a
print(f"Gradient: Normal {g_n.mean():.3f} -> Adjacent {g_a.mean():.3f}  p={wp(g_n,g_a):.4f}  dz={dz(g_n,g_a):.3f}")

print("\n=== 5. Per-hepatocyte normalised zone score (sensitivity) ===")
for name,n,a in [('Periportal','pp_n','pp_a'),('Pericentral','pc_n','pc_a')]:
    nh=df[n]/df['hep_n']; ah=df[a]/df['hep_a']
    print(f"{name}/hep: Normal {nh.mean():.2f} -> Adjacent {ah.mean():.2f}  p={wp(nh,ah):.4f}  dz={dz(nh,ah):.3f}")
