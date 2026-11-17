"""Build the browser data file for the drug-marker explorer from official GDSC1000 (v17) files.

Inputs (from the GDSC team's GitHub repositories):
  gdscdata/data/gdsc_nlme_stats.rda            fitted dose-response curves (gdscIC50 package)
  gdsctools/data/genomic_features_v17.csv.gz   binary genetic features (mutations, copy number, fusions, MSI)
  gdsctools/data/cosmic_info.csv.gz            cell line names and cancer type labels
  compounds_v17.csv                            drug names/targets from screened_compounds_rel_8.5.csv (cancerrxgene.org)
"""
import base64, json, pathlib, re
import numpy as np, pandas as pd, pyreadr

SRC = pathlib.Path(__file__).resolve().parent   # clone gdscdata and gdsctools into this folder
fits = pyreadr.read_r(f'{SRC}/gdscdata/data/gdsc_nlme_stats.rda')['gdsc_nlme_stats']
fits = fits.drop_duplicates(subset=['CL', 'drug'])[['CL', 'CELL_LINE_NAME', 'DRUG_ID_lib', 'maxc', 'IC50', 'auc', 'RMSE']]
fits['DRUG_ID'] = fits.DRUG_ID_lib.astype(int)
fits['COSMIC_ID'] = fits.CL.astype(int)
n_all = len(fits)
fits = fits[fits.RMSE <= 0.3]                      # GDSC release rule: drop poorly fitted curves
n_rmse = n_all - len(fits)
# one concentration range per drug: keep the range used for most cell lines
maj = (fits.groupby(['DRUG_ID', 'maxc']).size().reset_index(name='n')
       .sort_values(['DRUG_ID', 'n'], ascending=[True, False]).drop_duplicates('DRUG_ID'))
n_before = len(fits)
fits = fits.merge(maj[['DRUG_ID', 'maxc']], on=['DRUG_ID', 'maxc'])
n_minor = n_before - len(fits)

gf = pd.read_csv(f'{SRC}/gdsctools/gdsctools/data/genomic_features_v17.csv.gz')
ci = pd.read_csv(f'{SRC}/gdsctools/gdsctools/data/cosmic_info.csv.gz')
comp = pd.read_csv(f'{SRC}/compounds_v17.csv', keep_default_na=False)

cells = sorted(set(fits.COSMIC_ID) & set(gf.COSMIC_ID) & set(ci.COSMIC_ID))
ci = ci.set_index('COSMIC_ID').loc[cells].reset_index()
gf = gf.set_index('COSMIC_ID').loc[cells].reset_index()

TCGA = {
 'LUAD': 'Lung adenocarcinoma', 'LUSC': 'Lung squamous cell carcinoma', 'SCLC': 'Small cell lung cancer',
 'SKCM': 'Melanoma', 'BRCA': 'Breast cancer', 'COREAD': 'Colorectal cancer', 'HNSC': 'Head and neck cancer',
 'ESCA': 'Esophageal cancer', 'GBM': 'Glioblastoma', 'LGG': 'Lower-grade glioma',
 'DLBC': 'Diffuse large B-cell lymphoma', 'OV': 'Ovarian cancer', 'KIRC': 'Kidney cancer', 'NB': 'Neuroblastoma',
 'PAAD': 'Pancreatic cancer', 'LAML': 'Acute myeloid leukemia', 'ALL': 'Acute lymphoblastic leukemia',
 'STAD': 'Stomach cancer', 'MESO': 'Mesothelioma', 'BLCA': 'Bladder cancer', 'MM': 'Multiple myeloma',
 'LIHC': 'Liver cancer', 'THCA': 'Thyroid cancer', 'CESC': 'Cervical cancer', 'LCML': 'Chronic myeloid leukemia',
 'UCEC': 'Endometrial cancer', 'PRAD': 'Prostate cancer', 'MB': 'Medulloblastoma',
 'CLL': 'Chronic lymphocytic leukemia', 'ACC': 'Adrenocortical cancer'}
GDSC2 = {
 'ewings_sarcoma': 'Ewing sarcoma', 'osteosarcoma': 'Osteosarcoma', 'rhabdomyosarcoma': 'Rhabdomyosarcoma',
 'chondrosarcoma': 'Other sarcoma', 'fibrosarcoma': 'Other sarcoma', 'soft_tissue_other': 'Other sarcoma', 'bone_other': 'Other sarcoma',
 'Burkitt_lymphoma': 'Burkitt lymphoma', 'Hodgkin_lymphoma': 'Hodgkin lymphoma',
 'anaplastic_large_cell_lymphoma': 'Other lymphoma', 'lymphoid_neoplasm other': 'Other lymphoma',
 'B_cell_leukemia': 'Other leukemia', 'hairy_cell_leukaemia': 'Other leukemia', 'leukemia': 'Other leukemia',
 'lymphoblastic_leukemia': 'Other leukemia', 'T_cell_leukemia': 'Other leukemia', 'haematopoietic_neoplasm other': 'Other leukemia',
 'lung_NSCLC_large cell': 'Other lung cancer', 'lung_NSCLC_not specified': 'Other lung cancer', 'lung_NSCLC_carcinoid': 'Other lung cancer',
 'Lung_other': 'Other lung cancer', 'lung_NSCLC_adenocarcinoma': 'Other lung cancer',
 'ovary': 'Ovarian cancer', 'stomach': 'Stomach cancer', 'pancreas': 'Pancreatic cancer', 'head and neck': 'Head and neck cancer',
 'endometrium': 'Endometrial cancer', 'breast': 'Breast cancer', 'kidney': 'Kidney cancer', 'cervix': 'Cervical cancer',
 'prostate': 'Prostate cancer', 'biliary_tract': 'Biliary tract cancer',
 'testis': 'Other solid tumor', 'urogenital_system_other': 'Other solid tumor', 'skin_other': 'Other solid tumor'}
BLOOD = {'Diffuse large B-cell lymphoma', 'Acute myeloid leukemia', 'Acute lymphoblastic leukemia', 'Multiple myeloma',
         'Chronic myeloid leukemia', 'Chronic lymphocytic leukemia', 'Burkitt lymphoma', 'Hodgkin lymphoma',
         'Other lymphoma', 'Other leukemia'}

def ctype(row):
    t = row.TCGA if isinstance(row.TCGA, str) else ''
    if t and t != 'UNABLE TO CLASSIFY':
        return TCGA[t]
    return GDSC2[row.GDSC2]
ci['ctype'] = ci.apply(ctype, axis=1)
type_names = sorted(ci.ctype.unique())
type_idx = {t: i for i, t in enumerate(type_names)}

# drugs
comp = comp.set_index('DRUG_ID')
drug_ids = sorted(fits.DRUG_ID.unique())
dup_names = comp.loc[drug_ids].DRUG_NAME.value_counts()
dup_names = set(dup_names[dup_names > 1].index)
maxc = fits.groupby('DRUG_ID').maxc.first()

cidx = {c: i for i, c in enumerate(cells)}
NC, ND = len(cells), len(drug_ids)
ic = np.full((ND, NC), np.nan, dtype=np.float32)
au = np.full((ND, NC), 65535, dtype=np.uint16)
sub = fits[fits.COSMIC_ID.isin(cidx)]
for di, d in enumerate(drug_ids):
    s = sub[sub.DRUG_ID == d]
    ii = s.COSMIC_ID.map(cidx).values
    ic[di, ii] = s.IC50.values.astype(np.float32)
    au[di, ii] = np.round(s.auc.values * 10000).astype(np.uint16)

drugs = []
for d in drug_ids:
    r = comp.loc[d]
    twin = [int(x) for x in drug_ids if x != d and comp.loc[x].DRUG_NAME == r.DRUG_NAME]
    syn = '' if r.SYNONYMS in ('NA',) else r.SYNONYMS
    drugs.append({'id': int(d), 'name': r.DRUG_NAME, 'syn': syn, 'target': r.TARGET, 'path': r.TARGET_PATHWAY,
                  'maxc': float(maxc.loc[d]), 'twin': twin})

# features
DRIVERS = ["ERBB2","EGFR","MET","MYC","MYCN","CCND1","CCNE1","CDK4","CDK6","MDM2","KRAS","FGFR1","FGFR2","FGFR3",
           "PDGFRA","KIT","AKT1","AKT2","PIK3CA","PIK3CB","SOX2","TERT","MITF","AR","ESR1","NKX2-1","IL7R","MECOM",
           "GATA3","CDKN2A","PTEN","RB1","SMAD4","STK11","TP53","NF1","NF2","APC","WT1","BAP1","SMARCA4","ARID1A","FHIT",
           "WWOX","PTPRD","LRP1B","MAP2K4","CDH1","KDM6A","NOTCH1","FAT1","CIC","GNA11","GNAQ","RET","FAS","JAK2"]
feats = []
def add(f, col):
    idx = [i for i, v in enumerate(col) if v == 1]
    if len(idx) == 0:
        return
    f['cells'] = idx
    feats.append(f)
for c in gf.columns:
    if c.endswith('_mut'):
        g = c[:-4]
        if g in ('BCR-ABL', 'EWSR1-FLI1', 'EWSR1-X'):
            nice = {'BCR-ABL': 'BCR-ABL', 'EWSR1-FLI1': 'EWSR1-FLI1', 'EWSR1-X': 'EWSR1 (other partner)'}[g]
            add({'key': c, 'kind': 'fusion', 'gene': nice, 'genes': g.split('-'), 'label': f'{nice} fusion', 'note': '',
                 'altShort': f'{nice} fusion', 'refShort': 'No fusion',
                 'altName': f'lines with the {nice} fusion', 'refName': f'lines without the {nice} fusion',
                 'statusName': f'{nice} fusion status', 'has': f'has the {nice} fusion'}, gf[c].values)
        else:
            add({'key': c, 'kind': 'mut', 'gene': g, 'genes': [g], 'label': f'{g} mutation', 'note': '',
                 'altShort': f'{g} mutant', 'refShort': f'No {g} mutation',
                 'altName': f'{g}-mutant lines', 'refName': f'lines with no {g} mutation',
                 'statusName': f'{g} mutation status', 'has': f'is {g}-mutant'}, gf[c].values)
    elif c.startswith(('gain_', 'loss_')):
        m = re.match(r'(gain|loss)_cnaPANCAN(\d+)_\((.+)\)$', c)
        if not m:
            continue
        kind, reg, genes = m.group(1), m.group(2), m.group(3).split(',')
        head = next((g for g in DRIVERS if g in genes), None)
        if head is None:
            head = ', '.join(genes) if len(genes) <= 4 else f"{', '.join(genes[:3])} and {len(genes)-3} more"
        others = [g for g in genes if g != head]
        word = 'amplification' if kind == 'gain' else 'deletion'
        adj = 'amplified' if kind == 'gain' else 'deleted'
        note = f"Copy-number region {reg}" + (f", also contains {', '.join(others)}" if others and head in genes else '')
        add({'key': c, 'kind': 'amp' if kind == 'gain' else 'del', 'gene': head, 'genes': genes,
             'label': f'{head} {word}', 'note': note, 'region': int(reg),
             'altShort': f'{head} {adj}', 'refShort': f'Not {adj}',
             'altName': f'{head}-{adj} lines' if ',' not in head else f'lines with the {head} {word}',
             'refName': f'lines without the {head} {word}',
             'statusName': f'{head} {word} status', 'has': f'has the {head} {word}'}, gf[c].values)
add({'key': 'MSI', 'kind': 'msi', 'gene': 'MSI', 'genes': ['MSI'], 'label': 'Microsatellite instability (MSI-high)',
     'note': 'Mismatch repair deficiency', 'altShort': 'MSI-high', 'refShort': 'Not MSI-high',
     'altName': 'MSI-high lines', 'refName': 'lines without MSI-high status',
     'statusName': 'MSI status', 'has': 'is MSI-high'}, gf['MSI_FACTOR'].values)
# make duplicate labels unique
seen = {}
for f in feats:
    seen.setdefault(f['label'], []).append(f)
for lab, fl in seen.items():
    if len(fl) > 1:
        for f in fl:
            if 'region' in f:
                f['label'] = f"{lab} (region {f['region']})"
def b64(a):
    return base64.b64encode(a.tobytes()).decode()

data = {
  'meta': {
    'dataset': 'GDSC1000 (release v17)', 'screen': 'GDSC1',
    'fits_total': int(n_all), 'fits_removed_rmse': int(n_rmse), 'fits_removed_minor_range': int(n_minor),
    'pairs': int((~np.isnan(ic)).sum()), 'cells': NC, 'drugs': ND, 'features': len(feats),
    'extrapolated_fraction': float(np.nanmean(np.where(np.isnan(ic), np.nan, (ic > np.log(np.array([d['maxc'] for d in drugs]))[:, None]).astype(float))))
  },
  'types': type_names,
  'typeCounts': [int((ci.ctype == t).sum()) for t in type_names],
  'blood': [type_idx[t] for t in type_names if t in BLOOD],
  'cells': {'name': [ci.set_index('COSMIC_ID').loc[c].SAMPLE_NAME for c in cells], 'cosmic': [int(c) for c in cells],
            'type': [type_idx[t] for t in ci.ctype]},
  'drugs': drugs, 'features': feats,
  'ic50': b64(ic), 'auc': b64(au)
}
with open(SRC / 'gdsc_data.json', 'w') as fh:
    json.dump(data, fh, separators=(',', ':'))
print({k: v for k, v in data['meta'].items()})
print('types:', len(type_names))
print('file size MB:', round(len(json.dumps(data, separators=(",", ":")))/1e6, 2))
print('features by kind:', pd.Series([f['kind'] for f in feats]).value_counts().to_dict())
