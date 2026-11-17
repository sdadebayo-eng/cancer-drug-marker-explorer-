# Does the Gene Decide the Drug?

A teaching tool for BIOL 831. It tests whether a gene change in cancer cell lines predicts how well a drug works, using the public GDSC1000 dataset (release v17) from the Wellcome Sanger Institute.

`index.html` holds the whole tool: the data, the statistics and the page. The `media` folder holds the pictures and the video guide. The site needs no server code and no database.

## Put it online with GitHub Pages

1. Create a new public repository on GitHub, for example `drug-marker-explorer`.
2. Upload `index.html` and the `media` folder to the top level of the repository.
3. Open Settings, then Pages. Under "Build and deployment", pick "Deploy from a branch", choose `main` and `/ (root)`, and save.
4. After a minute or two the site is live at `https://<your-username>.github.io/drug-marker-explorer/`.

To update the site later, upload the new files with Add file, then Upload files. Files with the same name replace the old ones, and the site refreshes a minute or two after you commit.

Link to that address from the Program page of the class site. You can also open `index.html` straight from your computer by double-clicking it, as long as the `media` folder sits next to it. The fonts load from Google Fonts when you are online.

## What is in the `media` folder

| File | What it is |
|---|---|
| `guide.mp4` | The video guide, a narrated walkthrough of the Ask, Rank and How it works pages |
| `guide.vtt` | Captions for the video guide |
| `guide-poster.jpg` | The still shown before the video plays |
| `hero-cell.webp` | Thumbnail for the "Watch the guide" button at the top of every page |
| `dish-to-answer.webp` | The four lab steps, on the How it works page |
| `read-the-answer.webp` | How to read the dot chart, on the How it works page |
| `rank-funnel.webp` | How ranking finds the strongest links, on the How it works page |
| `not-patients.webp` | Why cell lines are not patients, in the Limits section |

## What is in the `source` folder

| File | What it does |
|---|---|
| `process.py` | Reads the GDSC files, removes poor curve fits, joins drugs, cell lines and markers, and writes `gdsc_data.json` |
| `compounds_v17.csv` | Drug names, targets and pathways for the 265 drug IDs, from `screened_compounds_rel_8.5.csv` on cancerrxgene.org |
| `stats.js` | The statistics: Mann-Whitney U test (exact for small groups), Hodges-Lehmann fold difference and 95% range, within-cancer-type model, Fisher exact test, Benjamini-Hochberg false discovery rate |
| `app.js` | The page logic: the pickers and their green and red marks, the answer, the charts, the checks, the Rank drugs page and the video guide |
| `template.html` | Page layout, styles and the How it works text |
| `build.py` | Puts the pieces together into `../index.html` |

## Example data file

`example-data/gdsc1000_v17_PLX4720_melanoma_example.csv` has 43 rows, one for each melanoma cell line tested with PLX-4720 (ID 1036), with BRAF and NRAS status, IC50, AUC and the highest dose. It is the clean example file for Block D of the worksheet.

## Rebuild the data yourself

```
cd source
git clone https://github.com/CancerRxGene/gdscdata
git clone https://github.com/CancerRxGene/gdsctools
pip install pandas numpy pyreadr
python process.py
python build.py
```

The rebuilt `index.html` should match the one in this folder exactly.

## Data sources

- Drug response: `gdsc_nlme_stats` in the gdscdata package (GPL-3), fitted by the GDSC team with the gdscIC50 package.
- Genetic markers and cell line details: `genomic_features_v17.csv.gz` and `cosmic_info.csv.gz` in gdsctools (BSD license).
- Drug names and targets: `screened_compounds_rel_8.5.csv`, cancerrxgene.org.

GDSC states that its data are freely available to the academic and medical communities. Cite Yang et al. 2013 (Nucleic Acids Res 41:D955) and Iorio et al. 2016 (Cell 166:740), and name the data version (GDSC1000 v17).

## Development log

The class grades individual work from the development log. Record this build session there: the date, the prompt, what was built, and what the team checked or changed afterwards. The illustrations were made by the team with an AI image tool. The video narration uses a synthetic voice (Kokoro text-to-speech).

For teaching only. Not for clinical use.
