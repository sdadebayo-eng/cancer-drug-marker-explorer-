"""Assemble index.html (a complete web page) from template.html, stats.js, app.js and gdsc_data.json."""
import pathlib
B = pathlib.Path(__file__).resolve().parent
frag = (B / 'template.html').read_text()
data = (B / 'gdsc_data.json').read_text()
stats = (B / 'stats.js').read_text()
app = (B / 'app.js').read_text()
assert '</script' not in data and '</script' not in stats and '</script' not in app
frag = frag.replace('/*DATA*/', data).replace('/*STATS*/', stats).replace('/*APP*/', app)
cut = frag.index('</style>') + len('</style>')
doc = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
       '<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">\n'
       + frag[:cut] + '\n</head>\n<body>\n' + frag[cut:] + '\n</body>\n</html>\n')
(B.parent / 'index.html').write_text(doc)
print('wrote', B.parent / 'index.html', round(len(doc.encode()) / 1e6, 2), 'MB')
