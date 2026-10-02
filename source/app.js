(function () {
  'use strict';
  var GS = window.GS;

  // ------------------------------------------------------------------ data
  var RAW = JSON.parse(document.getElementById('gdsc-data').textContent);
  function b64bytes(s) {
    var bin = atob(s), n = bin.length, u = new Uint8Array(n);
    for (var i = 0; i < n; i++) u[i] = bin.charCodeAt(i);
    return u;
  }
  var IC = new Float32Array(b64bytes(RAW.ic50).buffer);
  var AUC16 = new Uint16Array(b64bytes(RAW.auc).buffer);
  var NC = RAW.cells.name.length;
  var CELLNAME = RAW.cells.name, CELLTYPE = RAW.cells.type, COSMIC = RAW.cells.cosmic;
  var TYPES = RAW.types, NT = TYPES.length;
  var BLOOD = new Uint8Array(NT);
  RAW.blood.forEach(function (t) { BLOOD[t] = 1; });
  var TYPECOUNT = RAW.typeCounts;
  var LN10 = Math.LN10;

  var DRUGS = RAW.drugs, DRUG_BY_ID = {};
  DRUGS.forEach(function (d, i) {
    d.i = i; d.lnMax = Math.log(d.maxc); DRUG_BY_ID[d.id] = d;
    var n = 0, off = i * NC;
    for (var c = 0; c < NC; c++) if (IC[off + c] === IC[off + c]) n++;
    d.nTested = n;
    d.hay = (d.name + ' ' + d.syn + ' ' + d.target + ' ' + d.path + ' ' + d.id).toLowerCase();
  });
  var ALIAS = { ERBB2: 'HER2' };
  var FEATS = RAW.features, FEAT_BY_KEY = {};
  FEATS.forEach(function (f, i) {
    f.i = i; FEAT_BY_KEY[f.key] = f;
    f.mask = new Uint8Array(NC);
    f.cells.forEach(function (c) { f.mask[c] = 1; });
    var al = f.genes.map(function (g) { return ALIAS[g]; }).filter(Boolean);
    if (ALIAS[f.gene]) f.label = f.label.replace(f.gene, f.gene + ' (' + ALIAS[f.gene] + ')');
    f.hay = (f.label + ' ' + f.genes.join(' ') + ' ' + al.join(' ') + ' ' + f.kind + ' ' + (f.note || '')).toLowerCase();
    f.altOne = f.altName.replace(/lines\b/, 'line');
    f.refOne = f.refName.replace(/lines\b/, 'line');
  });

  // ------------------------------------------------------------------ helpers
  function $(id) { return document.getElementById(id); }
  function esc(s) {
    return String(s).replace(/[&<>"']/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c];
    });
  }
  // keeps hyphenated names such as PLX-4720 or BRAF-mutant on one line
  function nw(s) { return esc(s).replace(/(\S*-\S*)/g, '<span class="nw">$1</span>'); }
  var NF = new Intl.NumberFormat('en-US');
  function nf(n) { return NF.format(n); }
  function cap(s) { return s.charAt(0).toUpperCase() + s.slice(1); }
  function plural(n, one, many) { return n === 1 ? one : (many || one + 's'); }
  function fmtConc(v) {
    if (!isFinite(v)) return 'n/a';
    if (v >= 1000) return nf(Math.round(v));
    if (v >= 100) return String(Math.round(v));
    if (v >= 10) return v.toFixed(1);
    if (v >= 1) return v.toFixed(2);
    return String(Number(v.toPrecision(2)));
  }
  function fmtFold(f) {
    if (!isFinite(f)) return 'n/a';
    if (f >= 1000) return nf(Math.round(f));
    if (f >= 10) return String(Math.round(f));
    return f.toFixed(1);
  }
  function fmtAUC(a) { return isFinite(a) ? a.toFixed(2) : 'n/a'; }
  var SUP = { '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹', '-': '⁻' };
  function sup(s) { return String(s).split('').map(function (c) { return SUP[c] || c; }).join(''); }
  function fmtP(p) {
    if (!isFinite(p)) return 'n/a';
    if (p >= 0.0001) return String(Number(p.toPrecision(2)));
    if (p < 1e-300) return '< 10' + sup(-300);
    var e = Math.floor(Math.log10(p)), m = p / Math.pow(10, e);
    if (m >= 9.95) { m = 1; e += 1; }
    return m.toFixed(1) + ' × 10' + sup(e);
  }
  function pEq(p) { return (p < 1e-300 ? 'p ' : 'p = ') + fmtP(p); }
  function drugText(d) { var n = d.name.replace(/ \(-\)$/, ''); return /^[A-Z][a-z]+$/.test(n) ? n.toLowerCase() : n; }
  function typeNoun(name) {
    if (/^(Ewing|Burkitt|Hodgkin)/.test(name)) return name;
    return name.charAt(0).toLowerCase() + name.slice(1);
  }
  // "are NRAS-mutant" / "are not NRAS-mutant" / "have the ERBB2 amplification" / "do not have ..."
  function hasPlural(f, keep) {
    var h = f.has;
    if (h.indexOf('is ') === 0) return (keep ? 'are ' : 'are not ') + h.slice(3);
    if (h.indexOf('has ') === 0) return (keep ? 'have ' : 'do not have ') + h.slice(4);
    return h;
  }
  function hasOne(f, keep) {
    var h = f.has;
    if (keep) return h;
    if (h.indexOf('is ') === 0) return 'is not ' + h.slice(3);
    if (h.indexOf('has ') === 0) return 'does not have ' + h.slice(4);
    return h;
  }
  function nLines(n, f, alt) { return n + ' ' + (alt ? (n === 1 ? f.altOne : f.altName) : (n === 1 ? f.refOne : f.refName)); }
  function missText(R) {
    return nLines(R.missAlt, R.feat, 1) + ' ' + R.S.where + (R.F ? ' that ' + (R.missAlt === 1 ? hasOne(R.F.f, R.F.keep) : hasPlural(R.F.f, R.F.keep)) : '');
  }
  // "BRAF-mutant melanoma lines", "breast cancer lines with the BCR-ABL fusion"
  function groupNounT(f, t) {
    var tn = typeNoun(TYPES[t]);
    if (/^lines /.test(f.altName)) return tn + ' ' + f.altName;
    return f.altName.replace(/ lines$/, ' ' + tn + ' lines');
  }
  function groupNoun(f, S) { return S.pooled ? f.altName : groupNounT(f, S.type); }
  function otherNoun(S) { return S.pooled ? 'other lines of the same cancer type' : 'other ' + S.noun; }
  function foldX(m) { return m.fold >= 1 ? m.fold : 1 / m.fold; }
  function lessMore(m) { return m.fold >= 1 ? 'less' : 'more'; }
  function chanceWords(p) {
    if (p < 0.001) return 'Chance alone would produce a gap this big less than 1 time in 1,000.';
    if (p < 0.01) return 'Chance alone would produce a gap this big about ' + Math.max(1, Math.round(p * 1000)) + ' times in 1,000.';
    if (p < 0.05) return 'Chance alone would produce a gap this big about ' + Math.max(1, Math.round(p * 100)) + ' times in 100.';
    return 'Chance alone would produce a gap this big about ' + Math.round(p * 100) + ' times in 100.';
  }
  var LABEL = { sens: 'Gene change predicts sensitivity', res: 'Gene change predicts resistance', small: 'Small difference', none: 'No clear link' };

  // ------------------------------------------------------------------ sets of cancer types
  function setInfo(val) {
    var mask = new Uint8Array(NT), i, n = 0;
    if (val === 'all' || val === 'solid' || val === 'blood') {
      for (i = 0; i < NT; i++) {
        mask[i] = val === 'all' ? 1 : val === 'blood' ? BLOOD[i] : 1 - BLOOD[i];
        if (mask[i]) n += TYPECOUNT[i];
      }
      var lab = { all: 'All cancer types', solid: 'All solid tumors', blood: 'All blood cancers' }[val];
      var noun = { all: 'lines from all cancer types', solid: 'solid tumor lines', blood: 'blood cancer lines' }[val];
      var where = { all: 'across all cancer types', solid: 'across solid tumors', blood: 'across blood cancers' }[val];
      return { val: val, label: lab, pooled: true, mask: mask, noun: noun, where: where, nLines: n };
    }
    var k = +val.slice(2);
    mask[k] = 1;
    return { val: val, label: TYPES[k], pooled: false, mask: mask, type: k, noun: typeNoun(TYPES[k]) + ' lines', where: 'in ' + typeNoun(TYPES[k]), nLines: TYPECOUNT[k] };
  }
  function typeVal(name) { return 't:' + TYPES.indexOf(name); }
  function setVal(name) { return (name === 'all' || name === 'solid' || name === 'blood') ? name : typeVal(name); }

  // ------------------------------------------------------------------ core comparison
  function verdictOf(p, fold) {
    if (p < 0.05 && fold >= 2) return { code: 'sens' };
    if (p < 0.05 && fold <= 0.5) return { code: 'res' };
    if (p < 0.05) return { code: 'small' };
    return { code: 'none' };
  }

  function compare(drug, feat, S, opts) {
    opts = opts || {};
    var F = opts.filter || null;
    var off = drug.i * NC, rows = [], g1 = [], g0 = [], c, missAlt = 0;
    for (c = 0; c < NC; c++) {
      var v = IC[off + c];
      if (v !== v) {
        if (feat.mask[c] && S.mask[CELLTYPE[c]] && !(F && F.f.mask[c] !== F.keep)) missAlt++;
        continue;
      }
      var t = CELLTYPE[c];
      if (!S.mask[t]) continue;
      if (F && F.f.mask[c] !== F.keep) continue;
      var a16 = AUC16[off + c];
      var r = { c: c, v: v, a: a16 === 65535 ? NaN : a16 / 10000, alt: feat.mask[c], t: t };
      rows.push(r);
      (r.alt ? g1 : g0).push(r);
    }
    var R = { drug: drug, feat: feat, S: S, F: F, rows: rows, g1: g1, g0: g0, n1: g1.length, n0: g0.length, missAlt: missAlt };
    R.noun = S.noun + (F ? ' that ' + hasPlural(F.f, F.keep) : '');
    R.ext1 = g1.filter(function (r) { return r.v > drug.lnMax; }).length;
    R.ext0 = g0.filter(function (r) { return r.v > drug.lnMax; }).length;
    if (F && F.f === feat) { R.status = 'samefilter'; return R; }
    if (!rows.length) { R.status = 'untested'; return R; }
    if (!R.n1) { R.status = 'nomarker'; return R; }
    if (R.n1 < 3 || R.n0 < 3) { R.status = 'few'; return R; }
    R.status = 'ok';
    var x1 = g1.map(function (r) { return r.v; }), x0 = g0.map(function (r) { return r.v; });
    R.med1 = GS.median(x1); R.med0 = GS.median(x0);
    R.mw = GS.mannWhitney(x1, x0, opts.maxWork || 4e7);
    R.hl = GS.hodgesLehmann(x1, x0, opts.maxWork || 4e7);
    R.fold = Math.exp(R.hl.est);
    var a1 = g1.map(function (r) { return r.a; }).filter(isFinite), a0 = g0.map(function (r) { return r.a; }).filter(isFinite);
    R.aMed1 = GS.median(a1); R.aMed0 = GS.median(a0);
    R.mwA = (a1.length >= 3 && a0.length >= 3) ? GS.mannWhitney(a1, a0, opts.maxWork || 4e7) : null;
    if (S.pooled) {
      R.adj = GS.adjustedDiff(rows.map(function (r) { return r.v; }), rows.map(function (r) { return r.alt; }), rows.map(function (r) { return r.t; }));
    }
    if (S.pooled && R.adj.ok) {
      var b = R.adj.beta, se = R.adj.se;
      R.main = { fold: Math.exp(-b), p: R.adj.p, lo: Math.exp(-(b + 1.96 * se)), hi: Math.exp(-(b - 1.96 * se)), src: 'adj' };
    } else {
      R.main = { fold: R.fold, p: R.mw.p, lo: Math.exp(R.hl.lo), hi: Math.exp(R.hl.hi), src: S.pooled ? 'pooled' : 'mwu' };
    }
    R.verdict = verdictOf(R.main.p, R.main.fold);
    return R;
  }

  function acrossTypes(drug, feat, F) {
    var off = drug.i * NC, per = [], t, c;
    for (t = 0; t < NT; t++) per.push({ t: t, alt: [], ref: [] });
    for (c = 0; c < NC; c++) {
      var v = IC[off + c];
      if (v !== v) continue;
      if (F && F.f.mask[c] !== F.keep) continue;
      var row = { c: c, v: v };
      (feat.mask[c] ? per[CELLTYPE[c]].alt : per[CELLTYPE[c]].ref).push(row);
    }
    var out = per.filter(function (p) { return p.alt.length > 0; });
    out.forEach(function (p) {
      p.med1 = GS.median(p.alt.map(function (r) { return r.v; }));
      p.med0 = p.ref.length ? GS.median(p.ref.map(function (r) { return r.v; })) : NaN;
      if (p.alt.length >= 3 && p.ref.length >= 3) {
        var xa = p.alt.map(function (r) { return r.v; }), xr = p.ref.map(function (r) { return r.v; });
        p.fold = Math.exp(GS.hodgesLehmann(xa, xr, 0, true).est);
        p.p = GS.mannWhitney(xa, xr, 3e6).p;
      }
    });
    out.sort(function (a, b) { return ((b.alt.length >= 3) - (a.alt.length >= 3)) || (a.med1 - b.med1); });
    return out;
  }

  // The other marker that differs most between the two groups (Fisher's exact test).
  function coImbalance(R) {
    if (R.status !== 'ok' || R.S.pooled) return null;
    var best = null;
    FEATS.forEach(function (f2) {
      if (f2 === R.feat || f2.gene === R.feat.gene || (R.F && f2 === R.F.f)) return;
      var a = 0, c = 0, i;
      for (i = 0; i < R.n1; i++) a += f2.mask[R.g1[i].c];
      for (i = 0; i < R.n0; i++) c += f2.mask[R.g0[i].c];
      if (a + c < 3) return;
      var p = GS.fisherExact(a, R.n1 - a, c, R.n0 - c);
      if (!best || p < best.p) best = { f: f2, a: a, c: c, p: p };
    });
    return best && best.p < 0.001 ? best : null;
  }
  function twinAgreement(R, T) {
    if (R.status !== 'ok' || T.status !== 'ok') return 'na';
    var same = (R.main.fold >= 1) === (T.main.fold >= 1);
    if (R.verdict.code === T.verdict.code) return 'agree';
    if (R.main.p < 0.05 && T.main.p < 0.05 && !same) return 'disagree';
    if (same) return 'partial';
    return 'disagree';
  }
  function aucAgreement(R) {
    if (R.status !== 'ok' || !R.mwA) return 'na';
    var icSig = R.mw.p < 0.05, aSig = R.mwA.p < 0.05, same = (R.fold >= 1) === (R.aMed1 <= R.aMed0);
    if (icSig && aSig) return same ? 'agree' : 'disagree';
    if (!icSig && !aSig) return 'agree';
    return 'partial';
  }

  // ------------------------------------------------------------------ small SVG and icon helpers
  var ICON = {
    good: '<path d="M3.5 8.4l2.9 2.9 6.1-6.6"/>',
    warn: '<path d="M8 3.6v5.6"/><path d="M8 12.3v.2"/>',
    stop: '<path d="M8 3.6v5.6"/><path d="M8 12.3v.2"/>',
    no: '<path d="M4.6 4.6l6.8 6.8M11.4 4.6l-6.8 6.8"/>',
    dash: '<path d="M4 8h8"/>',
    info: '<path d="M8 7.2v5"/><path d="M8 4v.2"/>'
  };
  function svgIcon(k) { return '<svg viewBox="0 0 16 16" aria-hidden="true">' + ICON[k] + '</svg>'; }
  function iconDot(k, glyph) { return '<span class="ic k-' + k + '" aria-hidden="true">' + svgIcon(glyph || k) + '</span>'; }
  function L(x1, y1, x2, y2, cls) { return '<line x1="' + x1.toFixed(1) + '" y1="' + y1.toFixed(1) + '" x2="' + x2.toFixed(1) + '" y2="' + y2.toFixed(1) + '" class="' + cls + '"/>'; }
  function RECT(x, y, w, h, cls, rx) { return '<rect x="' + x.toFixed(1) + '" y="' + y.toFixed(1) + '" width="' + Math.max(0, w).toFixed(1) + '" height="' + Math.max(0, h).toFixed(1) + '"' + (rx ? ' rx="' + rx + '"' : '') + ' class="' + cls + '"/>'; }
  function T(x, y, str, anchor, cls) { return '<text x="' + x.toFixed(1) + '" y="' + y.toFixed(1) + '"' + (anchor && anchor !== 'start' ? ' text-anchor="' + anchor + '"' : '') + ' class="' + cls + '">' + esc(str) + '</text>'; }
  function logTicks(lo, hi) {
    var a = Math.floor(lo), b = Math.ceil(hi), t = [], e;
    var step = (b - a) > 8 ? 2 : 1;
    for (e = a; e <= b; e += step) t.push(e);
    if (b - a <= 2) {
      var extra = [];
      for (e = a - 1; e <= b; e++) [2, 5].forEach(function (m) { var v = e + Math.log10(m); if (v > lo && v < hi) extra.push(v); });
      t = t.concat(extra).sort(function (x, y) { return x - y; });
    }
    return t.filter(function (v) { return v >= lo - 1e-9 && v <= hi + 1e-9; });
  }
  function fmtTick(l10) {
    var v = Math.pow(10, l10);
    if (v >= 1) return nf(Math.round(v));
    return String(Number(v.toPrecision(1)));
  }
  function jitter(i) {
    var x = Math.sin(i * 12.9898 + 78.233) * 43758.5453;
    return x - Math.floor(x) - 0.5;
  }
  function wrapWords(str, max) {
    var words = String(str).split(' '), lines = [], cur = '';
    words.forEach(function (w) {
      if (!cur) cur = w;
      else if ((cur + ' ' + w).length <= max) cur += ' ' + w;
      else { lines.push(cur); cur = w; }
    });
    if (cur) lines.push(cur);
    return lines;
  }
  function niceDomain(mn, mx, pad) {
    var lo = mn - pad, hi = mx + pad, fl = Math.floor(mn), ce = Math.ceil(mx);
    if (mn - fl < 0.45) lo = Math.min(lo, fl - 0.08);
    if (ce - mx < 0.45) hi = Math.max(hi, ce + 0.08);
    if (hi - lo < 2) { var mid = (hi + lo) / 2; lo = mid - 1; hi = mid + 1; }
    return [lo, hi];
  }
  function boxWidth(id) { return Math.max(300, Math.min(920, Math.round($(id).clientWidth || 700))); }
  function svgPoint(svg, W, H, e) {
    var r = svg.getBoundingClientRect();
    return { x: (e.clientX - r.left) * W / r.width, y: (e.clientY - r.top) * H / r.height };
  }
  function nearest(points, p, maxD) {
    var best = null, bd = maxD * maxD;
    for (var i = 0; i < points.length; i++) {
      var dx = points[i].x - p.x, dy = points[i].y - p.y, d2 = dx * dx + dy * dy;
      if (d2 < bd) { bd = d2; best = points[i]; }
    }
    return best;
  }
  function tipShow(tip, box, W, x, y, html) {
    var sc = box.clientWidth / W, left = x * sc;
    left = Math.max(110, Math.min(box.clientWidth - 110, left));
    tip.innerHTML = html;
    tip.style.left = left + 'px';
    tip.style.top = (y * sc) + 'px';
    tip.hidden = false;
  }

  // ------------------------------------------------------------------ pickers
  var KIND_ORDER = ['mut', 'amp', 'del', 'fusion', 'msi'];
  var KIND_GROUP = { mut: 'Mutations', amp: 'Extra copies (amplifications)', del: 'Lost copies (deletions)', fusion: 'Fused genes', msi: 'DNA repair' };
  function shortText(s, n) { return s.length > n ? s.slice(0, n - 1).replace(/[\s,]+$/, '') + '…' : s; }
  var DRUG_ITEMS = DRUGS.slice().sort(function (a, b) {
    return a.name.toLowerCase().localeCompare(b.name.toLowerCase()) || a.id - b.id;
  }).map(function (d) {
    var screen = d.twin.length ? ' (screen ' + d.id + ')' : '';
    return { value: d.id, main: d.name, sub: shortText(d.target, 26) + screen, meta: shortText(d.target, 40) + screen, subLong: 'Targets ' + d.target + screen, hay: d.hay };
  });
  var GENE_ITEMS = FEATS.slice().sort(function (a, b) {
    return KIND_ORDER.indexOf(a.kind) - KIND_ORDER.indexOf(b.kind) || a.label.localeCompare(b.label);
  }).map(function (f) {
    return { value: f.key, main: f.label, sub: f.cells.length + ' lines', subLong: f.cells.length + ' of ' + nf(NC) + ' cell lines carry it', group: KIND_GROUP[f.kind], hay: f.hay };
  });
  var TYPE_ITEMS = (function () {
    var nBlood = 0, nSolid = 0, i;
    for (i = 0; i < NT; i++) { if (BLOOD[i]) nBlood += TYPECOUNT[i]; else nSolid += TYPECOUNT[i]; }
    var out = [
      { value: 'all', main: 'All cancer types', sub: nf(NC) + ' lines', subLong: nf(NC) + ' cell lines, compared within each type', group: 'Several types at once', hay: 'all cancer types pooled' },
      { value: 'solid', main: 'All solid tumors', sub: nf(nSolid) + ' lines', subLong: nf(nSolid) + ' cell lines, compared within each type', group: 'Several types at once', hay: 'all solid tumors pooled' },
      { value: 'blood', main: 'All blood cancers', sub: nf(nBlood) + ' lines', subLong: nf(nBlood) + ' cell lines, compared within each type', group: 'Several types at once', hay: 'all blood cancers leukemia lymphoma myeloma pooled' }
    ];
    TYPES.map(function (t, k) { return { t: t, k: k }; }).sort(function (a, b) { return a.t.localeCompare(b.t); }).forEach(function (o) {
      out.push({ value: 't:' + o.k, main: o.t, sub: TYPECOUNT[o.k] + ' lines', subLong: TYPECOUNT[o.k] + ' cell lines', group: 'One cancer type', hay: o.t.toLowerCase() });
    });
    return out;
  })();

  var PICKERS = {};
  // statusFn (optional) returns { value: { ok, txt } } for the current question, so each option can show
  // whether that choice leaves at least 3 cell lines on each side
  var ASK_LABELS = { yes: 'Enough lines to answer', no: 'Too few', only: 'Show only answerable', groupYes: 'Enough lines to answer', groupNo: 'Too few lines to answer', srYes: 'Enough lines. ', srNo: 'Too few lines. ', none: 'None of these has enough cell lines with your other two choices. Clear Show only answerable to see them all.' };
  var RANK_LABELS = { yes: 'Can be ranked', no: 'Nothing to rank', only: 'Show only ones with results', groupYes: 'Can be ranked', groupNo: 'Nothing to rank', srYes: 'Can be ranked. ', srNo: 'Nothing to rank. ', none: 'None of these has anything to rank with your other choice. Clear Show only ones with results to see them all.' };
  function Picker(key, items, onPick, statusFn, labels) {
    labels = labels || ASK_LABELS;
    var root = $('pk-' + key), btn = $('pk-' + key + '-btn'), pop = $('pk-' + key + '-pop'), q = $('pk-' + key + '-q');
    var list = $('pk-' + key + '-list'), empty = $('pk-' + key + '-empty'), emptyText = empty.textContent;
    var byVal = {}, statuses = null, onlyBox = null;
    items.forEach(function (it) { byVal[it.value] = it; });
    var P = { value: null, shown: [], active: -1, root: root };
    if (statusFn) {
      var tools = document.createElement('div');
      tools.className = 'pk-tools';
      tools.innerHTML = '<p class="pk-key"><span>' + iconDot('good') + esc(labels.yes) + '</span><span>' + iconDot('stop', 'no') + esc(labels.no) + '</span></p>' +
        '<label class="pk-only"><input type="checkbox" id="pk-' + key + '-only">' + esc(labels.only) + '</label>';
      q.insertAdjacentElement('afterend', tools);
      onlyBox = $('pk-' + key + '-only');
      onlyBox.addEventListener('change', function () { render(); q.focus(); });
    }
    P.set = function (value) {
      var it = byVal[value];
      if (!it) return;
      P.value = value;
      $('pk-' + key + '-main').textContent = it.main;
      $('pk-' + key + '-sub').textContent = it.subLong || it.sub || '';
    };
    function setActive(i, scroll) {
      var prev = list.querySelector('.pk-opt.active');
      if (prev) prev.classList.remove('active');
      P.active = i;
      if (i < 0) { q.removeAttribute('aria-activedescendant'); return; }
      var el = $('pk-' + key + '-o' + i);
      if (!el) return;
      el.classList.add('active');
      q.setAttribute('aria-activedescendant', el.id);
      if (scroll) el.scrollIntoView({ block: 'nearest' });
    }
    function render() {
      var words = q.value.trim().toLowerCase().split(/\s+/).filter(Boolean), html = '', lastGroup = null, sel = -1;
      var only = onlyBox && onlyBox.checked;
      P.shown = [];
      var matched = items.filter(function (it) { return !words.length || words.every(function (w) { return it.hay.indexOf(w) > -1; }); });
      if (statuses) {
        // answerable choices first, under their own heading, then the rest
        var okIt = matched.filter(function (it) { return statuses[it.value].ok; }), noIt = matched.filter(function (it) { return !statuses[it.value].ok; });
        matched = okIt.map(function (it) { return { it: it, g: labels.groupYes + ' (' + okIt.length + ')' }; })
          .concat(only ? [] : noIt.map(function (it) { return { it: it, g: labels.groupNo + ' (' + noIt.length + ')' }; }));
      } else {
        matched = matched.map(function (it) { return { it: it, g: it.group }; });
      }
      matched.forEach(function (m) {
        var it = m.it, st = statuses ? statuses[it.value] : null;
        if (m.g && m.g !== lastGroup) { html += '<li class="pk-group" role="presentation">' + esc(m.g) + '</li>'; lastGroup = m.g; }
        var idx = P.shown.length;
        P.shown.push(it);
        if (it.value === P.value) sel = idx;
        var right = st
          ? '<span class="pk-st ' + (st.ok ? 'st-ok' : 'st-no') + '">' + iconDot(st.ok ? 'good' : 'stop', st.ok ? 'good' : 'no') + '<span class="sr">' + (st.ok ? labels.srYes : labels.srNo) + '</span>' + esc(st.txt) + '</span>'
          : '<span class="pk-opt-sub">' + esc(it.sub) + '</span>';
        var meta = st && it.meta ? '<span class="pk-opt-meta">' + esc(it.meta) + '</span>' : '';
        html += '<li class="pk-opt' + (st && !st.ok ? ' no' : '') + '" role="option" id="pk-' + key + '-o' + idx + '" data-i="' + idx + '" aria-selected="' + (it.value === P.value) + '">' +
          '<span class="pk-opt-text"><span class="pk-opt-main">' + esc(it.main) + '</span>' + meta + '</span>' + right + '</li>';
      });
      list.innerHTML = html;
      empty.textContent = (only && !words.length) ? labels.none : emptyText;
      empty.hidden = P.shown.length > 0;
      setActive(sel >= 0 ? sel : (P.shown.length ? 0 : -1), true);
    }
    function open() {
      closeAll(key);
      if (statusFn) statuses = statusFn();
      pop.hidden = false;
      btn.setAttribute('aria-expanded', 'true');
      q.value = '';
      render();
      q.focus();
    }
    P.close = function (refocus) {
      if (pop.hidden) return;
      pop.hidden = true;
      btn.setAttribute('aria-expanded', 'false');
      if (refocus) btn.focus();
    };
    function choose(i) {
      var it = P.shown[i];
      if (!it) return;
      P.close(true);
      if (it.value !== P.value) { P.set(it.value); onPick(it.value); }
    }
    btn.addEventListener('click', function () { if (pop.hidden) open(); else P.close(false); });
    btn.addEventListener('keydown', function (e) {
      if (e.key === 'ArrowDown' || e.key === 'ArrowUp') { e.preventDefault(); open(); }
    });
    q.addEventListener('input', render);
    q.addEventListener('keydown', function (e) {
      if (e.key === 'ArrowDown') { e.preventDefault(); if (P.shown.length) setActive(Math.min(P.shown.length - 1, P.active + 1), true); }
      else if (e.key === 'ArrowUp') { e.preventDefault(); if (P.shown.length) setActive(Math.max(0, P.active - 1), true); }
      else if (e.key === 'PageDown') { e.preventDefault(); if (P.shown.length) setActive(Math.min(P.shown.length - 1, P.active + 8), true); }
      else if (e.key === 'PageUp') { e.preventDefault(); if (P.shown.length) setActive(Math.max(0, P.active - 8), true); }
      else if (e.key === 'Enter') { e.preventDefault(); if (P.active >= 0) choose(P.active); }
      else if (e.key === 'Escape') { e.preventDefault(); P.close(true); }
    });
    root.addEventListener('focusout', function (e) {
      if (!pop.hidden && e.relatedTarget && !root.contains(e.relatedTarget)) P.close(false);
    });
    pop.addEventListener('keydown', function (e) { if (e.key === 'Escape') { e.preventDefault(); P.close(true); } });
    list.addEventListener('mousedown', function (e) { e.preventDefault(); });
    list.addEventListener('click', function (e) {
      var o = e.target.closest('.pk-opt');
      if (o) choose(+o.getAttribute('data-i'));
    });
    list.addEventListener('mousemove', function (e) {
      var o = e.target.closest('.pk-opt');
      if (!o) return;
      var i = +o.getAttribute('data-i');
      if (i !== P.active) setActive(i, false);
    });
    PICKERS[key] = P;
    return P;
  }
  function closeAll(except) {
    Object.keys(PICKERS).forEach(function (k) { if (k !== except) PICKERS[k].close(false); });
  }
  document.addEventListener('pointerdown', function (e) {
    Object.keys(PICKERS).forEach(function (k) { if (!PICKERS[k].root.contains(e.target)) PICKERS[k].close(false); });
  });

  // ------------------------------------------------------------------ which choices can be answered
  var TESTED = {};
  function testedCells(d) {
    if (!TESTED[d.id]) {
      var off = d.i * NC, a = [];
      for (var c = 0; c < NC; c++) if (IC[off + c] === IC[off + c]) a.push(c);
      TESTED[d.id] = a;
    }
    return TESTED[d.id];
  }
  function stat(n1, n0) { return { ok: n1 >= 3 && n0 >= 3, txt: (n1 + n0) ? n1 + ' with, ' + n0 + ' without' : 'Not tested' }; }
  function drugStatuses() {
    var f = FEAT_BY_KEY[ST.gene], S = setInfo(ST.set), out = {};
    DRUGS.forEach(function (d) {
      var n1 = 0, n0 = 0;
      testedCells(d).forEach(function (c) { if (S.mask[CELLTYPE[c]]) { if (f.mask[c]) n1++; else n0++; } });
      out[d.id] = stat(n1, n0);
    });
    return out;
  }
  function geneStatuses() {
    var S = setInfo(ST.set), out = {};
    var cells = testedCells(DRUG_BY_ID[ST.drug]).filter(function (c) { return S.mask[CELLTYPE[c]]; });
    FEATS.forEach(function (f) {
      var n1 = 0;
      cells.forEach(function (c) { n1 += f.mask[c]; });
      out[f.key] = stat(n1, cells.length - n1);
    });
    return out;
  }
  function typeStatuses() {
    var f = FEAT_BY_KEY[ST.gene], n1 = [], n0 = [], out = {}, t;
    for (t = 0; t < NT; t++) { n1.push(0); n0.push(0); }
    testedCells(DRUG_BY_ID[ST.drug]).forEach(function (c) { if (f.mask[c]) n1[CELLTYPE[c]]++; else n0[CELLTYPE[c]]++; });
    var sums = { all: [0, 0], solid: [0, 0], blood: [0, 0] };
    for (t = 0; t < NT; t++) {
      out['t:' + t] = stat(n1[t], n0[t]);
      sums.all[0] += n1[t]; sums.all[1] += n0[t];
      var k = BLOOD[t] ? 'blood' : 'solid';
      sums[k][0] += n1[t]; sums[k][1] += n0[t];
    }
    ['all', 'solid', 'blood'].forEach(function (k) { out[k] = stat(sums[k][0], sums[k][1]); });
    return out;
  }

  // which rank choices have anything to rank: mirrors compareLite() exactly, including the
  // within-type model's needs for pooled sets (one cancer type with lines on both sides, and df >= 1)
  var NAME_IDX = {}, NNAMES = 0;
  DRUGS.forEach(function (d) { if (!(d.name in NAME_IDX)) NAME_IDX[d.name] = NNAMES++; });
  var ELIG = {};
  function eligFor(setVal) {
    if (ELIG[setVal]) return ELIG[setVal];
    var S = setInfo(setVal), ND = DRUGS.length, di, t;
    var nS = new Int32Array(ND), K = new Int32Array(ND), nT = new Int32Array(ND * NT);
    DRUGS.forEach(function (d, i) {
      testedCells(d).forEach(function (c) { var tt = CELLTYPE[c]; if (S.mask[tt]) { nS[i]++; nT[i * NT + tt]++; } });
      for (var u = 0; u < NT; u++) if (nT[i * NT + u]) K[i]++;
    });
    var byGene = {}, byDrug = new Int32Array(ND), alt = new Int32Array(ND * NT), n1 = new Int32Array(ND), seen = new Uint8Array(NNAMES);
    FEATS.forEach(function (f) {
      n1.fill(0); seen.fill(0);
      var touched = [];
      f.cells.forEach(function (c) {
        var tt = CELLTYPE[c];
        if (!S.mask[tt]) return;
        if (touched.indexOf(tt) < 0) touched.push(tt);
        for (var j = 0; j < ND; j++) { var v = IC[j * NC + c]; if (v === v) { n1[j]++; alt[j * NT + tt]++; } }
      });
      var cnt = 0;
      for (di = 0; di < ND; di++) {
        var a = n1[di], ok = a >= 3 && nS[di] - a >= 3;
        if (ok && S.pooled) {
          var inf = false;
          for (var k = 0; k < touched.length; k++) { var x = alt[di * NT + touched[k]]; if (x >= 1 && nT[di * NT + touched[k]] - x >= 1) { inf = true; break; } }
          ok = inf && nS[di] - K[di] - 1 >= 1;
        }
        if (ok) {
          byDrug[di]++;
          var ni = NAME_IDX[DRUGS[di].name];
          if (!seen[ni]) { seen[ni] = 1; cnt++; }
        }
      }
      touched.forEach(function (u) { for (var j = 0; j < ND; j++) alt[j * NT + u] = 0; });
      byGene[f.key] = cnt;
    });
    var byDrugId = {};
    DRUGS.forEach(function (d, i) { byDrugId[d.id] = byDrug[i]; });
    ELIG[setVal] = { byGene: byGene, byDrug: byDrugId };
    return ELIG[setVal];
  }
  function stat2(n, noun) { return { ok: n >= 1, txt: n ? n + ' ' + plural(n, noun) : 'None' }; }
  function rankGeneStatuses() {
    var e = eligFor(RK.set).byGene, out = {};
    FEATS.forEach(function (f) { out[f.key] = stat2(e[f.key], 'drug'); });
    return out;
  }
  function rankDrugStatuses() {
    var e = eligFor(RK.set).byDrug, out = {};
    DRUGS.forEach(function (d) { out[d.id] = stat2(e[d.id], 'gene change'); });
    return out;
  }
  function rankTypeStatuses() {
    var t, sets = ['all', 'solid', 'blood'], members = { all: [], solid: [], blood: [] }, counts = {}, names = {};
    for (t = 0; t < NT; t++) { members.all.push(t); members[BLOOD[t] ? 'blood' : 'solid'].push(t); counts['t:' + t] = 0; names['t:' + t] = {}; }
    sets.forEach(function (k) { counts[k] = 0; names[k] = {}; });
    function poolOk(n1t, n0t, mem) {
      var a = 0, b = 0, K = 0, inf = false;
      mem.forEach(function (u) { a += n1t[u]; b += n0t[u]; if (n1t[u] + n0t[u]) K++; if (n1t[u] >= 1 && n0t[u] >= 1) inf = true; });
      return a >= 3 && b >= 3 && inf && a + b - K - 1 >= 1;
    }
    function add(n1t, n0t, label) {
      for (var u = 0; u < NT; u++) if (n1t[u] >= 3 && n0t[u] >= 3) names['t:' + u][label] = 1;
      sets.forEach(function (k) { if (poolOk(n1t, n0t, members[k])) names[k][label] = 1; });
    }
    if (RK.mode === 'drugs') {
      var f = FEAT_BY_KEY[RK.gene];
      DRUGS.forEach(function (d) {
        var n1t = new Int32Array(NT), n0t = new Int32Array(NT);
        testedCells(d).forEach(function (c) { if (f.mask[c]) n1t[CELLTYPE[c]]++; else n0t[CELLTYPE[c]]++; });
        add(n1t, n0t, d.name);
      });
    } else {
      var d = DRUG_BY_ID[RK.drug], nT = new Int32Array(NT), off = d.i * NC;
      testedCells(d).forEach(function (c) { nT[CELLTYPE[c]]++; });
      FEATS.forEach(function (f) {
        var n1t = new Int32Array(NT), n0t = new Int32Array(NT);
        f.cells.forEach(function (c) { if (IC[off + c] === IC[off + c]) n1t[CELLTYPE[c]]++; });
        for (var u = 0; u < NT; u++) n0t[u] = nT[u] - n1t[u];
        add(n1t, n0t, f.key);
      });
    }
    var out = {}, noun = RK.mode === 'drugs' ? 'drug' : 'gene change';
    Object.keys(names).forEach(function (k) { out[k] = stat2(Object.keys(names[k]).length, noun); });
    return out;
  }

  // ------------------------------------------------------------------ state and run
  var ST = { drug: 1036, gene: 'BRAF_mut', set: typeVal('Melanoma'), filter: null, last: null, typesAll: false };

  function run() {
    var d = DRUG_BY_ID[ST.drug], f = FEAT_BY_KEY[ST.gene], S = setInfo(ST.set), F = ST.filter;
    var R = compare(d, f, S, { filter: F });
    R.twins = d.twin.map(function (id) { return compare(DRUG_BY_ID[id], f, S, { filter: F }); });
    R.across = acrossTypes(d, f, F);
    R.co = coImbalance(R);
    ST.last = R;
    renderFilterNote(R);
    markExamples();
    renderAnswer(R);
    renderDose();
    renderChecks(R);
    renderTypes();
    renderDetails(R);
  }
  function setQuestion(drug, gene, set) {
    ST.drug = drug; ST.gene = gene; ST.set = set; ST.filter = null; ST.typesAll = false;
    PICKERS.drug.set(drug); PICKERS.gene.set(gene); PICKERS.type.set(set);
    run();
  }

  function renderFilterNote(R) {
    var el = $('filterNote');
    if (!R.F) { el.hidden = true; el.innerHTML = ''; return; }
    el.hidden = false;
    el.innerHTML = '<span>Leaving out ' + esc(R.F.f.altName) + ' from both groups.</span><button type="button" class="btn" data-undo="1">Show all lines again</button>';
  }

  // ------------------------------------------------------------------ the answer
  var BADGE = {
    sens: ['good', 'good', LABEL.sens], res: ['good', 'good', LABEL.res], small: ['warn', 'warn', LABEL.small], none: ['no', 'dash', LABEL.none],
    nomarker: ['stop', 'stop', 'Not enough cell lines'], few: ['stop', 'stop', 'Not enough cell lines'],
    untested: ['stop', 'stop', 'Drug not tested here'], samefilter: ['stop', 'stop', 'Check your question']
  };
  function askedText(R) {
    var d = drugText(R.drug), S = R.S, g = groupNoun(R.feat, S);
    var s = (S.pooled ? cap(S.where) + ', do ' + g : 'Do ' + g) + ' respond differently to ' + d + ' than ' + otherNoun(S) + '?';
    if (R.F) s += ' Lines that ' + hasPlural(R.F.f, 1) + ' are left out.';
    return s;
  }
  function renderAnswer(R) {
    var d = drugText(R.drug), f = R.feat, S = R.S, head = '', basis = [], code = R.status === 'ok' ? R.verdict.code : R.status;
    var b = BADGE[code];
    $('asked').innerHTML = '<b>Your question.</b> ' + nw(askedText(R));
    var badge = $('badge');
    badge.className = 'badge k-' + b[0];
    badge.innerHTML = svgIcon(b[1]) + esc(b[2]);
    if (R.status === 'ok') {
      var m = R.main, x = foldX(m);
      if (code === 'none') {
        head = cap(f.statusName) + ' made no clear difference to how ' + (S.pooled ? 'lines responded to ' + d + ' within each cancer type' : S.noun + ' responded to ' + d) + '.';
      } else {
        head = cap(groupNoun(f, S)) + ' needed ' + fmtFold(x) + ' times ' + lessMore(m) + ' ' + d + ' than ' + otherNoun(S) + '.';
        if (code === 'small') head += ' That gap is too small to count as a marker.';
      }
      basis.push('The tool compared ' + nLines(R.n1, f, 1) + ' against ' + nLines(R.n0, f, 0) + (S.pooled ? ', comparing lines only within the same cancer type' : '') + (R.F ? ', after leaving out ' + R.F.f.altName : '') + '.');
    } else if (R.status === 'nomarker') {
      head = 'None of the ' + R.n0 + ' ' + R.noun + ' tested with ' + d + ' ' + f.has + '.';
      if (R.missAlt) basis.push('The dataset has ' + missText(R) + ', but ' + (R.missAlt === 1 ? 'it has no' : 'none of them has a') + ' usable result for ' + d + '. GDSC did not test every drug on every cell line, so this empty group says nothing about how those lines would respond.');
      basis.push('With one group empty there is nothing to compare, so the tool gives no answer. Try All cancer types, or another cancer type.');
    } else if (R.status === 'few') {
      var small1 = R.n1 < 3;
      var nSmall = small1 ? R.n1 : R.n0;
      head = 'Only ' + (small1 ? nLines(R.n1, f, 1) : nLines(R.n0, f, 0)) + ' ' + S.where + ' ' + (R.F ? (nSmall === 1 ? 'is' : 'are') + ' left after leaving out ' + R.F.f.altName : (nSmall === 1 ? 'was' : 'were') + ' tested with ' + d) + '.';
      basis.push('The tool needs at least 3 cell lines in each group before it gives an answer, because one unusual line could decide the result.');
      if (small1 && R.missAlt) basis.push('The dataset holds ' + R.missAlt + ' more ' + missText(R).replace(/^\d+ /, '') + ' with no usable result for ' + d + '. GDSC did not test every drug on every cell line.');
    } else if (R.status === 'untested') {
      head = cap(d) + ' was not tested on any ' + R.noun + ' in this dataset.';
      basis.push(R.F ? 'Leaving out ' + R.F.f.altName + ' may have removed every line. Show all lines again to check.' : 'Pick another cancer type, or choose All cancer types.');
    } else {
      head = 'The gene change and the lines left out are the same.';
      basis.push('Show all lines again, or pick a different gene change.');
    }
    $('headline').innerHTML = nw(head);
    $('basis').innerHTML = basis.map(function (t) { return '<p>' + nw(t) + '</p>'; }).join('');
  }

  // ------------------------------------------------------------------ dose chart
  var DOSE = { pts: [], W: 700, H: 300 };
  function renderDose() {
    var R = ST.last, fig = $('doseFig'), svg = $('dose');
    $('doseTip').hidden = true;
    if (R.status === 'untested' || R.status === 'samefilter' || !R.rows.length) { fig.hidden = true; return; }
    fig.hidden = false;
    var d = R.drug, f = R.feat, S = R.S;
    $('doseTitle').textContent = 'How much ' + drugText(d) + ' each ' + (S.pooled ? '' : typeNoun(TYPES[S.type]) + ' ') + 'cell line needed';
    $('doseLegend').innerHTML =
      '<li><span class="sw sw-gene"></span>' + esc(f.altShort) + '</li>' +
      '<li><span class="sw sw-no"></span>' + esc(f.refShort) + '</li>' +
      '<li><span class="sw sw-est"></span>Estimated, past the highest dose</li>' +
      '<li><span class="sw sw-med"></span>Median</li>';
    var W = boxWidth('doseBox');
    var groups = [{ rows: R.g1, cls: 'gene', name: f.altShort }, { rows: R.g0, cls: 'no', name: f.refShort }];
    var longest = Math.max(f.altShort.length, f.refShort.length);
    var stacked = W < 560 || longest > 20;
    var LW = stacked ? 0 : 172, PR = 14, TOP = 34, AXH = 70;
    var fs = W < 420 ? 13 : 14;
    groups.forEach(function (g) {
      var n = g.rows.length;
      g.n = n;
      g.med = n ? GS.median(g.rows.map(function (r) { return r.v; })) : NaN;
      g.nTxt = n + ' ' + plural(n, 'line') + (n > 0 && n < 3 ? ', too few to compare' : '');
      g.medTxt = n ? 'median ' + fmtConc(Math.exp(g.med)) + ' µM' : '';
      g.nameLines = stacked ? wrapWords(g.name, Math.max(18, Math.floor((W - 20) / (fs * 0.62)))) : [g.name];
      g.labH = stacked ? g.nameLines.length * 18 + 28 : 0;
      g.rh = stacked ? g.labH + 62 : 78;
    });
    var H = TOP + groups[0].rh + groups[1].rh + AXH;
    var x0 = stacked ? 10 : LW + 10, x1 = W - PR;
    var vals = R.rows.map(function (r) { return r.v / LN10; }), lmax = d.lnMax / LN10;
    var dom = niceDomain(Math.min(Math.min.apply(null, vals), lmax), Math.max(Math.max.apply(null, vals), lmax), 0.25), lo = dom[0], hi = dom[1];
    function X(l) { return x0 + (l - lo) / (hi - lo) * (x1 - x0); }
    DOSE.W = W; DOSE.H = H; DOSE.pts = [];
    svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    var s = '', yTop = TOP - 4, yBot = TOP + groups[0].rh + groups[1].rh, xm = X(lmax);
    s += RECT(xm, yTop, x1 - xm, yBot - yTop, 'zone');
    logTicks(lo, hi).forEach(function (t) {
      var x = X(t);
      s += L(x, yTop, x, yBot, 'grid');
      s += T(x, yBot + 18, fmtTick(t), 'middle', 't');
    });
    s += L(xm, yTop - 10, xm, yBot, 'edge');
    var lab = 'Highest dose tested, ' + fmtConc(d.maxc) + ' µM', labW = lab.length * 6.6;
    if (xm - x0 >= labW + 8) s += T(xm - 6, TOP - 16, lab, 'end', 't t-2');
    else if (x1 - xm >= labW + 8) s += T(xm + 6, TOP - 16, lab, 'start', 't t-2');
    else s += T(x0, TOP - 16, lab, 'start', 't t-2');
    var top = TOP;
    groups.forEach(function (g, gi) {
      if (gi === 1) s += L(0, top, W, top, 'grid');
      var bandTop = top + (stacked ? g.labH : 10), bandBot = top + g.rh - 10, cy = (bandTop + bandBot) / 2, amp = Math.max(4, (bandBot - bandTop) / 2 - 6);
      if (stacked) {
        g.nameLines.forEach(function (ln, i) { s += T(x0, top + 18 + i * 18, ln, 'start', 't t-l t-ink t-b'); });
        s += T(x0, top + 18 + g.nameLines.length * 18, g.nTxt + (g.medTxt ? ', ' + g.medTxt : ''), 'start', 't');
      } else {
        s += T(0, cy - 8, g.name, 'start', 't t-l t-ink t-b');
        s += T(0, cy + 10, g.nTxt, 'start', 't');
        if (g.medTxt) s += T(0, cy + 27, g.medTxt, 'start', 't');
      }
      var n = g.n, rad = n > 300 ? 3 : n > 120 ? 3.8 : 5;
      g.rows.slice().sort(function (a, b) { return (b.v > d.lnMax) - (a.v > d.lnMax); }).forEach(function (r) {
        var x = X(r.v / LN10), y = cy + jitter(r.c) * 2 * amp, est = r.v > d.lnMax;
        s += '<circle class="dot dot-' + g.cls + (est ? ' est' : '') + '" data-i="' + DOSE.pts.length + '" cx="' + x.toFixed(1) + '" cy="' + y.toFixed(1) + '" r="' + rad + '"/>';
        DOSE.pts.push({ x: x, y: y, r: r, g: gi });
      });
      if (n) { var xmed = X(g.med / LN10); s += L(xmed, cy - amp - 6, xmed, cy + amp + 6, 'med'); }
      else s += T((x0 + x1) / 2, cy + 5, 'No cell lines in this group', 'middle', 't t-l');
      top += g.rh;
    });
    s += L(x0, yBot, x1, yBot, 'axis');
    s += T((x0 + x1) / 2, yBot + 40, W < 520 ? 'IC50 (µM, log scale)' : 'IC50, the dose that halves growth (µM, log scale)', 'middle', 't t-2');
    s += T(x0, yBot + 62, '← less drug needed', 'start', 't');
    s += T(x1, yBot + 62, 'more drug needed →', 'end', 't');
    svg.innerHTML = s;
    var notes = [], ext = R.ext1 + R.ext0;
    if (ext) notes.push(ext + ' of ' + R.rows.length + ' lines never dropped to half their growth, even at ' + fmtConc(d.maxc) + ' µM, so their IC50 is an estimate.');
    if (S.pooled) notes.push('Lines from different cancer types are shown together here. The answer compares lines only within the same cancer type.');
    if (R.status === 'nomarker' || R.status === 'few') notes.push('Each group needs at least 3 lines before the tool gives an answer.');
    $('doseNote').textContent = notes.join(' ');
  }
  function doseHover(e) {
    var svg = $('dose'), tip = $('doseTip'), p = svgPoint(svg, DOSE.W, DOSE.H, e), best = nearest(DOSE.pts, p, 18);
    var prev = svg.querySelector('.dot.hl');
    if (prev) prev.classList.remove('hl');
    if (!best) { tip.hidden = true; return; }
    var el = svg.querySelector('[data-i="' + DOSE.pts.indexOf(best) + '"]');
    if (el) el.classList.add('hl');
    var R = ST.last, r = best.r, est = r.v > R.drug.lnMax;
    var h = '<b>' + esc(CELLNAME[r.c]) + '</b><br>IC50 ' + (est ? 'about ' : '') + fmtConc(Math.exp(r.v)) + ' µM' + (est ? ', estimated' : '') +
      '<br>' + esc(r.alt ? R.feat.altShort : R.feat.refShort) + (R.S.pooled ? '<br>' + esc(TYPES[r.t]) : '');
    tipShow(tip, $('doseBox'), DOSE.W, best.x, best.y, h);
  }

  // ------------------------------------------------------------------ checks
  function renderChecks(R) {
    var block = $('checksBlock');
    if (R.status !== 'ok') { block.hidden = true; return; }
    block.hidden = false;
    var f = R.feat, S = R.S, m = R.main, x = foldX(m), dn = drugText(R.drug), items = [];
    if (x >= 2) items.push({ k: 'good', t: 'Big enough to matter', d: 'Lines with the gene change needed ' + fmtFold(x) + ' times ' + lessMore(m) + ' drug. The tool asks for a gap of at least 2 times.' });
    else items.push({ k: 'no', t: 'Too small to matter', d: 'The gap is ' + fmtFold(x) + ' times, under the 2 times the tool asks for.' });
    if (m.p < 0.05) items.push({ k: 'good', t: 'Unlikely to be chance', d: pEq(m.p) + ', below the 0.05 cutoff. ' + chanceWords(m.p) });
    else items.push({ k: 'no', t: 'Could be chance', d: pEq(m.p) + ', above the 0.05 cutoff. ' + chanceWords(m.p) });
    var mn = Math.min(R.n1, R.n0);
    if (mn >= 5) items.push({ k: 'good', t: 'Enough cell lines', d: cap(nLines(R.n1, f, 1)) + ' and ' + nLines(R.n0, f, 0) + '.' });
    else items.push({ k: 'warn', t: 'A small group', d: 'The smaller group has only ' + mn + ' lines, so one unusual line could change the answer.' });
    if (S.pooled) {
      if (R.adj.ok) items.push({ k: 'good', t: 'Compared within each cancer type', d: 'Cancer types differ in how they respond to many drugs, so the tool compares lines only within the same type. ' + R.adj.informative + ' cancer ' + plural(R.adj.informative, 'type has', 'types have') + ' lines on both sides.' });
      else items.push({ k: 'warn', t: 'Could not compare within cancer types', d: 'The gene change never varies within a cancer type here, so cancer type may explain part of the gap.' });
    }
    R.twins.forEach(function (Tw) {
      var a = twinAgreement(R, Tw), tx = Tw.status === 'ok' ? fmtFold(foldX(Tw.main)) + ' times ' + lessMore(Tw.main) + ' drug, ' + pEq(Tw.main.p) : '';
      if (a === 'agree') items.push({ k: 'good', t: 'A second screen agrees', d: cap(dn) + ' was screened twice. The other screen (ID ' + Tw.drug.id + ') gave the same answer, ' + tx + '.' });
      else if (a === 'partial') items.push({ k: 'warn', t: 'A second screen only partly agrees', d: 'The other screen of ' + dn + ' (ID ' + Tw.drug.id + ') points the same way but gives a different answer, ' + tx + '.' });
      else if (a === 'disagree') items.push({ k: 'warn', t: 'A second screen disagrees', d: 'The other screen of ' + dn + ' (ID ' + Tw.drug.id + ') gives a different answer, ' + tx + '. Treat this result as unconfirmed.' });
      else items.push({ k: 'info', t: 'A second screen exists', d: 'The other screen of ' + dn + ' (ID ' + Tw.drug.id + ') has too few lines here to check.' });
    });
    if (R.co) {
      var co = R.co, left1 = R.n1 - co.a, left0 = R.n0 - co.c, canLeave = left1 >= 3 && left0 >= 3;
      var coTxt = cap(co.f.label) + ' is found in ' + co.a + ' of ' + R.n1 + ' ' + f.altName + ' and ' + co.c + ' of ' + R.n0 + ' ' + f.refName + '. It could explain part of the gap.';
      if (!canLeave) coTxt += ' Leaving those lines out would leave only ' + (left0 < 3 ? nLines(left0, f, 0) : nLines(left1, f, 1)) + ', too few to compare.';
      items.push({
        k: 'warn', t: 'Another gene change differs between the groups', d: coTxt,
        btn: canLeave ? '<button type="button" class="act-btn" data-leave="' + esc(co.f.key) + '">Leave out ' + esc(co.f.altName) + '</button>' : ''
      });
    }
    if (R.F) items.push({ k: 'info', t: cap(R.F.f.altName) + ' are left out', d: 'Both groups now exclude them, so they cannot explain the gap.', btn: '<button type="button" class="act-btn" data-undo="1">Show all lines again</button>' });
    if (R.ext1 / R.n1 >= 0.5 && R.ext0 / R.n0 >= 0.5) items.push({ k: 'warn', t: 'Many values are estimates', d: 'Most lines in both groups never dropped to half their growth, even at the highest dose tested.' });
    var ag = aucAgreement(R);
    if (ag === 'disagree' || ag === 'partial') items.push({ k: 'warn', t: 'A second score does not back this up', d: 'The area under the dose curve (AUC), a score that needs no estimates, ' + (ag === 'disagree' ? 'points the other way.' : 'does not pass the same cutoff.') });
    $('checks').innerHTML = items.map(function (it) {
      var glyph = it.k === 'no' ? 'no' : it.k;
      return '<li class="check">' + iconDot(it.k, glyph) + '<div><p class="check-t">' + nw(it.t) + '</p><p class="check-d">' + nw(it.d) + '</p>' + (it.btn || '') + '</div></li>';
    }).join('');
  }

  // ------------------------------------------------------------------ other cancer types
  var TYPESC = { pts: [], W: 700, H: 300 };
  function gapLabel(r) {
    if (!isFinite(r.fold)) return { txt: 'too few lines', cls: 't' };
    var v = verdictOf(r.p, r.fold).code, x = r.fold >= 1 ? r.fold : 1 / r.fold, dir = r.fold >= 1 ? 'less' : 'more';
    if (v === 'sens' || v === 'res') return { txt: fmtFold(x) + '× ' + dir, cls: 't t-l t-ink t-b' };
    if (v === 'small') return { txt: fmtFold(x) + '× ' + dir + ', small', cls: 't t-2' };
    return { txt: 'no clear gap', cls: 't' };
  }
  function renderTypes() {
    var R = ST.last, d = R.drug, f = R.feat, rows = R.across, fig = $('typesFig'), svg = $('types'), dn = drugText(d);
    $('typesTip').hidden = true;
    if (!rows.length) {
      $('typesSub').textContent = 'No cell line tested with ' + dn + ' ' + f.has + ', so there is nothing to compare across cancer types.';
      fig.hidden = true;
      return;
    }
    fig.hidden = false;
    var main = rows.filter(function (r) { return r.alt.length >= 3; });
    var selT = R.S.pooled ? -1 : R.S.type;
    var nDef = Math.min(rows.length, Math.max(4, Math.min(main.length, 8)));
    var shown = ST.typesAll ? rows.slice() : rows.slice(0, nDef);
    if (!ST.typesAll && selT >= 0 && shown.every(function (r) { return r.t !== selT; })) {
      var selRow = rows.filter(function (r) { return r.t === selT; })[0];
      if (selRow) shown.push(selRow);
    }
    // the sentence above the chart
    var sub = '', best = main.slice().sort(function (a, b) { return a.med1 - b.med1; })[0];
    if (!R.S.pooled) {
      var sel = rows.filter(function (r) { return r.t === selT; })[0];
      if (sel && best && sel.alt.length >= 3 && best.t !== sel.t) {
        var ratio = Math.exp(sel.med1 - best.med1);
        sub = cap(groupNounT(f, sel.t)) + ' needed a median of ' + fmtConc(Math.exp(sel.med1)) + ' µM of ' + dn + '. ' + cap(groupNounT(f, best.t)) + ' needed ' + fmtConc(Math.exp(best.med1)) + ' µM, ' + fmtFold(ratio) + ' times less.';
        if (ratio >= 3) sub += ' The same gene change does not guarantee the same response in every tissue.';
      } else if (sel && best && best.t === sel.t && main.length > 1) {
        sub = cap(groupNounT(f, sel.t)) + ' needed less ' + dn + ' than lines with the same gene change from any other cancer type.';
      }
    }
    if (!sub && main.length > 1) {
      var srt = main.slice().sort(function (a, b) { return a.med1 - b.med1; }), b0 = srt[0], b1 = srt[srt.length - 1];
      sub = 'Among cancer types with at least 3 ' + f.altName + ', the median dose runs from ' + fmtConc(Math.exp(b0.med1)) + ' µM in ' + typeNoun(TYPES[b0.t]) + ' to ' + fmtConc(Math.exp(b1.med1)) + ' µM in ' + typeNoun(TYPES[b1.t]) + '.';
    }
    if (!sub) sub = 'Only one cancer type has 3 or more ' + f.altName + ' tested with ' + dn + ', so there is little to compare across tissues.';
    $('typesSub').innerHTML = nw(sub);
    $('typesLegend').innerHTML = '<li><span class="sw sw-gene"></span>Median, ' + esc(f.altShort) + '</li><li><span class="sw sw-no"></span>Median, ' + esc(f.refShort) + '</li><li><span class="sw sw-zone"></span>Past the highest dose tested</li>';

    var W = boxWidth('typesBox'), stacked = W < 700;
    var LW = stacked ? 0 : 236, RW = stacked ? 0 : 128, RH = stacked ? 62 : 48, TOP = 8, AXH = 50;
    var H = TOP + shown.length * RH + AXH;
    var x0 = stacked ? 12 : LW + 8, x1 = stacked ? W - 12 : W - RW - 8;
    var vals = [], lmax = d.lnMax / LN10;
    shown.forEach(function (r) { vals.push(r.med1 / LN10); if (isFinite(r.med0)) vals.push(r.med0 / LN10); });
    vals.push(lmax);
    var dom = niceDomain(Math.min.apply(null, vals), Math.max.apply(null, vals), 0.3), lo = dom[0], hi = dom[1];
    function X(l) { return x0 + (l - lo) / (hi - lo) * (x1 - x0); }
    TYPESC.W = W; TYPESC.H = H; TYPESC.pts = [];
    svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    var s = '', yBot = TOP + shown.length * RH, xm = X(lmax);
    s += RECT(xm, TOP, x1 - xm, yBot - TOP, 'zone');
    shown.forEach(function (r, i) { if (r.t === selT) s += RECT(0, TOP + i * RH + 2, W, RH - 4, 'row-hl', 8); });
    logTicks(lo, hi).forEach(function (t) {
      var x = X(t);
      s += L(x, TOP, x, yBot, 'grid');
      s += T(x, yBot + 18, fmtTick(t), 'middle', 't');
    });
    s += L(xm, TOP, xm, yBot, 'edge');
    shown.forEach(function (r, i) {
      var top = TOP + i * RH, cy = stacked ? top + 42 : top + RH / 2, g = gapLabel(r);
      var countTxt = r.alt.length + ' with, ' + r.ref.length + ' without';
      if (stacked) {
        s += T(x0, top + 18, TYPES[r.t], 'start', 't t-l t-ink t-b');
        s += T(x1, top + 18, g.txt, 'end', g.cls);
      } else {
        s += T(8, cy - 3, TYPES[r.t], 'start', 't t-l t-ink t-b');
        s += T(8, cy + 14, countTxt, 'start', 't');
        s += T(W - 8, cy + 5, g.txt, 'end', g.cls);
      }
      if (isFinite(r.med0)) s += L(X(r.med0 / LN10), cy, X(r.med1 / LN10), cy, 'link');
      if (isFinite(r.med0)) {
        var xr = X(r.med0 / LN10);
        s += '<circle class="dot dot-no" cx="' + xr.toFixed(1) + '" cy="' + cy.toFixed(1) + '" r="6.5"/>';
        TYPESC.pts.push({ x: xr, y: cy, r: r });
      }
      var xa = X(r.med1 / LN10);
      s += '<circle class="dot dot-gene" cx="' + xa.toFixed(1) + '" cy="' + cy.toFixed(1) + '" r="6.5"/>';
      TYPESC.pts.push({ x: xa, y: cy, r: r });
    });
    s += L(x0, yBot, x1, yBot, 'axis');
    s += T((x0 + x1) / 2, yBot + 40, 'Median IC50 (µM, log scale)', 'middle', 't t-2');
    svg.innerHTML = s;
    var tg = $('typesToggle');
    if (rows.length > nDef) {
      tg.hidden = false;
      tg.textContent = ST.typesAll ? 'Show fewer cancer types' : 'Show all ' + rows.length + ' cancer types';
    } else tg.hidden = true;
  }
  function typesHover(e) {
    var svg = $('types'), tip = $('typesTip'), p = svgPoint(svg, TYPESC.W, TYPESC.H, e), best = nearest(TYPESC.pts, p, 20);
    if (!best) { tip.hidden = true; return; }
    var r = best.r, f = ST.last.feat;
    var h = '<b>' + esc(TYPES[r.t]) + '</b><br>' + esc(f.altShort) + ', median ' + fmtConc(Math.exp(r.med1)) + ' µM (' + r.alt.length + ' ' + plural(r.alt.length, 'line') + ')' +
      (isFinite(r.med0) ? '<br>' + esc(f.refShort) + ', median ' + fmtConc(Math.exp(r.med0)) + ' µM (' + r.ref.length + ' ' + plural(r.ref.length, 'line') + ')' : '');
    tipShow(tip, $('typesBox'), TYPESC.W, best.x, best.y, h);
  }

  // ------------------------------------------------------------------ more detail
  function renderDetails(R) {
    var f = R.feat, d = R.drug, h = '';
    if (R.status === 'ok' || R.status === 'few' || R.status === 'nomarker') {
      var m1 = R.n1 ? Math.exp(GS.median(R.g1.map(function (r) { return r.v; }))) : NaN;
      var m0 = R.n0 ? Math.exp(GS.median(R.g0.map(function (r) { return r.v; }))) : NaN;
      var a1 = R.n1 ? GS.median(R.g1.map(function (r) { return r.a; }).filter(isFinite)) : NaN;
      var a0 = R.n0 ? GS.median(R.g0.map(function (r) { return r.a; }).filter(isFinite)) : NaN;
      h += '<thead><tr><th></th><th class="num">' + esc(f.altShort) + '</th><th class="num">' + esc(f.refShort) + '</th></tr></thead><tbody>';
      h += '<tr><td>Cell lines</td><td class="num">' + R.n1 + '</td><td class="num">' + R.n0 + '</td></tr>';
      h += '<tr><td>Median IC50 (µM)</td><td class="num">' + fmtConc(m1) + '</td><td class="num">' + fmtConc(m0) + '</td></tr>';
      h += '<tr><td>Past the highest dose (' + fmtConc(d.maxc) + ' µM), estimated</td><td class="num">' + R.ext1 + ' of ' + R.n1 + '</td><td class="num">' + R.ext0 + ' of ' + R.n0 + '</td></tr>';
      h += '<tr><td>Median AUC (area under the dose curve)</td><td class="num">' + fmtAUC(a1) + '</td><td class="num">' + fmtAUC(a0) + '</td></tr>';
      if (R.status === 'ok') {
        var lo = Math.exp(R.hl.lo), hi = Math.exp(R.hl.hi), x = R.fold >= 1 ? R.fold : 1 / R.fold;
        var rng = R.fold >= 1 ? fmtFold(lo) + ' to ' + fmtFold(hi) : fmtFold(1 / hi) + ' to ' + fmtFold(1 / lo);
        h += '<tr><td>Typical fold difference, all pairs of lines (95% range)</td><td class="num" colspan="2">' + fmtFold(x) + ' times ' + (R.fold >= 1 ? 'less' : 'more') + ' drug (' + rng + ')</td></tr>';
        h += '<tr><td>Mann-Whitney U test on IC50</td><td class="num" colspan="2">' + pEq(R.mw.p) + (R.mw.method === 'exact' ? ', exact' : ', normal approximation') + '</td></tr>';
        if (R.mwA) h += '<tr><td>Mann-Whitney U test on AUC</td><td class="num" colspan="2">' + pEq(R.mwA.p) + '</td></tr>';
        if (R.S.pooled) {
          if (R.adj.ok) h += '<tr><td>Within cancer types (' + R.adj.informative + ' ' + plural(R.adj.informative, 'type') + ' with lines on both sides)</td><td class="num" colspan="2">' + fmtFold(foldX(R.main)) + ' times ' + lessMore(R.main) + ' drug, ' + pEq(R.adj.p) + '</td></tr>';
          else h += '<tr><td>Within cancer types</td><td class="num" colspan="2">Not possible here</td></tr>';
        }
      }
      h += '</tbody>';
    } else {
      h = '<tbody><tr><td>No numbers to show for this question.</td></tr></tbody>';
    }
    $('numTable').innerHTML = h;
    // by cancer type
    var t = '<thead><tr><th>Cancer type</th><th class="num">With</th><th class="num">Without</th><th class="num">Median with (µM)</th><th class="num">Median without (µM)</th><th class="num">Gap</th><th class="num">p</th></tr></thead><tbody>';
    if (!R.across.length) t += '<tr><td colspan="7">No tested line carries this gene change.</td></tr>';
    R.across.forEach(function (r) {
      var gap = isFinite(r.fold) ? fmtFold(r.fold >= 1 ? r.fold : 1 / r.fold) + '× ' + (r.fold >= 1 ? 'less' : 'more') : 'too few';
      t += '<tr><td>' + esc(TYPES[r.t]) + '</td><td class="num">' + r.alt.length + '</td><td class="num">' + r.ref.length + '</td><td class="num">' + fmtConc(Math.exp(r.med1)) + '</td><td class="num">' + (isFinite(r.med0) ? fmtConc(Math.exp(r.med0)) : 'n/a') + '</td><td class="num">' + gap + '</td><td class="num">' + (isFinite(r.p) ? fmtP(r.p) : '') + '</td></tr>';
    });
    $('typeTable').innerHTML = t + '</tbody>';
    // every cell line
    var rows = R.rows.slice().sort(function (a, b) { return a.v - b.v; });
    $('lineTitle').textContent = 'Every cell line (' + rows.length + ')';
    var lh = '<thead><tr><th>Cell line</th><th>Cancer type</th><th>Gene change</th><th class="num">IC50 (µM)</th><th>Measured or estimated</th></tr></thead><tbody>';
    if (!rows.length) lh += '<tr><td colspan="5">No cell lines for this question.</td></tr>';
    rows.forEach(function (r) {
      lh += '<tr><td>' + esc(CELLNAME[r.c]) + '</td><td>' + esc(TYPES[r.t]) + '</td><td>' + (r.alt ? 'Yes' : 'No') + '</td><td class="num">' + fmtConc(Math.exp(r.v)) + '</td><td>' + (r.v > d.lnMax ? 'Estimated' : 'Measured') + '</td></tr>';
    });
    $('lineTable').innerHTML = lh + '</tbody>';
    $('copyLines').hidden = !rows.length;
  }
  function linesCSV(R) {
    var f = R.feat, d = R.drug, out = ['cell_line,cosmic_id,cancer_type,' + f.key.replace(/,/g, ';') + ',drug,drug_id,ln_ic50,ic50_uM,auc,max_dose_uM,within_tested_doses'];
    R.rows.slice().sort(function (a, b) { return a.v - b.v; }).forEach(function (r) {
      out.push(['"' + CELLNAME[r.c] + '"', COSMIC[r.c], '"' + TYPES[r.t] + '"', r.alt, '"' + d.name + '"', d.id, r.v.toFixed(4), Math.exp(r.v).toPrecision(4), isFinite(r.a) ? r.a.toFixed(4) : '', d.maxc, r.v > d.lnMax ? 'no' : 'yes'].join(','));
    });
    return out.join('\n');
  }
  function copyText(text, toastEl) {
    function done(msg) { toastEl.textContent = msg; toastEl.hidden = false; setTimeout(function () { toastEl.hidden = true; }, 4000); }
    function fallback() {
      var old = toastEl.parentNode.parentNode.querySelector('.copy-box');
      if (old) old.remove();
      var ta = document.createElement('textarea');
      ta.className = 'copy-box';
      ta.value = text; ta.setAttribute('readonly', '');
      toastEl.parentNode.insertAdjacentElement('afterend', ta);
      ta.focus(); ta.select();
      done('Copying was blocked here. The table is selected below, so press Ctrl+C or Cmd+C.');
    }
    try {
      navigator.clipboard.writeText(text).then(function () { done('Copied ' + (text.split('\n').length - 1) + ' rows as comma-separated text'); }, fallback);
    } catch (e) { fallback(); }
  }

  // ------------------------------------------------------------------ examples
  var ERBB2 = 'gain_cnaPANCAN301_(CDK12,ERBB2,MED24)';
  var EXAMPLES = [
    { t: 'Known answer: BRAF in melanoma', drug: 1036, gene: 'BRAF_mut', set: 'Melanoma' },
    { t: 'Same mutation, colon cancer', drug: 1036, gene: 'BRAF_mut', set: 'Colorectal cancer' },
    { t: 'Control drug with no BRAF target', drug: 1005, gene: 'BRAF_mut', set: 'Melanoma' },
    { t: 'HER2 in breast cancer', drug: 119, gene: ERBB2, set: 'Breast cancer' },
    { t: 'EGFR in lung cancer', drug: 1010, gene: 'EGFR_mut', set: 'Lung adenocarcinoma' },
    { t: 'A gene that predicts resistance', drug: 1047, gene: 'TP53_mut', set: 'all' },
    { t: 'Too little data to answer', drug: 1, gene: 'EGFR_mut', set: 'Lung adenocarcinoma' }
  ];
  function renderExamples() {
    $('examples').innerHTML = EXAMPLES.map(function (ex, i) {
      return '<button type="button" class="ex-btn" data-ex="' + i + '" aria-pressed="false">' + esc(ex.t) + '</button>';
    }).join('');
    $('examples').addEventListener('click', function (e) {
      var b = e.target.closest('[data-ex]');
      if (!b) return;
      var ex = EXAMPLES[+b.getAttribute('data-ex')];
      setQuestion(ex.drug, ex.gene, setVal(ex.set));
    });
  }
  function markExamples() {
    document.querySelectorAll('.ex-btn').forEach(function (b) {
      var ex = EXAMPLES[+b.getAttribute('data-ex')];
      var on = !ST.filter && ex.drug === ST.drug && ex.gene === ST.gene && setVal(ex.set) === ST.set;
      b.setAttribute('aria-pressed', on ? 'true' : 'false');
    });
  }

  // ------------------------------------------------------------------ how it works
  var KA = [
    { drug: 1036, gene: 'BRAF_mut', set: 'Melanoma', expect: 'sens', q: 'PLX-4720 and a BRAF mutation in melanoma', why: 'BRAF inhibitors work in BRAF-mutant melanoma (Chapman 2011, Garnett 2012)' },
    { drug: 1371, gene: 'BRAF_mut', set: 'Melanoma', expect: 'sens', q: 'PLX-4720, second screen, and a BRAF mutation in melanoma', why: 'The same drug, screened a second time' },
    { drug: 1373, gene: 'BRAF_mut', set: 'Melanoma', expect: 'sens', q: 'Dabrafenib and a BRAF mutation in melanoma', why: 'Dabrafenib is an approved BRAF inhibitor for BRAF V600 melanoma' },
    { drug: 1010, gene: 'EGFR_mut', set: 'Lung adenocarcinoma', expect: 'sens', q: 'Gefitinib and an EGFR mutation in lung adenocarcinoma', why: 'EGFR-mutant lung cancers respond to gefitinib (Lynch 2004)' },
    { drug: 119, gene: ERBB2, set: 'Breast cancer', expect: 'sens', q: 'Lapatinib and ERBB2 (HER2) amplification in breast cancer', why: 'HER2-amplified breast cancers respond to lapatinib (Geyer 2006)' },
    { drug: 1013, gene: 'BCR-ABL_mut', set: 'all', expect: 'sens', q: 'Nilotinib and the BCR-ABL fusion, all cancer types', why: 'Nilotinib blocks the BCR-ABL kinase (Druker 2001 shows this for imatinib)' },
    { drug: 1047, gene: 'TP53_mut', set: 'all', expect: 'res', q: 'Nutlin-3a and a TP53 mutation, all cancer types', why: 'Nutlin-3a needs working p53 to kill cells (Vassilev 2004)' },
    { drug: 1005, gene: 'BRAF_mut', set: 'Melanoma', expect: 'none', q: 'Cisplatin and a BRAF mutation in melanoma', why: 'Cisplatin damages DNA and has no BRAF target, so no link is expected' }
  ];
  function renderKnownAnswers() {
    var pass = 0, h = '<thead><tr><th>Question</th><th>Expected</th><th>The tool found</th><th>Match</th></tr></thead><tbody>';
    KA.forEach(function (k, i) {
      var R = compare(DRUG_BY_ID[k.drug], FEAT_BY_KEY[k.gene], setInfo(setVal(k.set)));
      var got = R.status === 'ok' ? R.verdict.code : R.status, ok = got === k.expect;
      if (ok) pass++;
      var says = R.status === 'ok' ? (ok ? '' : LABEL[got] + ', ') + fmtFold(foldX(R.main)) + ' times ' + lessMore(R.main) + ' drug, ' + pEq(R.main.p) : 'No answer';
      h += '<tr data-ka="' + i + '" tabindex="0"><td>' + esc(k.q) + '<span class="why">' + esc(k.why) + '</span></td><td>' + LABEL[k.expect] + '</td><td>' + esc(says) + '</td><td>' +
        (ok ? '<span class="yes">' + iconDot('good') + 'Yes</span>' : '<span class="yes" style="color:var(--stop)">' + iconDot('stop') + 'No</span>') + '</td></tr>';
    });
    $('kaTable').innerHTML = h + '</tbody>';
    $('kaSummary').textContent = 'The tool gets ' + pass + ' of these ' + KA.length + ' known cases right. Seven are drug and gene pairs with a proven link, and one is a control drug where no link is expected.';
    function openKA(tr) {
      var k = KA[+tr.getAttribute('data-ka')];
      showView('ask');
      setQuestion(k.drug, k.gene, setVal(k.set));
      window.scrollTo({ top: 0, behavior: 'smooth' });
    }
    $('kaTable').addEventListener('click', function (e) { var tr = e.target.closest('[data-ka]'); if (tr) openKA(tr); });
    $('kaTable').addEventListener('keydown', function (e) { if (e.key === 'Enter') { var tr = e.target.closest('[data-ka]'); if (tr) openKA(tr); } });
  }
  var REFS = [
    ['Bauer DF. Constructing confidence sets using rank statistics. J Am Stat Assoc. 1972;67(339):687-690.', 'https://doi.org/10.1080/01621459.1972.10481279'],
    ['Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B. 1995;57(1):289-300.', 'https://doi.org/10.1111/j.2517-6161.1995.tb02031.x'],
    ['Chapman PB, Hauschild A, Robert C, et al. Improved survival with vemurafenib in melanoma with BRAF V600E mutation. N Engl J Med. 2011;364(26):2507-2516.', 'https://doi.org/10.1056/NEJMoa1103782'],
    ['Cokelaer T, Chen E, Iorio F, et al. GDSCTools for mining pharmacogenomic interactions in cancer. Bioinformatics. 2018;34(7):1226-1228.', 'https://doi.org/10.1093/bioinformatics/btx744'],
    ['Druker BJ, Talpaz M, Resta DJ, et al. Efficacy and safety of a specific inhibitor of the BCR-ABL tyrosine kinase in chronic myeloid leukemia. N Engl J Med. 2001;344(14):1031-1037.', 'https://doi.org/10.1056/NEJM200104053441401'],
    ['Garnett MJ, Edelman EJ, Heidorn SJ, et al. Systematic identification of genomic markers of drug sensitivity in cancer cells. Nature. 2012;483(7391):570-575.', 'https://doi.org/10.1038/nature11005'],
    ['Geyer CE, Forster J, Lindquist D, et al. Lapatinib plus capecitabine for HER2-positive advanced breast cancer. N Engl J Med. 2006;355(26):2733-2743.', 'https://doi.org/10.1056/NEJMoa064320'],
    ['Haibe-Kains B, El-Hachem N, Birkbak NJ, et al. Inconsistency in large pharmacogenomic studies. Nature. 2013;504(7480):389-393.', 'https://doi.org/10.1038/nature12831'],
    ['Hodges JL, Lehmann EL. Estimates of location based on rank tests. Ann Math Stat. 1963;34(2):598-611.', 'https://doi.org/10.1214/aoms/1177704172'],
    ['Iorio F, Knijnenburg TA, Vis DJ, et al. A landscape of pharmacogenomic interactions in cancer. Cell. 2016;166(3):740-754.', 'https://doi.org/10.1016/j.cell.2016.06.017'],
    ['Lynch TJ, Bell DW, Sordella R, et al. Activating mutations in the epidermal growth factor receptor underlying responsiveness of non-small-cell lung cancer to gefitinib. N Engl J Med. 2004;350(21):2129-2139.', 'https://doi.org/10.1056/NEJMoa040938'],
    ['Mann HB, Whitney DR. On a test of whether one of two random variables is stochastically larger than the other. Ann Math Stat. 1947;18(1):50-60.', 'https://doi.org/10.1214/aoms/1177730491'],
    ['Prahallad A, Sun C, Huang S, et al. Unresponsiveness of colon cancer to BRAF(V600E) inhibition through feedback activation of EGFR. Nature. 2012;483(7387):100-103.', 'https://doi.org/10.1038/nature10868'],
    ['Vassilev LT, Vu BT, Graves B, et al. In vivo activation of the p53 pathway by small-molecule antagonists of MDM2. Science. 2004;303(5659):844-848.', 'https://doi.org/10.1126/science.1092472'],
    ['Vis DJ, Bombardelli L, Lightfoot H, et al. Multilevel models improve precision and speed of IC50 estimates. Pharmacogenomics. 2016;17(7):691-700.', 'https://doi.org/10.2217/pgs.16.15'],
    ['Yang W, Soares J, Greninger P, et al. Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells. Nucleic Acids Res. 2013;41(D1):D955-D961.', 'https://doi.org/10.1093/nar/gks1111']
  ];
  function renderRefs() {
    $('refs').innerHTML = REFS.map(function (r) { return '<li>' + esc(r[0]) + ' <a href="' + r[1] + '" target="_blank" rel="noopener">' + esc(r[1].replace('https://', '')) + '</a></li>'; }).join('');
  }

  // ------------------------------------------------------------------ rank drugs
  // fast version of compare() for ranking: no confidence range, capped exact work
  function compareLite(drug, feat, S) {
    var off = drug.i * NC, x1 = [], x0 = [], y = [], a = [], st = [];
    for (var c = 0; c < NC; c++) {
      var v = IC[off + c];
      if (v !== v) continue;
      var t = CELLTYPE[c];
      if (!S.mask[t]) continue;
      var m = feat.mask[c];
      (m ? x1 : x0).push(v);
      if (S.pooled) { y.push(v); a.push(m); st.push(t); }
    }
    if (x1.length < 3 || x0.length < 3) return null;
    if (S.pooled) {
      var ad = GS.adjustedDiff(y, a, st);
      if (!ad.ok) return null;
      return { n1: x1.length, n0: x0.length, fold: Math.exp(-ad.beta), p: ad.p };
    }
    var mw = GS.mannWhitney(x1, x0, 3e5);
    return { n1: x1.length, n0: x0.length, fold: Math.exp(GS.hodgesLehmann(x1, x0, 0, true).est), p: mw.p };
  }
  var RK = { mode: 'drugs', gene: 'BRAF_mut', drug: 1036, set: typeVal('Melanoma'), res: null, allBetter: false, allWorse: false };
  function setRankMode(m) {
    RK.mode = m; RK.allBetter = RK.allWorse = false;
    $('rk-mode-drugs').setAttribute('aria-checked', m === 'drugs' ? 'true' : 'false');
    $('rk-mode-genes').setAttribute('aria-checked', m === 'genes' ? 'true' : 'false');
    $('pk-rkgene').hidden = m !== 'drugs';
    $('pk-rkdrug').hidden = m !== 'genes';
    closeAll(null);
    runRank();
  }
  function runRank() {
    var S = setInfo(RK.set), out = [];
    if (RK.mode === 'drugs') {
      var f = FEAT_BY_KEY[RK.gene];
      DRUGS.forEach(function (d) { var r = compareLite(d, f, S); if (r) { r.drug = d; r.feat = f; out.push(r); } });
    } else {
      var d = DRUG_BY_ID[RK.drug];
      FEATS.forEach(function (f2) { var r = compareLite(d, f2, S); if (r) { r.drug = d; r.feat = f2; out.push(r); } });
    }
    var q = out.length ? GS.bh(out.map(function (r) { return r.p; })) : [];
    out.forEach(function (r, i) {
      r.i = i; r.q = q[i];
      r.big = r.fold >= 2 || r.fold <= 0.5;
      r.clear = r.q < 0.05 && r.big;
      r.possible = !r.clear && r.p < 0.05 && r.big;
    });
    RK.res = { S: S, out: out };
    renderRank();
  }
  function rankRows(rows, all, better) {
    if (!rows.length) return '<li class="rank-empty">' + (better ? 'None in this group.' : 'None in this group.') + '</li>';
    var maxX = Math.max.apply(null, rows.map(function (r) { return foldX(r); })), lmax = Math.log10(Math.max(maxX, 2.5));
    return (all ? rows : rows.slice(0, 10)).map(function (r, k) {
      var x = foldX(r), w = Math.max(4, Math.round(Math.log10(x) / lmax * 100));
      var isD = RK.mode === 'drugs';
      var name = isD ? r.drug.name : r.feat.label;
      var meta = (isD ? 'Targets ' + shortText(r.drug.target, 34) + (r.drug.twin.length ? ' (screen ' + r.drug.id + ')' : '') + ', ' : '') + r.n1 + ' lines with and ' + r.n0 + ' without';
      var tag = r.clear ? '<span class="rank-tag k-good">' + svgIcon('good') + 'Clear link</span>' : '<span class="rank-tag k-warn">' + svgIcon('warn') + 'Possible link</span>';
      return '<li><button type="button" class="rank-row" data-open="' + r.i + '"><span class="rank-n">' + (k + 1) + '</span>' +
        '<span class="rank-name"><b>' + nw(name) + '</b><span class="rank-meta">' + nw(meta) + '</span></span>' +
        '<span class="rank-fold"><span>' + fmtFold(x) + ' times ' + (r.fold >= 1 ? 'less' : 'more') + ' drug</span><span class="rank-bar" aria-hidden="true"><i style="width:' + w + '%"></i></span></span>' + tag + '</button></li>';
    }).join('');
  }
  function renderRank() {
    var res = RK.res, S = res.S, out = res.out, isD = RK.mode === 'drugs';
    var f = FEAT_BY_KEY[RK.gene], d = DRUG_BY_ID[RK.drug], dn = drugText(d);
    function order(a, b) { return (b.clear - a.clear) || (foldX(b) - foldX(a)); }
    var better = out.filter(function (r) { return r.p < 0.05 && r.fold >= 2; }).sort(order);
    var worse = out.filter(function (r) { return r.p < 0.05 && r.fold <= 0.5; }).sort(order);
    // count drugs by name, since 15 drugs were screened twice
    function nUnique(rows) { var seen = {}; rows.forEach(function (r) { seen[isD ? r.drug.name : r.feat.key] = 1; }); return Object.keys(seen).length; }
    var nClear = nUnique(out.filter(function (r) { return r.clear; })), nPoss = nUnique(out.filter(function (r) { return r.possible; })), nAll = nUnique(out), nDup = out.length - nAll;
    $('rkAsked').innerHTML = '<b>Your question.</b> ' + nw(isD
      ? (S.pooled ? cap(S.where) + ', which drugs work differently on ' + f.altName + '?' : 'Which drugs work differently on ' + groupNoun(f, S) + ' than on ' + otherNoun(S) + '?')
      : (S.pooled ? cap(S.where) + ', which gene changes make lines respond differently to ' + dn + '?' : 'Which gene changes make ' + S.noun + ' respond differently to ' + dn + '?'));
    var head, basis = [];
    if (!out.length) {
      head = isD ? 'No drug can be tested for ' + f.label + ' ' + S.where + '.' : 'No gene change can be tested with ' + dn + ' ' + S.where + '.';
      basis.push(isD ? 'Every drug needs at least 3 lines with the gene change and 3 without it. Try another cancer type, or All cancer types.' : 'Each gene change needs at least 3 tested lines with it and 3 without it. Try another cancer type, or All cancer types.');
    } else {
      head = nClear + ' ' + (isD ? plural(nClear, 'drug') : 'gene ' + plural(nClear, 'change')) + ' ' + (nClear === 1 ? 'shows' : 'show') + ' a clear link' + (isD ? ' with ' + f.label : ' with ' + dn) + ' ' + S.where + '.';
      basis.push('The tool tested every ' + (isD ? 'drug' : 'gene change') + ' with at least 3 lines on each side, ' + nAll + ' in all' + (isD && nDup ? ' (' + nDup + ' of them screened twice)' : '') + (S.pooled ? ', comparing lines only within the same cancer type' : '') + '. It used a stricter cutoff because it ran so many tests at once, so a clear link is unlikely to be luck.' +
        (nPoss ? ' ' + nPoss + ' more ' + (nPoss === 1 ? 'shows' : 'show') + ' a possible link that could be luck.' : ''));
    }
    $('rkHead').innerHTML = nw(head);
    $('rkBasis').innerHTML = basis.map(function (t) { return '<p>' + nw(t) + '</p>'; }).join('');
    $('rkBetterTitle').innerHTML = nw(isD ? 'Drugs that work better on ' + f.altName : 'Gene changes that make lines more sensitive to ' + dn);
    $('rkWorseTitle').innerHTML = nw(isD ? 'Drugs that work worse on ' + f.altName : 'Gene changes that make lines less sensitive to ' + dn);
    $('rkBetter').innerHTML = rankRows(better, RK.allBetter, true);
    $('rkWorse').innerHTML = rankRows(worse, RK.allWorse, false);
    var mb = $('rkBetterMore'), mw = $('rkWorseMore');
    mb.hidden = better.length <= 10; mb.textContent = RK.allBetter ? 'Show the top 10 only' : 'Show all ' + better.length;
    mw.hidden = worse.length <= 10; mw.textContent = RK.allWorse ? 'Show the top 10 only' : 'Show all ' + worse.length;
  }

  // ------------------------------------------------------------------ views
  function showView(name) {
    if (name !== 'how' && name !== 'rank') name = 'ask';
    ['ask', 'rank', 'how'].forEach(function (v) {
      $('view-' + v).hidden = v !== name;
      $('tab-' + v).setAttribute('aria-current', v === name ? 'page' : 'false');
    });
    closeAll(null);
    if (name === 'ask' && ST.last) { renderDose(); renderTypes(); }
    if (name === 'rank' && !RK.res) runRank();
    try { history.replaceState(null, '', '#' + name); } catch (e) { /* some viewers block this */ }
  }

  // ------------------------------------------------------------------ init
  function init() {
    Picker('drug', DRUG_ITEMS, function (v) { ST.drug = v; ST.filter = null; ST.typesAll = false; run(); }, drugStatuses);
    Picker('gene', GENE_ITEMS, function (v) { ST.gene = v; ST.filter = null; ST.typesAll = false; run(); }, geneStatuses);
    Picker('type', TYPE_ITEMS, function (v) { ST.set = v; ST.filter = null; ST.typesAll = false; run(); }, typeStatuses);
    Picker('rkgene', GENE_ITEMS, function (v) { RK.gene = v; RK.allBetter = RK.allWorse = false; runRank(); }, rankGeneStatuses, RANK_LABELS);
    Picker('rkdrug', DRUG_ITEMS, function (v) { RK.drug = v; RK.allBetter = RK.allWorse = false; runRank(); }, rankDrugStatuses, RANK_LABELS);
    Picker('rktype', TYPE_ITEMS, function (v) { RK.set = v; RK.allBetter = RK.allWorse = false; runRank(); }, rankTypeStatuses, RANK_LABELS);
    PICKERS.rkgene.set(RK.gene); PICKERS.rkdrug.set(RK.drug); PICKERS.rktype.set(RK.set);
    $('rankForm').addEventListener('submit', function (e) { e.preventDefault(); });
    ['drugs', 'genes'].forEach(function (m) {
      $('rk-mode-' + m).addEventListener('click', function () { setRankMode(m); });
    });
    $('rkBetterMore').addEventListener('click', function () { RK.allBetter = !RK.allBetter; renderRank(); });
    $('rkWorseMore').addEventListener('click', function () { RK.allWorse = !RK.allWorse; renderRank(); });
    function openRankRow(e) {
      var b = e.target.closest('[data-open]');
      if (!b) return;
      var r = RK.res.out[+b.getAttribute('data-open')];
      showView('ask');
      setQuestion(r.drug.id, r.feat.key, RK.set);
      window.scrollTo({ top: 0, behavior: 'smooth' });
    }
    $('rkBetter').addEventListener('click', openRankRow);
    $('rkWorse').addEventListener('click', openRankRow);
    PICKERS.drug.set(ST.drug); PICKERS.gene.set(ST.gene); PICKERS.type.set(ST.set);
    $('question').addEventListener('submit', function (e) { e.preventDefault(); });
    renderExamples();
    function onAction(e) {
      var lv = e.target.closest('[data-leave]'), un = e.target.closest('[data-undo]');
      if (lv) { ST.filter = { f: FEAT_BY_KEY[lv.getAttribute('data-leave')], keep: 0 }; run(); $('answer').scrollIntoView({ behavior: 'smooth', block: 'start' }); }
      else if (un) { ST.filter = null; run(); }
    }
    $('checks').addEventListener('click', onAction);
    $('filterNote').addEventListener('click', onAction);
    ['pointermove', 'pointerdown'].forEach(function (ev) {
      $('dose').addEventListener(ev, doseHover);
      $('types').addEventListener(ev, typesHover);
    });
    $('dose').addEventListener('pointerleave', function () { $('doseTip').hidden = true; var p = $('dose').querySelector('.dot.hl'); if (p) p.classList.remove('hl'); });
    $('types').addEventListener('pointerleave', function () { $('typesTip').hidden = true; });
    $('typesToggle').addEventListener('click', function () { ST.typesAll = !ST.typesAll; renderTypes(); });
    $('copyLines').addEventListener('click', function () { if (ST.last) copyText(linesCSV(ST.last), $('copyToast')); });
    document.querySelectorAll('.view-tab').forEach(function (b) {
      b.addEventListener('click', function () { showView(b.getAttribute('data-view')); });
    });
    renderKnownAnswers();
    renderRefs();
    run();
    var h = (location.hash || '').replace('#', '');
    showView(h === 'how' || h === 'about' || h === 'methods' ? 'how' : (h === 'rank' || h === 'scan') ? 'rank' : 'ask');
    var rt;
    window.addEventListener('resize', function () {
      clearTimeout(rt);
      rt = setTimeout(function () { if (!$('view-ask').hidden && ST.last) { renderDose(); renderTypes(); } }, 150);
    });
  }
  init();
})();
