/* Statistics for the drug-marker explorer. No dependencies. Exposes window.GS (or module.exports in Node). */
(function (root) {
  'use strict';

  function sorted(a) { return Float64Array.from(a).sort(); }

  function median(a) {
    var s = sorted(a), n = s.length;
    if (!n) return NaN;
    return n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2;
  }

  // Quantile of a sorted array, linear interpolation (same as R type 7 / numpy default)
  function quantileSorted(s, q) {
    var n = s.length;
    if (!n) return NaN;
    var h = (n - 1) * q, lo = Math.floor(h), hi = Math.ceil(h);
    return s[lo] + (h - lo) * (s[hi] - s[lo]);
  }

  // log-gamma, Lanczos approximation (g = 7, n = 9)
  var LG = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61503916999185, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  function lgamma(x) {
    if (x < 0.5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * x))) - lgamma(1 - x);
    x -= 1;
    var a = LG[0], t = x + 7.5;
    for (var i = 1; i < 9; i++) a += LG[i] / (x + i);
    return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
  }

  // Continued fraction for the incomplete beta function (modified Lentz)
  function betacf(a, b, x) {
    var MAXIT = 300, EPS = 3e-16, FPMIN = 1e-300;
    var qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    d = 1 / d;
    var h = d;
    for (var m = 1; m <= MAXIT; m++) {
      var m2 = 2 * m, aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d;
      var del = d * c;
      h *= del;
      if (Math.abs(del - 1) < EPS) break;
    }
    return h;
  }

  // Regularized incomplete beta I_x(a, b)
  function ibeta(x, a, b) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    var lbt = lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1 - x);
    if (x < (a + 1) / (a + b + 2)) return Math.exp(lbt) * betacf(a, b, x) / a;
    return 1 - Math.exp(lbt) * betacf(b, a, 1 - x) / b;
  }

  // Two-sided p-value for Student's t
  function tTwoSidedP(t, df) {
    if (!isFinite(t)) return 0;
    return ibeta(df / (df + t * t), df / 2, 0.5);
  }

  // Complementary error function with small relative error in the tails (Numerical Recipes erfcc)
  function erfc(x) {
    var z = Math.abs(x), t = 1 / (1 + 0.5 * z);
    var r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
      t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
      t * (-0.82215223 + t * 0.17087277)))))))));
    return x >= 0 ? r : 2 - r;
  }
  function normTwoSidedP(z) { return erfc(Math.abs(z) / Math.SQRT2); }

  // Exact cumulative distribution P(U <= v), v = 0..uMax, of the Mann-Whitney U statistic for group sizes m and n
  // (no ties). Counts partitions of v into at most a parts, each part at most b (a = smaller group, b = larger).
  // Uses additions only, so it stays accurate far into the tail. Returns null when the work would exceed maxWork.
  function exactCdfArray(m, n, uMax, maxWork) {
    var a = Math.min(m, n), b = Math.max(m, n);
    uMax = Math.max(0, Math.min(Math.floor(uMax + 1e-9), a * b));
    if (b * a * (uMax + 1) > maxWork) return null;
    var W = uMax + 1, prev = new Float64Array((a + 1) * W), cur = new Float64Array((a + 1) * W), tmp, k, v;
    for (k = 0; k <= a; k++) prev[k * W] = 1;            // largest part 0: only the empty sum
    for (var j = 1; j <= b; j++) {
      cur.fill(0);
      cur[0] = 1;                                        // zero parts: only v = 0
      for (k = 1; k <= a; k++) {
        var base = k * W, baseK1 = (k - 1) * W;
        for (v = 0; v <= uMax; v++) {
          var val = cur[baseK1 + v];
          if (v >= k) val += prev[base + v - k];
          cur[base + v] = val;
        }
      }
      tmp = prev; prev = cur; cur = tmp;
    }
    var logC = lgamma(a + b + 1) - lgamma(a + 1) - lgamma(b + 1), out = new Float64Array(W), sum = 0, off = a * W;
    for (v = 0; v <= uMax; v++) { sum += prev[off + v]; out[v] = Math.min(1, Math.exp(Math.log(sum) - logC)); }
    return out;
  }

  // Exact lower tail P(U <= u). Returns null when too slow to compute.
  function exactLowerTail(m, n, u, maxWork) {
    u = Math.floor(u + 1e-9);
    if (u < 0) return 0;
    if (u >= m * n) return 1;
    var cdf = exactCdfArray(m, n, u, maxWork);
    return cdf ? cdf[u] : null;
  }

  // Hodges-Lehmann estimate of the shift x0 - x1 (median of all pairwise differences) with the
  // distribution-free 95% confidence interval that inverts the Mann-Whitney test (same rule as R's wilcox.test).
  function hodgesLehmann(x1, x0, maxWork, noCI) {
    var m = x1.length, n = x0.length, N = m * n, d = new Float64Array(N), k = 0, i, j;
    for (i = 0; i < m; i++) for (j = 0; j < n; j++) d[k++] = x0[j] - x1[i];
    d.sort();
    var est = N % 2 ? d[(N - 1) / 2] : (d[N / 2 - 1] + d[N / 2]) / 2;
    if (noCI) return { est: est };
    var mu = N / 2, sd = Math.sqrt(N * (m + n + 1) / 12), qu = null, method = 'exact';
    var uMax = Math.min(Math.floor(mu), Math.ceil(mu - 1.2 * sd) + 2);
    var cdf = exactCdfArray(m, n, uMax, maxWork || 4e7);
    if (cdf) { for (var u = 0; u < cdf.length; u++) if (cdf[u] >= 0.025) { qu = u; break; } }
    if (qu === null) { method = 'normal'; qu = Math.round(mu - 1.959964 * sd); }
    if (qu < 1) qu = 1;
    return { est: est, lo: d[qu - 1], hi: d[N - qu], method: method };
  }

  // Mann-Whitney U test, two-sided.
  // U1 counts pairs where the group-1 value is larger (ties count one half).
  // cles = probability that a random group-1 value is smaller than a random group-0 value.
  function mannWhitney(x1, x0, maxWork) {
    var n1 = x1.length, n0 = x0.length, N = n1 + n0, i, j, k;
    var all = new Array(N);
    for (i = 0; i < n1; i++) all[i] = [x1[i], 1];
    for (i = 0; i < n0; i++) all[n1 + i] = [x0[i], 0];
    all.sort(function (p, q) { return p[0] - q[0]; });
    var R1 = 0, tieTerm = 0, ties = false;
    i = 0;
    while (i < N) {
      j = i;
      while (j + 1 < N && all[j + 1][0] === all[i][0]) j++;
      var r = (i + j) / 2 + 1, t = j - i + 1;
      if (t > 1) { ties = true; tieTerm += t * t * t - t; }
      for (k = i; k <= j; k++) if (all[k][1] === 1) R1 += r;
      i = j + 1;
    }
    var U1 = R1 - n1 * (n1 + 1) / 2, mn = n1 * n0, mu = mn / 2;
    var res = { n1: n1, n0: n0, U1: U1, cles: 1 - U1 / mn, p: NaN, method: '' };
    if (!ties) {
      var tail = exactLowerTail(n1, n0, Math.min(U1, mn - U1), maxWork || 4e7);
      if (tail !== null) { res.p = Math.min(1, 2 * tail); res.method = 'exact'; return res; }
    }
    var sigma = Math.sqrt(mn / 12 * ((N + 1) - tieTerm / (N * (N - 1))));
    var z = sigma > 0 ? Math.max(0, Math.abs(U1 - mu) - 0.5) / sigma : 0;
    res.p = sigma > 0 ? Math.min(1, normTwoSidedP(z)) : 1;
    res.method = 'normal';
    return res;
  }

  function mulberry32(a) {
    return function () {
      a |= 0; a = a + 0x6D2B79F5 | 0;
      var t = Math.imul(a ^ a >>> 15, 1 | a);
      t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
      return ((t ^ t >>> 14) >>> 0) / 4294967296;
    };
  }
  function hashString(s) {
    var h = 2166136261;
    for (var i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 16777619); }
    return h >>> 0;
  }

  // Percentile bootstrap interval for median(x0) - median(x1). Resamples each group separately.
  function bootstrapMedianDiff(x1, x0, B, seed) {
    var rnd = mulberry32(seed), n1 = x1.length, n0 = x0.length;
    var s1 = new Float64Array(n1), s0 = new Float64Array(n0), d = new Float64Array(B), i, b;
    function med(s) { s.sort(); var n = s.length; return n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2; }
    for (b = 0; b < B; b++) {
      for (i = 0; i < n1; i++) s1[i] = x1[Math.floor(rnd() * n1)];
      for (i = 0; i < n0; i++) s0[i] = x0[Math.floor(rnd() * n0)];
      d[b] = med(s0) - med(s1);
    }
    d.sort();
    return [quantileSorted(d, 0.025), quantileSorted(d, 0.975)];
  }

  // Difference between groups adjusted for cancer type: y = type mean + beta * marker + error.
  // Fitted by removing each type's mean from y and from the marker (fixed-effects least squares).
  // beta is the average within-type difference (marker minus no marker) in y.
  function adjustedDiff(y, alt, strata) {
    var n = y.length, sums = {}, i, s;
    for (i = 0; i < n; i++) {
      s = sums[strata[i]] || (sums[strata[i]] = { n: 0, y: 0, x: 0 });
      s.n++; s.y += y[i]; s.x += alt[i];
    }
    var K = 0, informative = 0;
    for (var key in sums) {
      K++;
      var g = sums[key];
      if (g.x > 0 && g.x < g.n) informative++;
    }
    var yt = new Float64Array(n), xt = new Float64Array(n), sxx = 0, sxy = 0;
    for (i = 0; i < n; i++) {
      s = sums[strata[i]];
      yt[i] = y[i] - s.y / s.n; xt[i] = alt[i] - s.x / s.n;
      sxx += xt[i] * xt[i]; sxy += xt[i] * yt[i];
    }
    var df = n - K - 1;
    if (sxx < 1e-9 || df < 1) return { ok: false, K: K, informative: informative };
    var beta = sxy / sxx, sse = 0;
    for (i = 0; i < n; i++) { var e = yt[i] - beta * xt[i]; sse += e * e; }
    var se = Math.sqrt(sse / df / sxx), t = beta / se;
    return { ok: true, beta: beta, se: se, t: t, df: df, p: tTwoSidedP(t, df), K: K, informative: informative };
  }

  function lchoose(n, k) { return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1); }

  // Fisher's exact test, two-sided, for the 2x2 table [[a, b], [c, d]] (same rule as SciPy)
  function fisherExact(a, b, c, d) {
    var r1 = a + b, n = a + b + c + d, c1 = a + c;
    function lp(x) { return lchoose(c1, x) + lchoose(n - c1, r1 - x) - lchoose(n, r1); }
    var lo = Math.max(0, r1 - (n - c1)), hi = Math.min(r1, c1), l0 = lp(a), p = 0;
    for (var x = lo; x <= hi; x++) { var l = lp(x); if (l <= l0 + 1e-7) p += Math.exp(l); }
    return Math.min(1, p);
  }

  // Benjamini-Hochberg adjusted p-values (false discovery rate)
  function bh(p) {
    var m = p.length, idx = p.map(function (_, i) { return i; }), q = new Array(m), run = 1;
    idx.sort(function (a, b) { return p[a] - p[b]; });
    for (var r = m - 1; r >= 0; r--) {
      run = Math.min(run, p[idx[r]] * m / (r + 1));
      q[idx[r]] = Math.min(1, run);
    }
    return q;
  }

  var GS = {
    median: median, sorted: sorted, quantileSorted: quantileSorted, lgamma: lgamma, ibeta: ibeta,
    tTwoSidedP: tTwoSidedP, normTwoSidedP: normTwoSidedP, exactLowerTail: exactLowerTail, exactCdfArray: exactCdfArray, hodgesLehmann: hodgesLehmann,
    mannWhitney: mannWhitney, bootstrapMedianDiff: bootstrapMedianDiff, adjustedDiff: adjustedDiff, fisherExact: fisherExact,
    bh: bh, mulberry32: mulberry32, hashString: hashString
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = GS;
  else root.GS = GS;
})(typeof window !== 'undefined' ? window : this);
