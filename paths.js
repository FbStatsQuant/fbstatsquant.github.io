// Hero: interactive fan chart of simulated geometric Brownian motion paths.
// Hover: shows quantiles at that time. Click: re-simulate. Static if reduced motion.
(function () {
  var hero = document.querySelector('.hero');
  var canvas = document.getElementById('paths');
  if (!hero || !canvas || !canvas.getContext) return;
  var ctx = canvas.getContext('2d');
  var reduce = window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches;

  var N = 160, STEPS = 240;
  var w, h, x0, paths, q, progress = 1, mouse = null, raf = 0, animating = false;

  function randn() {
    var u = 0, v = 0;
    while (u === 0) u = Math.random();
    while (v === 0) v = Math.random();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  function simulate() {
    var mu = 0.06, sigma = 0.38, dt = 1 / STEPS, i, k;
    paths = [];
    for (i = 0; i < N; i++) {
      var x = 0, p = [0];
      for (k = 1; k <= STEPS; k++) {
        x += (mu - 0.5 * sigma * sigma) * dt + sigma * Math.sqrt(dt) * randn();
        p.push(x); // log price
      }
      paths.push(p);
    }
    q = [];
    for (k = 0; k <= STEPS; k++) {
      var col = paths.map(function (p) { return p[k]; }).sort(function (a, b) { return a - b; });
      q.push({
        p05: col[Math.floor(0.05 * N)], p25: col[Math.floor(0.25 * N)],
        p50: col[Math.floor(0.50 * N)], p75: col[Math.floor(0.75 * N)],
        p95: col[Math.floor(0.95 * N)]
      });
    }
  }

  function resize() {
    var r = hero.getBoundingClientRect();
    var dpr = window.devicePixelRatio || 1;
    w = r.width; h = r.height;
    canvas.width = w * dpr; canvas.height = h * dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    x0 = w > 820 ? w * 0.60 : 0;
  }

  function css(name, fallback) {
    var v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
    return v || fallback;
  }

  var RANGE = 1.15; // log-price half-range shown vertically
  function X(k) { return x0 + (k / STEPS) * (w - x0); }
  function Y(v) { return h * 0.5 - (Math.max(-RANGE, Math.min(RANGE, v)) / RANGE) * h * 0.42; }

  function band(lo, hi, upto, fill) {
    ctx.beginPath();
    var k;
    for (k = 0; k <= upto; k++) ctx[k ? 'lineTo' : 'moveTo'](X(k), Y(q[k][hi]));
    for (k = upto; k >= 0; k--) ctx.lineTo(X(k), Y(q[k][lo]));
    ctx.closePath();
    ctx.fillStyle = fill;
    ctx.fill();
  }

  function draw() {
    var a = css('--path', '31, 95, 214'), b = css('--path2', '124, 58, 237');
    var upto = Math.max(2, Math.floor(progress * STEPS));
    var i, k;
    ctx.clearRect(0, 0, w, h);

    var g = ctx.createLinearGradient(x0, 0, w, 0);
    g.addColorStop(0, 'rgb(' + a + ')');
    g.addColorStop(1, 'rgb(' + b + ')');

    // uncertainty bands
    band('p05', 'p95', upto, 'rgba(' + a + ',0.07)');
    band('p25', 'p75', upto, 'rgba(' + a + ',0.11)');

    // paths
    var hi = -1, best = 1e9;
    if (mouse && mouse.x > x0) {
      var km = Math.min(STEPS, Math.max(0, Math.round(((mouse.x - x0) / (w - x0)) * STEPS)));
      if (km <= upto) {
        for (i = 0; i < N; i++) {
          var d = Math.abs(Y(paths[i][km]) - mouse.y);
          if (d < best) { best = d; hi = i; }
        }
      }
    }
    ctx.lineWidth = 1;
    ctx.strokeStyle = g;
    ctx.globalAlpha = 0.16;
    for (i = 0; i < N; i++) {
      if (i === hi) continue;
      ctx.beginPath();
      for (k = 0; k <= upto; k++) ctx[k ? 'lineTo' : 'moveTo'](X(k), Y(paths[i][k]));
      ctx.stroke();
    }
    ctx.globalAlpha = 1;

    // median
    ctx.beginPath();
    for (k = 0; k <= upto; k++) ctx[k ? 'lineTo' : 'moveTo'](X(k), Y(q[k].p50));
    ctx.lineWidth = 2.2;
    ctx.strokeStyle = g;
    ctx.stroke();

    // highlighted path
    if (hi >= 0 && best < 80) {
      ctx.beginPath();
      for (k = 0; k <= upto; k++) ctx[k ? 'lineTo' : 'moveTo'](X(k), Y(paths[hi][k]));
      ctx.lineWidth = 1.8;
      ctx.strokeStyle = 'rgb(' + b + ')';
      ctx.stroke();
    }

    // hover readout
    if (mouse && mouse.x > x0 && progress >= 1) {
      var kk = Math.min(STEPS, Math.max(0, Math.round(((mouse.x - x0) / (w - x0)) * STEPS)));
      var cx = X(kk), c = q[kk];
      ctx.strokeStyle = 'rgba(' + a + ',0.35)';
      ctx.lineWidth = 1;
      ctx.setLineDash([4, 4]);
      ctx.beginPath(); ctx.moveTo(cx, 0); ctx.lineTo(cx, h); ctx.stroke();
      ctx.setLineDash([]);
      ['p05', 'p25', 'p50', 'p75', 'p95'].forEach(function (key) {
        ctx.beginPath();
        ctx.arc(cx, Y(c[key]), key === 'p50' ? 4 : 3, 0, 6.2832);
        ctx.fillStyle = 'rgb(' + (key === 'p50' ? b : a) + ')';
        ctx.fill();
      });
      var t = (kk / STEPS).toFixed(2);
      var txt = 't=' + t + '   median ' + Math.exp(c.p50).toFixed(2) +
                '   50%: ' + Math.exp(c.p25).toFixed(2) + '–' + Math.exp(c.p75).toFixed(2) +
                '   90%: ' + Math.exp(c.p05).toFixed(2) + '–' + Math.exp(c.p95).toFixed(2);
      ctx.font = '600 12px ui-monospace, SFMono-Regular, Menlo, monospace';
      var tw = ctx.measureText(txt).width + 20;
      var bx = Math.min(w - tw - 10, Math.max(x0, cx - tw / 2)), by = 14;
      ctx.fillStyle = css('--surface', '#fff');
      ctx.globalAlpha = 0.92;
      ctx.fillRect(bx, by, tw, 26);
      ctx.globalAlpha = 1;
      ctx.strokeStyle = 'rgba(' + a + ',0.4)';
      ctx.strokeRect(bx + 0.5, by + 0.5, tw - 1, 25);
      ctx.fillStyle = css('--text', '#111');
      ctx.fillText(txt, bx + 10, by + 17);
    }
  }

  function schedule() {
    if (raf) return;
    raf = requestAnimationFrame(function () { raf = 0; draw(); });
  }

  function play() {
    if (reduce) { progress = 1; draw(); return; }
    progress = 0; animating = true;
    (function frame() {
      progress = Math.min(1, progress + 0.012);
      draw();
      if (progress < 1) requestAnimationFrame(frame); else animating = false;
    })();
  }

  hero.addEventListener('mousemove', function (e) {
    var r = hero.getBoundingClientRect();
    mouse = { x: e.clientX - r.left, y: e.clientY - r.top };
    if (!animating) schedule();
  });
  hero.addEventListener('mouseleave', function () { mouse = null; if (!animating) schedule(); });
  hero.addEventListener('click', function (e) {
    if (e.target.closest('a, button')) return;
    simulate(); play();
  });

  var timer;
  window.addEventListener('resize', function () {
    clearTimeout(timer);
    timer = setTimeout(function () { resize(); draw(); }, 120);
  });
  if (window.matchMedia) {
    var m = matchMedia('(prefers-color-scheme: dark)');
    (m.addEventListener ? m.addEventListener.bind(m, 'change') : m.addListener.bind(m))(schedule);
  }

  resize(); simulate(); play();
})();
