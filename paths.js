// Hero background: simulated geometric Brownian motion paths, drawn once (static if reduced motion).
(function () {
  var canvas = document.getElementById('paths');
  if (!canvas || !canvas.getContext) return;
  var ctx = canvas.getContext('2d');
  var reduce = window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches;
  var w, h, paths, t;
  var N = 90, STEPS = 220;

  function randn() {
    var u = 0, v = 0;
    while (u === 0) u = Math.random();
    while (v === 0) v = Math.random();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  function build() {
    paths = [];
    var mu = 0.05, sigma = 0.35, dt = 1 / STEPS;
    for (var i = 0; i < N; i++) {
      var s = 1, p = [1];
      for (var k = 1; k <= STEPS; k++) {
        s *= Math.exp((mu - 0.5 * sigma * sigma) * dt + sigma * Math.sqrt(dt) * randn());
        p.push(s);
      }
      paths.push(p);
    }
  }

  function resize() {
    var r = canvas.getBoundingClientRect();
    var dpr = window.devicePixelRatio || 1;
    w = r.width; h = r.height;
    canvas.width = w * dpr; canvas.height = h * dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }

  function color() {
    var v = getComputedStyle(document.documentElement).getPropertyValue('--path').trim();
    return v || '31, 95, 214';
  }

  function draw(progress) {
    ctx.clearRect(0, 0, w, h);
    var rgb = color();
    var upto = Math.max(2, Math.floor(progress * STEPS));
    var lo = 0.35, hi = 2.6; // vertical range in units of S0
    ctx.lineWidth = 1;
    for (var i = 0; i < N; i++) {
      var p = paths[i];
      ctx.beginPath();
      for (var k = 0; k <= upto; k++) {
        var x = (k / STEPS) * w;
        var y = h - ((Math.min(Math.max(p[k], lo), hi) - lo) / (hi - lo)) * h;
        if (k === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }
      ctx.strokeStyle = 'rgba(' + rgb + ',' + (0.05 + 0.10 * ((i % 7) / 6)) + ')';
      ctx.stroke();
    }
  }

  function start() {
    resize(); build();
    if (reduce) { draw(1); return; }
    t = 0;
    (function frame() {
      t += 0.006;
      draw(Math.min(t, 1));
      if (t < 1) requestAnimationFrame(frame);
    })();
  }

  var timer;
  window.addEventListener('resize', function () {
    clearTimeout(timer);
    timer = setTimeout(function () { resize(); draw(1); }, 150);
  });
  start();
})();
