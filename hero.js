// Hero: MCMC samplers exploring a 2D Bayesian posterior.
// Chains start scattered and converge on the modes (burn-in, live).
// Move the cursor to drag one mode of the posterior; click to release fresh chains.
(function () {
  var hero = document.querySelector('.hero');
  var canvas = document.getElementById('paths');
  if (!hero || !canvas || !canvas.getContext) return;
  var ctx = canvas.getContext('2d');
  var reduce = window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches;

  var GW = 260, GH = 156;
  var dens = document.createElement('canvas'); dens.width = GW; dens.height = GH;
  var dctx = dens.getContext('2d');
  var img = dctx.createImageData(GW, GH);
  var trail = document.createElement('canvas');
  var tctx = trail.getContext('2d');

  var w = 0, h = 0, px0 = 0, pw = 0;
  var colA = [31, 95, 214], colB = [124, 58, 237];
  var mouse = null, visible = true, raf = 0, dirty = true;
  var prop = 0, accd = 0, tick = 0;
  var accEl = document.getElementById('acc'), nEl = document.getElementById('nchains');

  function copy(m) { return { x: m.x, y: m.y, sx: m.sx, sy: m.sy, rho: m.rho, wt: m.wt }; }
  var HOME = [
    { x: 0.30, y: 0.42, sx: 0.11, sy: 0.17, rho: 0.55, wt: 1.0 },
    { x: 0.72, y: 0.68, sx: 0.15, sy: 0.10, rho: -0.45, wt: 0.75 },
    { x: 0.62, y: 0.22, sx: 0.09, sy: 0.09, rho: 0.0, wt: 0.65 }
  ];
  var modes = HOME.map(copy);
  var follow = { x: HOME[2].x, y: HOME[2].y };

  function density(u, v) {
    var s = 0;
    for (var i = 0; i < 3; i++) {
      var m = modes[i], dx = (u - m.x) / m.sx, dy = (v - m.y) / m.sy;
      var q = (dx * dx - 2 * m.rho * dx * dy + dy * dy) / (2 * (1 - m.rho * m.rho));
      s += m.wt * Math.exp(-q);
    }
    return s;
  }

  function randn() {
    var u = 0, v = 0;
    while (u === 0) u = Math.random();
    while (v === 0) v = Math.random();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  var NC = 42, chains = [];
  function scatter(c) { c.u = Math.random(); c.v = Math.random(); c.p = density(c.u, c.v); }
  function initChains() {
    chains = [];
    for (var i = 0; i < NC; i++) { var c = {}; scatter(c); chains.push(c); }
  }

  function X(u) { return px0 + u * pw; }
  function Y(v) { return v * h; }

  function step() {
    var rgb = 'rgba(' + colB.join(',') + ',0.55)';
    tctx.strokeStyle = rgb;
    tctx.lineWidth = 1;
    for (var i = 0; i < chains.length; i++) {
      var c = chains[i];
      for (var r = 0; r < 2; r++) {
        var nu = c.u + randn() * 0.045, nv = c.v + randn() * 0.045;
        if (nu < 0 || nu > 1 || nv < 0 || nv > 1) continue;
        var np = density(nu, nv);
        prop++;
        if (Math.random() * c.p < np) {
          accd++;
          tctx.beginPath(); tctx.moveTo(X(c.u), Y(c.v)); tctx.lineTo(X(nu), Y(nv)); tctx.stroke();
          c.u = nu; c.v = nv; c.p = np;
        }
      }
      if (Math.random() < 0.0015) scatter(c); // occasional fresh chain: burn-in on display
    }
  }

  function updateModes() {
    var tx = HOME[2].x, ty = HOME[2].y;
    if (mouse && mouse.x >= px0 - 40) {
      tx = Math.min(1, Math.max(0, (mouse.x - px0) / pw));
      ty = Math.min(1, Math.max(0, mouse.y / h));
    }
    var m = modes[2], dx = tx - m.x, dy = ty - m.y;
    if (Math.abs(dx) + Math.abs(dy) > 0.0004) { m.x += dx * 0.07; m.y += dy * 0.07; dirty = true; }
    if (dirty) { for (var i = 0; i < chains.length; i++) chains[i].p = density(chains[i].u, chains[i].v); }
  }

  function renderDensity() {
    var d = img.data, i = 0, a = colA, b = colB;
    for (var j = 0; j < GH; j++) {
      var v = (j + 0.5) / GH;
      for (var k = 0; k < GW; k++) {
        var u = (k + 0.5) / GW;
        var p = density(u, v); if (p > 1) p = 1;
        var lvl = p * 8, band = lvl - Math.floor(lvl);
        var al = Math.pow(p, 0.75) * 0.40;
        if (p > 0.04 && band < 0.11) al += 0.20 * Math.min(1, p * 3); // contour ridges
        d[i++] = a[0] + (b[0] - a[0]) * u;
        d[i++] = a[1] + (b[1] - a[1]) * u;
        d[i++] = a[2] + (b[2] - a[2]) * u;
        d[i++] = al * 255;
      }
    }
    dctx.putImageData(img, 0, 0);
    dirty = false;
  }

  function draw() {
    ctx.clearRect(0, 0, w, h);
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(dens, px0, 0, pw, h);
    ctx.drawImage(trail, 0, 0, w, h);
    for (var i = 0; i < chains.length; i++) {
      var cx = X(chains[i].u), cy = Y(chains[i].v);
      ctx.fillStyle = 'rgba(' + colB.join(',') + ',0.18)';
      ctx.beginPath(); ctx.arc(cx, cy, 6, 0, 6.2832); ctx.fill();
      ctx.fillStyle = 'rgb(' + colB.join(',') + ')';
      ctx.beginPath(); ctx.arc(cx, cy, 2.2, 0, 6.2832); ctx.fill();
    }
  }

  function frame() {
    raf = 0;
    if (!visible || document.hidden) return;
    tctx.globalCompositeOperation = 'destination-out';
    tctx.fillStyle = 'rgba(0,0,0,0.05)';
    tctx.fillRect(0, 0, w, h);
    tctx.globalCompositeOperation = 'source-over';
    updateModes();
    if (dirty) renderDensity();
    step();
    draw();
    if (++tick % 20 === 0 && accEl && prop) {
      accEl.textContent = (accd / prop).toFixed(2);
      prop *= 0.5; accd *= 0.5;
    }
    raf = requestAnimationFrame(frame);
  }
  function start() { if (!raf && !reduce) raf = requestAnimationFrame(frame); }

  function readColors() {
    function rgb(name, fb) {
      var v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
      var p = (v || fb).split(',').map(Number);
      return p.length === 3 ? p : fb.split(',').map(Number);
    }
    colA = rgb('--path', '31, 95, 214');
    colB = rgb('--path2', '124, 58, 237');
    dirty = true;
  }

  function resize() {
    var r = hero.getBoundingClientRect(), dpr = window.devicePixelRatio || 1;
    w = r.width; h = r.height;
    px0 = w > 900 ? w * 0.46 : 0; pw = w - px0;
    canvas.width = w * dpr; canvas.height = h * dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    trail.width = w * dpr; trail.height = h * dpr;
    tctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    dirty = true;
  }

  function staticFrame() {
    for (var n = 0; n < 260; n++) step();
    renderDensity(); draw();
  }

  hero.addEventListener('mousemove', function (e) {
    var r = hero.getBoundingClientRect();
    mouse = { x: e.clientX - r.left, y: e.clientY - r.top };
  });
  hero.addEventListener('mouseleave', function () { mouse = null; });
  hero.addEventListener('click', function (e) {
    if (e.target.closest('a, button')) return;
    chains.forEach(scatter);
  });
  if ('IntersectionObserver' in window) {
    new IntersectionObserver(function (es) {
      visible = es[0].isIntersecting;
      if (visible) start();
    }).observe(hero);
  }
  document.addEventListener('visibilitychange', function () { if (!document.hidden) start(); });

  var timer;
  window.addEventListener('resize', function () {
    clearTimeout(timer);
    timer = setTimeout(function () { resize(); if (reduce) staticFrame(); }, 120);
  });
  if (window.matchMedia) {
    var mq = matchMedia('(prefers-color-scheme: dark)');
    var onScheme = function () { readColors(); if (reduce) staticFrame(); };
    if (mq.addEventListener) mq.addEventListener('change', onScheme); else mq.addListener(onScheme);
  }

  document.addEventListener('themechange', function () { readColors(); if (reduce) staticFrame(); });
  if (nEl) nEl.textContent = NC;

  readColors(); resize(); initChains();
  if (reduce) staticFrame(); else start();
})();
