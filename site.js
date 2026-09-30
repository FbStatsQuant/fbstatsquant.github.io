// UI behaviours: theme toggle, scroll progress, nav state, reveals, card spotlight, copy email.
(function () {
  var root = document.documentElement;

  // Theme toggle (persists per visitor; falls back to OS preference)
  var toggle = document.querySelector('.theme-toggle');
  function current() {
    var t = root.getAttribute('data-theme');
    if (t) return t;
    return window.matchMedia && matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
  }
  if (toggle) toggle.addEventListener('click', function () {
    var next = current() === 'dark' ? 'light' : 'dark';
    root.setAttribute('data-theme', next);
    try { localStorage.setItem('theme', next); } catch (e) {}
    document.dispatchEvent(new CustomEvent('themechange'));
  });

  // Scroll progress + nav border
  var bar = document.querySelector('.progress');
  var nav = document.querySelector('.nav');
  var ticking = false;
  function onScroll() {
    ticking = false;
    var max = root.scrollHeight - innerHeight;
    if (bar) bar.style.transform = 'scaleX(' + (max > 0 ? scrollY / max : 0) + ')';
    if (nav) nav.classList.toggle('scrolled', scrollY > 8);
  }
  addEventListener('scroll', function () { if (!ticking) { ticking = true; requestAnimationFrame(onScroll); } }, { passive: true });
  onScroll();

  // Reveal on scroll
  var items = document.querySelectorAll('.reveal');
  if ('IntersectionObserver' in window) {
    var io = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (e.isIntersecting) { e.target.classList.add('in'); io.unobserve(e.target); }
      });
    }, { threshold: 0.1, rootMargin: '0px 0px -40px 0px' });
    items.forEach(function (el) { io.observe(el); });
  } else {
    items.forEach(function (el) { el.classList.add('in'); });
  }

  // Active nav link
  var links = {};
  document.querySelectorAll('.nav nav a').forEach(function (a) { links[a.getAttribute('href')] = a; });
  if ('IntersectionObserver' in window) {
    var so = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (!e.isIntersecting) return;
        Object.keys(links).forEach(function (k) { links[k].classList.remove('active'); });
        var l = links['#' + e.target.id];
        if (l) {
          l.classList.add('active');
          if (l.scrollIntoView && l.parentNode.scrollWidth > l.parentNode.clientWidth) {
            l.parentNode.scrollTo({ left: l.offsetLeft - 40, behavior: 'smooth' });
          }
        }
      });
    }, { rootMargin: '-45% 0px -50% 0px' });
    document.querySelectorAll('main section[id]').forEach(function (s) { so.observe(s); });
  }

  // Cursor spotlight on cards
  document.querySelectorAll('.paper, .proj').forEach(function (el) {
    el.addEventListener('pointermove', function (e) {
      var r = el.getBoundingClientRect();
      el.style.setProperty('--mx', (e.clientX - r.left) + 'px');
      el.style.setProperty('--my', (e.clientY - r.top) + 'px');
    });
  });

  // Copy email on click (mailto still fires; covers machines with no mail app)
  var toast = document.getElementById('toast');
  var t;
  function show(msg) {
    if (!toast) return;
    toast.textContent = msg;
    toast.classList.add('show');
    clearTimeout(t);
    t = setTimeout(function () { toast.classList.remove('show'); }, 2400);
  }
  document.querySelectorAll('[data-copy]').forEach(function (a) {
    a.addEventListener('click', function () {
      var text = a.getAttribute('data-copy');
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(text).then(
          function () { show('✓ Copied ' + text); },
          function () { show(text); }
        );
      } else {
        show(text);
      }
    });
  });

  var y = document.getElementById('year');
  if (y) y.textContent = new Date().getFullYear();
})();
