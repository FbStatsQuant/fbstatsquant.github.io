// Small UI behaviours: scroll reveal, active nav link, copy-email toast.
(function () {
  // Scroll reveal
  var items = document.querySelectorAll('.reveal');
  if ('IntersectionObserver' in window) {
    var io = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (e.isIntersecting) { e.target.classList.add('in'); io.unobserve(e.target); }
      });
    }, { threshold: 0.12 });
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
        if (e.isIntersecting) {
          Object.keys(links).forEach(function (k) { links[k].classList.remove('active'); });
          var l = links['#' + e.target.id];
          if (l) l.classList.add('active');
        }
      });
    }, { rootMargin: '-45% 0px -50% 0px' });
    document.querySelectorAll('main section[id]').forEach(function (s) { so.observe(s); });
  }

  // Copy email on click (mailto still fires; this covers machines with no mail app)
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
          function () { show('Copied ' + text); },
          function () { show(text); }
        );
      } else {
        show(text);
      }
    });
  });
})();
