/* ROAR-DB · table.js
 *
 * One implementation of the table behaviour the whole site needs, replacing the
 * seven near-identical inline copies that used to live in the templates.
 *
 * Applies to every <table class="data"> on the page:
 *   · click (or Enter/Space) on a header sorts by that column, toggling asc/desc
 *   · the active column shows its direction, and reports it through aria-sort
 *   · numbers sort as numbers, using data-v when the cell carries one
 *   · empty cells and dashes always sink to the bottom, in both directions
 *   · long tables get a capped, scrollable viewport so the sticky header works
 *
 * Opt-in extras, declared in the template:
 *   <table class="data" data-filter>                 text filter + "showing N of M"
 *   <table class="data" data-filter="placeholder">   ... with a custom placeholder
 *   <th data-facet>            that column also becomes a value filter
 *   <th data-facet="Label">    ... with an explicit legend
 *   <table data-nosort> / <th data-nosort>           exclude from sorting
 *
 * No dependencies.
 */
(function () {
  'use strict';

  /* A cell holding only one of these counts as "no value": it must never
     outrank a real measurement, whichever direction the column is sorted. */
  var BLANK = /^(|[‐-―-]|n\/?a|none|null)$/i;
  var NUMERIC = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/;
  var LONG_TABLE = 24;          /* rows beyond which a table gets its own scroll box */

  function text(el) {
    return (el && el.textContent || '').replace(/\s+/g, ' ').trim();
  }

  /* Sort key for one cell.
     data-v wins when present, because the templates use it to expose the raw
     number behind a formatted cell ("6.0 %", "1,809", "4.20x"). Otherwise the
     visible text is parsed after stripping thousands separators, percent signs
     and multiplication crosses; scientific notation (7.2e-33) is kept, which
     matters because the statistics pages print p values that way. */
  function keyOf(td) {
    if (!td) return { blank: true, n: null, s: '' };
    var raw = (td.dataset && td.dataset.v != null && td.dataset.v !== '')
      ? td.dataset.v : text(td);
    var s = String(raw).trim();
    if (BLANK.test(s)) return { blank: true, n: null, s: '' };
    var t = s.replace(/[\s,  ]/g, '').replace(/[%×x]+$/i, '');
    var n = NUMERIC.test(t) ? parseFloat(t) : NaN;
    return { blank: false, n: isNaN(n) ? null : n, s: s.toLowerCase() };
  }

  function compare(a, b, dir) {
    if (a.blank && b.blank) return 0;
    if (a.blank) return 1;                 /* blanks last, both directions */
    if (b.blank) return -1;
    if (a.n !== null && b.n !== null) return (a.n - b.n) * dir;
    if (a.n !== null) return -dir;         /* numbers before free text */
    if (b.n !== null) return dir;
    return a.s.localeCompare(b.s) * dir;
  }

  /* ------------------------------------------------------------------ sorting */

  function headerRow(table) {
    var thead = table.tHead;
    if (!thead || !thead.rows.length) return null;
    /* Residue-signature tables stack two header rows; sorting them would
       separate a variant from its own residues, so they are left alone. */
    if (thead.rows.length > 1) return null;
    var row = thead.rows[0];
    for (var i = 0; i < row.cells.length; i++) {
      if (row.cells[i].colSpan > 1 || row.cells[i].rowSpan > 1) return null;
    }
    return row;
  }

  /* Rows that do not line up with the header - a "No entries." cell spanning the
     whole width, or a sub-total - are parked at the bottom and never reordered. */
  function sortableRows(tbody, ncols) {
    var keep = [], park = [];
    for (var i = 0; i < tbody.rows.length; i++) {
      var r = tbody.rows[i], ok = r.cells.length >= ncols;
      for (var j = 0; ok && j < r.cells.length; j++) {
        if (r.cells[j].colSpan > 1) ok = false;
      }
      (ok ? keep : park).push(r);
    }
    return { keep: keep, park: park };
  }

  function enableSort(table) {
    if (table.hasAttribute('data-nosort')) return;
    var row = headerRow(table);
    if (!row || !table.tBodies.length) return;
    var tbody = table.tBodies[0];
    if (tbody.rows.length < 2) return;
    var ncols = row.cells.length;

    Array.prototype.forEach.call(row.cells, function (th, index) {
      if (th.hasAttribute('data-nosort')) return;
      th.classList.add('sortable');
      th.setAttribute('aria-sort', 'none');
      if (!th.hasAttribute('tabindex')) th.tabIndex = 0;
      if (!th.hasAttribute('role')) th.setAttribute('role', 'button');
      if (!th.title) th.title = 'Sort by ' + (text(th) || 'this column');

      function run() {
        /* Default to descending for numeric columns: on this site those are
           counts and shares, where the interesting end is the top. */
        var first = th.getAttribute('aria-sort') === 'none';
        var numeric = th.classList.contains('num');
        var dir;
        if (first) dir = numeric ? -1 : 1;
        else dir = th.getAttribute('aria-sort') === 'ascending' ? -1 : 1;

        var split = sortableRows(tbody, ncols);
        var decorated = split.keep.map(function (r, i) {
          return { row: r, key: keyOf(r.cells[index]), i: i };
        });
        decorated.sort(function (a, b) {
          return compare(a.key, b.key, dir) || (a.i - b.i);   /* stable */
        });
        var frag = document.createDocumentFragment();
        decorated.forEach(function (d) { frag.appendChild(d.row); });
        split.park.forEach(function (r) { frag.appendChild(r); });
        tbody.appendChild(frag);

        Array.prototype.forEach.call(row.cells, function (other) {
          if (other !== th) {
            other.setAttribute('aria-sort', 'none');
            other.classList.remove('sorted');
          }
        });
        th.setAttribute('aria-sort', dir === 1 ? 'ascending' : 'descending');
        th.classList.add('sorted');
      }

      th.addEventListener('click', run);
      th.addEventListener('keydown', function (e) {
        if (e.key === 'Enter' || e.key === ' ' || e.key === 'Spacebar') {
          e.preventDefault();
          run();
        }
      });
    });
  }

  /* ----------------------------------------------------------------- filtering */

  function facetValue(td) {
    if (!td) return '';
    if (td.dataset && td.dataset.f != null) return td.dataset.f.trim();
    return text(td);
  }

  function enableFilter(table) {
    if (!table.hasAttribute('data-filter')) return;
    var row = headerRow(table);
    if (!row || !table.tBodies.length) return;
    var tbody = table.tBodies[0];
    var ncols = row.cells.length;
    var rows = sortableRows(tbody, ncols).keep;
    if (rows.length < 4) return;

    /* Cache each row's searchable text once. */
    var haystack = rows.map(function (r) { return text(r).toLowerCase(); });

    var facets = [];
    Array.prototype.forEach.call(row.cells, function (th, index) {
      if (!th.hasAttribute('data-facet')) return;
      var counts = Object.create(null), order = [];
      rows.forEach(function (r) {
        var v = facetValue(r.cells[index]);
        if (!v || BLANK.test(v)) return;
        if (!(v in counts)) { counts[v] = 0; order.push(v); }
        counts[v]++;
      });
      if (order.length < 2 || order.length > 16) return;
      order.sort(function (a, b) {
        return (counts[b] - counts[a]) || a.localeCompare(b);
      });
      facets.push({
        index: index,
        label: th.getAttribute('data-facet') || text(th),
        values: order,
        counts: counts,
        active: ''
      });
    });

    /* ---- toolbar ---- */
    var bar = document.createElement('div');
    bar.className = 'tbar';

    var field = document.createElement('div');
    field.className = 'tbar__field';
    var input = document.createElement('input');
    input.type = 'search';
    input.className = 'tbar__q';
    var ph = table.getAttribute('data-filter');
    input.placeholder = (ph && ph !== 'true' && ph !== '') ? ph : 'Filter rows…';
    input.setAttribute('aria-label', input.placeholder);
    var count = document.createElement('span');
    count.className = 'tbar__count';
    count.setAttribute('role', 'status');
    count.setAttribute('aria-live', 'polite');
    var clear = document.createElement('button');
    clear.type = 'button';
    clear.className = 'tbar__clear';
    clear.textContent = 'Reset';
    clear.hidden = true;
    field.appendChild(input);
    field.appendChild(count);
    field.appendChild(clear);
    bar.appendChild(field);

    facets.forEach(function (f) {
      var group = document.createElement('div');
      group.className = 'tbar__facet';
      var legend = document.createElement('span');
      legend.className = 'tbar__legend';
      legend.textContent = f.label;
      group.appendChild(legend);

      if (f.values.length > 7) {
        /* Too many for chips: a select keeps the toolbar one line tall. */
        var sel = document.createElement('select');
        sel.className = 'tbar__select';
        sel.setAttribute('aria-label', f.label);
        var all = document.createElement('option');
        all.value = '';
        all.textContent = 'all (' + f.values.length + ')';
        sel.appendChild(all);
        f.values.forEach(function (v) {
          var o = document.createElement('option');
          o.value = v;
          o.textContent = v + ' · ' + f.counts[v];
          sel.appendChild(o);
        });
        sel.addEventListener('change', function () {
          f.active = sel.value;
          apply();
        });
        f.reset = function () { sel.value = ''; };
        group.appendChild(sel);
      } else {
        f.chips = [];
        f.values.forEach(function (v) {
          var b = document.createElement('button');
          b.type = 'button';
          b.className = 'facetchip';
          b.setAttribute('aria-pressed', 'false');
          b.innerHTML = '';
          b.appendChild(document.createTextNode(v));
          var n = document.createElement('span');
          n.className = 'facetchip__n';
          n.textContent = f.counts[v];
          b.appendChild(n);
          b.addEventListener('click', function () {
            f.active = (f.active === v) ? '' : v;
            apply();
          });
          f.chips.push({ el: b, value: v });
          group.appendChild(b);
        });
        f.reset = function () { f.active = ''; };
      }
      bar.appendChild(group);
    });

    var wrap = table.closest('.tablewrap') || table;
    wrap.parentNode.insertBefore(bar, wrap);

    function apply() {
      var terms = input.value.toLowerCase().split(/\s+/).filter(Boolean);
      var shown = 0;
      rows.forEach(function (r, i) {
        var ok = true;
        for (var t = 0; ok && t < terms.length; t++) {
          if (haystack[i].indexOf(terms[t]) === -1) ok = false;
        }
        for (var f = 0; ok && f < facets.length; f++) {
          if (facets[f].active && facetValue(r.cells[facets[f].index]) !== facets[f].active) {
            ok = false;
          }
        }
        r.hidden = !ok;
        /* hidden alone is overridden by `display: table-row` in some sheets */
        r.style.display = ok ? '' : 'none';
        if (ok) shown++;
      });

      facets.forEach(function (f) {
        if (!f.chips) return;
        f.chips.forEach(function (c) {
          var on = f.active === c.value;
          c.el.classList.toggle('on', on);
          c.el.setAttribute('aria-pressed', on ? 'true' : 'false');
        });
      });

      var filtering = terms.length > 0 || facets.some(function (f) { return !!f.active; });
      count.textContent = filtering
        ? 'showing ' + shown.toLocaleString() + ' of ' + rows.length.toLocaleString()
        : rows.length.toLocaleString() + ' rows';
      count.classList.toggle('on', filtering);
      clear.hidden = !filtering;

      var none = table.querySelector('tr.tbar__none');
      if (shown === 0 && !none) {
        none = tbody.insertRow();
        none.className = 'tbar__none';
        var td = none.insertCell();
        td.colSpan = ncols;
        td.className = 'empty';
        td.textContent = 'No row matches this filter.';
      } else if (none) {
        none.style.display = shown === 0 ? '' : 'none';
      }
    }

    input.addEventListener('input', apply);
    clear.addEventListener('click', function () {
      input.value = '';
      facets.forEach(function (f) { f.active = ''; if (f.reset) f.reset(); });
      apply();
      input.focus();
    });
    apply();
  }

  /* --------------------------------------------- sticky headers that stick */

  /* `position: sticky` resolves against the nearest scrolling ancestor. Because
     .tablewrap is an overflow container, a sticky header inside it pins to the
     top of that box - which, with the box as tall as its content, is no pin at
     all. Giving long tables a capped viewport turns the wrapper into a real
     scroll area, and the header then behaves as intended. */
  function stickyHeader(table) {
    var thead = table.tHead;
    if (!thead || !thead.rows.length) return;
    /* Screen readers need to know these cells head their columns. */
    Array.prototype.forEach.call(thead.querySelectorAll('th'), function (th) {
      if (!th.hasAttribute('scope')) th.setAttribute('scope', 'col');
    });
    var wrap = table.closest('.tablewrap');
    var body = table.tBodies.length ? table.tBodies[0] : null;
    if (wrap && body && body.rows.length > LONG_TABLE) {
      wrap.classList.add('tablewrap--tall');
    }
    /* Stacked header rows each need their own offset. */
    if (thead.rows.length > 1) {
      var top = 0;
      Array.prototype.forEach.call(thead.rows, function (r) {
        Array.prototype.forEach.call(r.cells, function (c) { c.style.top = top + 'px'; });
        top += r.offsetHeight;
      });
    }
  }

  /* ------------------------------------------------- card grids (home page) */

  /* The home page presents the enzyme types as cards grouped by chemical
     family rather than as one table, so it needs the same find-a-type
     behaviour applied to a card grid. The container declares which data
     attributes are facets; cards carry the values.

       <div data-cardfilter="Find a type by name, substrate or product"
            data-card-facets="reaction:Reaction,sclass:Substrate class,group:RO group">
         ... <a class="structcard" data-reaction="hydroxylation" ...>
  */
  function enableCardFilter(root) {
    var cards = Array.prototype.slice.call(root.querySelectorAll('[data-card]'));
    if (cards.length < 4) return;

    var spec = (root.getAttribute('data-card-facets') || '').split(',');
    var facets = [];
    spec.forEach(function (part) {
      part = part.trim();
      if (!part) return;
      var bits = part.split(':');
      var attr = bits[0].trim();
      if (!attr) return;
      var counts = Object.create(null), order = [];
      cards.forEach(function (c) {
        var v = (c.getAttribute('data-' + attr) || '').trim();
        if (!v || BLANK.test(v)) return;
        if (!(v in counts)) { counts[v] = 0; order.push(v); }
        counts[v]++;
      });
      if (order.length < 2 || order.length > 16) return;
      order.sort(function (a, b) { return (counts[b] - counts[a]) || a.localeCompare(b); });
      facets.push({
        attr: attr,
        label: (bits[1] || attr).trim(),
        values: order,
        counts: counts,
        active: ''
      });
    });

    var haystack = cards.map(function (c) {
      return ((c.getAttribute('data-k') || '') + ' ' + text(c)).toLowerCase();
    });
    /* Sections that end up with no visible card are hidden too, so the page
       does not keep a run of empty family headings. */
    var sections = Array.prototype.slice.call(root.querySelectorAll('[data-cardsection]'));

    var bar = document.createElement('div');
    bar.className = 'tbar tbar--cards';
    var field = document.createElement('div');
    field.className = 'tbar__field';
    var input = document.createElement('input');
    input.type = 'search';
    input.className = 'tbar__q';
    var ph = root.getAttribute('data-cardfilter');
    input.placeholder = (ph && ph !== 'true') ? ph : 'Filter…';
    input.setAttribute('aria-label', input.placeholder);
    var count = document.createElement('span');
    count.className = 'tbar__count';
    count.setAttribute('role', 'status');
    count.setAttribute('aria-live', 'polite');
    var clear = document.createElement('button');
    clear.type = 'button';
    clear.className = 'tbar__clear';
    clear.textContent = 'Reset';
    clear.hidden = true;
    field.appendChild(input);
    field.appendChild(count);
    field.appendChild(clear);
    bar.appendChild(field);

    facets.forEach(function (f) {
      var group = document.createElement('div');
      group.className = 'tbar__facet';
      var legend = document.createElement('span');
      legend.className = 'tbar__legend';
      legend.textContent = f.label;
      group.appendChild(legend);
      f.chips = [];
      f.values.forEach(function (v) {
        var b = document.createElement('button');
        b.type = 'button';
        b.className = 'facetchip';
        b.setAttribute('aria-pressed', 'false');
        b.appendChild(document.createTextNode(v));
        var n = document.createElement('span');
        n.className = 'facetchip__n';
        n.textContent = f.counts[v];
        b.appendChild(n);
        b.addEventListener('click', function () {
          f.active = (f.active === v) ? '' : v;
          apply();
        });
        f.chips.push({ el: b, value: v });
        group.appendChild(b);
      });
      bar.appendChild(group);
    });

    root.insertBefore(bar, root.firstChild);

    function apply() {
      var terms = input.value.toLowerCase().split(/\s+/).filter(Boolean);
      var shown = 0;
      cards.forEach(function (c, i) {
        var ok = true;
        for (var t = 0; ok && t < terms.length; t++) {
          if (haystack[i].indexOf(terms[t]) === -1) ok = false;
        }
        for (var f = 0; ok && f < facets.length; f++) {
          if (facets[f].active &&
              (c.getAttribute('data-' + facets[f].attr) || '').trim() !== facets[f].active) {
            ok = false;
          }
        }
        c.style.display = ok ? '' : 'none';
        if (ok) shown++;
      });
      sections.forEach(function (s) {
        var any = s.querySelector('[data-card]:not([style*="display: none"])');
        s.style.display = any ? '' : 'none';
      });
      facets.forEach(function (f) {
        f.chips.forEach(function (c) {
          var on = f.active === c.value;
          c.el.classList.toggle('on', on);
          c.el.setAttribute('aria-pressed', on ? 'true' : 'false');
        });
      });
      var filtering = terms.length > 0 || facets.some(function (f) { return !!f.active; });
      count.textContent = filtering
        ? 'showing ' + shown + ' of ' + cards.length + ' types'
        : cards.length + ' types';
      count.classList.toggle('on', filtering);
      clear.hidden = !filtering;
    }

    input.addEventListener('input', apply);
    clear.addEventListener('click', function () {
      input.value = '';
      facets.forEach(function (f) { f.active = ''; });
      apply();
      input.focus();
    });
    apply();
  }

  function init() {
    var tables = document.querySelectorAll('table.data');
    Array.prototype.forEach.call(tables, function (t) {
      try {
        enableSort(t);
        enableFilter(t);
        stickyHeader(t);
      } catch (e) {
        /* A single malformed table must not take the rest of the page down. */
        if (window.console && console.warn) console.warn('table.js', e);
      }
    });
    Array.prototype.forEach.call(document.querySelectorAll('[data-cardfilter]'), function (r) {
      try {
        enableCardFilter(r);
      } catch (e) {
        if (window.console && console.warn) console.warn('table.js cards', e);
      }
    });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
