/* ============================================
   Publications — Search, Filter, Sort, Render
   ============================================ */

(function () {
  'use strict';

  var pubList = document.getElementById('pubList');
  var pubEmpty = document.getElementById('pubEmpty');
  var pubCount = document.getElementById('pubCount');
  var pubFilters = document.getElementById('pubFilters');
  var searchInput = document.getElementById('pubSearchInput');
  var sortSelect = document.getElementById('pubSortSelect');

  if (!pubList) return;

  var publications = [];
  var activeFilter = 'all';
  var activeRole = 'all';

  // --- Load Data ---
  fetch('data/publications.json')
    .then(function (res) { return res.json(); })
    .then(function (data) {
      publications = data;
      buildFilterChips();
      renderPublications();
    })
    .catch(function (err) {
      console.error('Publications load error:', err);
      pubList.innerHTML = '<div class="pub-empty">Failed to load publications (' + (err.message || err) + '). Please try refreshing.</div>';
    });

  // --- Build Filter Chips from Tags ---
  function buildFilterChips() {
    var tagSet = {};
    publications.forEach(function (pub) {
      if (pub.tags) {
        pub.tags.forEach(function (tag) { tagSet[tag] = true; });
      }
    });

    var tags = Object.keys(tagSet).sort();
    // Pick the most common/interesting tags to show as chips
    var displayTags = ['GWAS', 'microbiome', 'neurogenetics', 'bioinformatics', 'COVID-19', 'psychiatric-genetics', 'cancer', 'sleep-genetics', 'drug-discovery', 'statistical-methods'];
    var availableTags = displayTags.filter(function (t) { return tagSet[t]; });

    // Role filter buttons
    var roleFilters = document.createElement('div');
    roleFilters.className = 'pub-role-filters';
    var roles = [
      { value: 'all', label: 'All' },
      { value: 'first', label: 'First / co-first author' },
      { value: 'named', label: 'Named author' },
      { value: 'consortium', label: 'Consortium' }
    ];
    roles.forEach(function (role) {
      var btn = document.createElement('button');
      btn.className = 'pub-role-chip' + (role.value === 'all' ? ' active' : '');
      btn.setAttribute('data-role', role.value);
      btn.textContent = role.label;
      btn.addEventListener('click', function () {
        roleFilters.querySelectorAll('.pub-role-chip').forEach(function (c) {
          c.classList.remove('active');
        });
        btn.classList.add('active');
        activeRole = role.value;
        renderPublications();
      });
      roleFilters.appendChild(btn);
    });
    pubFilters.parentNode.insertBefore(roleFilters, pubFilters);

    availableTags.forEach(function (tag) {
      var chip = document.createElement('button');
      chip.className = 'pub-filter-chip';
      chip.setAttribute('data-filter', tag);
      chip.textContent = tag.replace(/-/g, ' ');
      pubFilters.appendChild(chip);
    });

    // Event listeners for tag chips
    pubFilters.querySelectorAll('.pub-filter-chip').forEach(function (chip) {
      chip.addEventListener('click', function () {
        pubFilters.querySelectorAll('.pub-filter-chip').forEach(function (c) {
          c.classList.remove('active');
        });
        chip.classList.add('active');
        activeFilter = chip.getAttribute('data-filter');
        renderPublications();
      });
    });
  }

  // --- Search ---
  if (searchInput) {
    var debounceTimer;
    searchInput.addEventListener('input', function () {
      clearTimeout(debounceTimer);
      debounceTimer = setTimeout(renderPublications, 200);
    });
  }

  // --- Sort ---
  if (sortSelect) {
    sortSelect.addEventListener('change', renderPublications);
  }

  // --- Filter & Sort Logic ---
  function getFilteredPublications() {
    var query = searchInput ? searchInput.value.toLowerCase().trim() : '';
    var sort = sortSelect ? sortSelect.value : 'year-desc';

    var filtered = publications.filter(function (pub) {
      // Role filter
      if (activeRole !== 'all') {
        var role = pub.role || 'consortium';
        if (activeRole === 'first') {
          if (role !== 'first_author' && role !== 'co_first_author') return false;
        } else if (activeRole === 'named') {
          if (role !== 'named_author') return false;
        } else if (activeRole === 'consortium') {
          if (role !== 'consortium') return false;
        }
      }

      // Tag filter
      if (activeFilter !== 'all') {
        if (!pub.tags || pub.tags.indexOf(activeFilter) === -1) return false;
      }

      // Search
      if (query) {
        var searchable = [
          pub.title,
          pub.authors,
          pub.journal,
          pub.abstract || '',
          (pub.tags || []).join(' '),
          String(pub.year)
        ].join(' ').toLowerCase();

        return searchable.indexOf(query) !== -1;
      }

      return true;
    });

    // Sort
    filtered.sort(function (a, b) {
      switch (sort) {
        case 'year-asc': return a.year - b.year;
        case 'title-asc': return a.title.localeCompare(b.title);
        case 'year-desc':
        default: return b.year - a.year;
      }
    });

    return filtered;
  }

  // --- Render ---
  function renderPublications() {
    var filtered = getFilteredPublications();

    pubCount.textContent = 'Showing ' + filtered.length + ' of ' + publications.length + ' publications';

    if (filtered.length === 0) {
      pubList.innerHTML = '';
      pubEmpty.style.display = 'block';
      return;
    }

    pubEmpty.style.display = 'none';

    pubList.innerHTML = filtered.map(function (pub) {
      var authors = highlightAuthor(pub.authors);
      var links = buildLinks(pub);
      var tags = (pub.tags || []).map(function (t) {
        return '<span class="pub-tag">' + escapeHtml(t.replace(/-/g, ' ')) + '</span>';
      }).join('');

      var abstractHtml = '';
      if (pub.abstract) {
        abstractHtml = '<button class="pub-abstract-toggle" onclick="this.nextElementSibling.classList.toggle(\'expanded\');this.textContent=this.nextElementSibling.classList.contains(\'expanded\')?\'\u25B2 Hide abstract\':\'\u25BC Show abstract\'">' +
          '\u25BC Show abstract</button>' +
          '<div class="pub-abstract">' + escapeHtml(pub.abstract) + '</div>';
      }

      var roleClass = (pub.role === 'first_author' || pub.role === 'co_first_author') ? ' pub-card-first' :
                      pub.role === 'named_author' ? ' pub-card-named' : '';

      return '<div class="pub-card' + roleClass + '">' +
        '<span class="pub-year-badge">' + pub.year + '</span>' +
        '<h3 class="pub-title">' + escapeHtml(pub.title) + '</h3>' +
        '<p class="pub-authors">' + authors + '</p>' +
        '<p class="pub-journal">' + escapeHtml(pub.journal) + (pub.volume ? ', ' + escapeHtml(pub.volume) : '') + '</p>' +
        (tags ? '<div class="pub-tags">' + tags + '</div>' : '') +
        '<div class="pub-links">' + links + '</div>' +
        abstractHtml +
        '</div>';
    }).join('');
  }

  // --- Helpers ---
  function highlightAuthor(authors) {
    if (!authors) return '';
    var escaped = escapeHtml(authors);
    // Highlight Sathirapongsasuti in author lists
    return escaped.replace(/(Sathirapongsasuti\s*(?:JF|J|F)?)/gi, '<span class="highlight">$1</span>');
  }

  function buildLinks(pub) {
    var links = [];
    if (pub.doi) {
      links.push('<a href="https://doi.org/' + encodeURIComponent(pub.doi) + '" target="_blank" class="pub-link">DOI</a>');
    }
    if (pub.pubmed_id) {
      links.push('<a href="https://pubmed.ncbi.nlm.nih.gov/' + encodeURIComponent(pub.pubmed_id) + '/" target="_blank" class="pub-link">PubMed</a>');
    }
    if (pub.pmcid) {
      links.push('<a href="https://www.ncbi.nlm.nih.gov/pmc/articles/' + encodeURIComponent(pub.pmcid) + '/" target="_blank" class="pub-link">PMC</a>');
    }
    return links.join('');
  }

  function escapeHtml(str) {
    var div = document.createElement('div');
    div.appendChild(document.createTextNode(str));
    return div.innerHTML;
  }
})();
