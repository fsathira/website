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

  // --- Load Data ---
  fetch('data/publications.json')
    .then(function (res) { return res.json(); })
    .then(function (data) {
      publications = data;
      buildFilterChips();
      renderPublications();
    })
    .catch(function (err) {
      pubList.innerHTML = '<div class="pub-empty">Failed to load publications. Please try refreshing.</div>';
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
    var displayTags = ['genomics', 'microbiome', 'bioinformatics', 'drug-discovery', 'statistical-methods', 'AI-drug-development'];
    var availableTags = displayTags.filter(function (t) { return tagSet[t]; });

    availableTags.forEach(function (tag) {
      var chip = document.createElement('button');
      chip.className = 'pub-filter-chip';
      chip.setAttribute('data-filter', tag);
      chip.textContent = tag.replace(/-/g, ' ');
      pubFilters.appendChild(chip);
    });

    // Event listeners for all chips
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

      return '<div class="pub-card">' +
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
    return links.join('');
  }

  function escapeHtml(str) {
    var div = document.createElement('div');
    div.appendChild(document.createTextNode(str));
    return div.innerHTML;
  }
})();
