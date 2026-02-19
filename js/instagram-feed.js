/* ============================================
   Instagram Feed + News — Updates page
   ============================================ */

(function () {
  'use strict';

  // =============================================
  // INSTAGRAM FEED CONFIGURATION
  // =============================================
  // To enable live Instagram feeds, you need an Instagram Graph API access token.
  //
  // Setup steps:
  // 1. Create a Facebook Developer account at https://developers.facebook.com
  // 2. Create a new app and add the "Instagram Basic Display" product
  // 3. Add your Instagram account as a test user
  // 4. Generate a long-lived access token
  // 5. Set the token below
  //
  // Alternatively, add specific post URLs to featuredPosts for manual embedding.
  // =============================================

  var INSTAGRAM_CONFIG = {
    accounts: [
      {
        username: 'fah_shine',
        displayName: 'Fah Sathirapongsasuti',
        label: 'Personal',
        description: 'Life adventures, travel, and behind-the-scenes moments.',
        gradient: 'linear-gradient(135deg, #833ab4, #fd1d1d, #fcb045)',
        profileUrl: 'https://instagram.com/fah_shine',
        accessToken: '', // Set for live feed
        featuredPosts: [
          // Add Instagram post URLs here for manual embedding, e.g.:
          // 'https://www.instagram.com/p/POST_ID/'
        ]
      },
      {
        username: 'fahshinewine',
        displayName: 'Sunset Cellars',
        label: 'Wine',
        description: 'Winemaking journey from vineyard to bottle in Sonoma County.',
        gradient: 'linear-gradient(135deg, #7c3aed, #a855f7, #c084fc)',
        profileUrl: 'https://instagram.com/fahshinewine',
        accessToken: '', // Set for live feed
        featuredPosts: []
      }
    ],
    postsPerAccount: 6
  };

  // =============================================
  // NEWS SECTION
  // =============================================
  var newsGrid = document.getElementById('newsGrid');

  if (newsGrid) {
    fetch('data/news.json')
      .then(function (res) { return res.json(); })
      .then(function (articles) {
        renderNews(articles);
        setupNewsFilters(articles);
      })
      .catch(function () {
        newsGrid.innerHTML = '<p style="color:var(--text-muted);">Failed to load news articles.</p>';
      });
  }

  function renderNews(articles) {
    if (!newsGrid) return;

    // Category color map
    var categoryColors = {
      'Science': '#3b82f6',
      'Wine': '#a855f7',
      'Personal': '#10b981'
    };

    // Placeholder gradient backgrounds for article images
    var placeholderGradients = {
      'wine-interview': 'linear-gradient(135deg, #7c3aed, #a855f7, #c084fc)',
      'wine-enthusiast': 'linear-gradient(135deg, #dc2626, #ef4444, #f87171)',
      'biotech': 'linear-gradient(135deg, #2563eb, #3b82f6, #60a5fa)',
      'harvard': 'linear-gradient(135deg, #991b1b, #b91c1c, #dc2626)',
      'endpoints': 'linear-gradient(135deg, #0d9488, #14b8a6, #2dd4bf)',
      'podcast': 'linear-gradient(135deg, #ea580c, #f97316, #fb923c)',
      'aaas': 'linear-gradient(135deg, #1e40af, #2563eb, #3b82f6)',
      'rhythm': 'linear-gradient(135deg, #7c3aed, #6366f1, #818cf8)',
      'dress': 'linear-gradient(135deg, #1e3a8a, #c4b454, #f5f5dc)'
    };

    newsGrid.innerHTML = articles.map(function (article) {
      var gradient = placeholderGradients[article.image_placeholder] || 'linear-gradient(135deg, #374151, #4b5563)';
      var dateStr = formatDate(article.date);

      return '<article class="news-card" data-category="' + escapeAttr(article.category) + '">' +
        '<div class="news-card-image" style="background:' + gradient + ';display:flex;align-items:center;justify-content:center;">' +
          '<span style="font-size:2.5rem;opacity:0.3;color:#fff;">' + getCategoryIcon(article.category) + '</span>' +
        '</div>' +
        '<div class="news-card-body">' +
          '<div class="news-card-meta">' +
            '<span class="news-card-source">' + escapeHtml(article.source) + '</span>' +
            '<span class="news-card-date">' + dateStr + '</span>' +
            '<span class="news-card-category">' + escapeHtml(article.category) + '</span>' +
          '</div>' +
          '<h3 class="news-card-title"><a href="' + escapeAttr(article.url) + '" target="_blank">' + escapeHtml(article.title) + '</a></h3>' +
          '<p class="news-card-excerpt">' + escapeHtml(article.excerpt) + '</p>' +
          '<a href="' + escapeAttr(article.url) + '" target="_blank" class="news-card-link">Read article &rarr;</a>' +
        '</div>' +
        '</article>';
    }).join('');
  }

  function setupNewsFilters(articles) {
    var filterBtns = document.querySelectorAll('[data-news-filter]');

    filterBtns.forEach(function (btn) {
      btn.addEventListener('click', function () {
        filterBtns.forEach(function (b) { b.classList.remove('active'); });
        btn.classList.add('active');

        var filter = btn.getAttribute('data-news-filter');
        var cards = newsGrid.querySelectorAll('.news-card');

        cards.forEach(function (card) {
          if (filter === 'all' || card.getAttribute('data-category') === filter) {
            card.style.display = '';
          } else {
            card.style.display = 'none';
          }
        });
      });
    });
  }

  // =============================================
  // INSTAGRAM FEEDS
  // =============================================
  var feedsContainer = document.getElementById('instagramFeeds');

  if (feedsContainer) {
    INSTAGRAM_CONFIG.accounts.forEach(function (account) {
      var feedEl = createFeedElement(account);
      feedsContainer.appendChild(feedEl);

      if (account.accessToken) {
        loadInstagramPosts(account, feedEl);
      } else if (account.featuredPosts && account.featuredPosts.length > 0) {
        loadManualEmbeds(account, feedEl);
      }
      // Otherwise placeholder is already showing
    });
  }

  function createFeedElement(account) {
    var feed = document.createElement('div');
    feed.className = 'instagram-feed';

    var initial = account.username.charAt(0).toUpperCase();

    feed.innerHTML =
      '<div class="instagram-feed-header">' +
        '<div class="instagram-feed-profile">' +
          '<div class="instagram-avatar">' + initial + '</div>' +
          '<div>' +
            '<div class="instagram-username">@' + escapeHtml(account.username) + '</div>' +
            '<div class="instagram-label">' + escapeHtml(account.label) + '</div>' +
          '</div>' +
        '</div>' +
        '<a href="' + escapeAttr(account.profileUrl) + '" target="_blank" class="instagram-follow-btn">Follow</a>' +
      '</div>' +
      '<div class="instagram-feed-content">' +
        renderPlaceholder(account) +
      '</div>' +
      '<div class="instagram-feed-footer">' +
        '<a href="' + escapeAttr(account.profileUrl) + '" target="_blank">View all on Instagram &rarr;</a>' +
      '</div>';

    return feed;
  }

  function renderPlaceholder(account) {
    var grad = account.gradient || 'linear-gradient(135deg, #374151, #4b5563)';
    var desc = account.description ? '<p class="instagram-placeholder-desc">' + escapeHtml(account.description) + '</p>' : '';

    return '<div class="instagram-placeholder">' +
      '<div class="instagram-placeholder-banner" style="background:' + grad + ';">' +
        '<svg viewBox="0 0 24 24" width="48" height="48" fill="none" stroke="rgba(255,255,255,0.6)" stroke-width="1.5"><rect x="2" y="2" width="20" height="20" rx="5"/><circle cx="12" cy="12" r="5"/><circle cx="17.5" cy="6.5" r="1.5" fill="rgba(255,255,255,0.6)" stroke="none"/></svg>' +
        '<div class="instagram-placeholder-banner-text">' +
          '<span class="instagram-placeholder-handle">@' + escapeHtml(account.username) + '</span>' +
        '</div>' +
      '</div>' +
      desc +
      '<a href="' + escapeAttr(account.profileUrl) + '" target="_blank" class="btn btn-sm instagram-cta-btn">View on Instagram &rarr;</a>' +
    '</div>';
  }

  function loadInstagramPosts(account, feedEl) {
    var fields = 'id,caption,media_type,media_url,thumbnail_url,permalink,timestamp';
    var url = 'https://graph.instagram.com/me/media?fields=' + fields +
      '&limit=' + INSTAGRAM_CONFIG.postsPerAccount +
      '&access_token=' + encodeURIComponent(account.accessToken);

    fetch(url)
      .then(function (res) { return res.json(); })
      .then(function (data) {
        if (data.data && data.data.length > 0) {
          renderInstagramGrid(data.data, feedEl);
        }
      })
      .catch(function () {
        // Token expired or error — keep placeholder
      });
  }

  function renderInstagramGrid(posts, feedEl) {
    var contentEl = feedEl.querySelector('.instagram-feed-content');

    var html = '<div class="instagram-feed-grid">';

    posts.forEach(function (post) {
      var imgUrl = post.media_type === 'VIDEO' ? (post.thumbnail_url || '') : (post.media_url || '');
      var caption = post.caption ? post.caption.substring(0, 100) : '';

      html += '<a href="' + escapeAttr(post.permalink) + '" target="_blank" class="instagram-post">' +
        '<img src="' + escapeAttr(imgUrl) + '" alt="' + escapeAttr(caption) + '" loading="lazy">' +
        '<div class="instagram-post-overlay">' +
          '<span>' + (post.media_type === 'VIDEO' ? '\u25B6' : '') + '</span>' +
        '</div>' +
      '</a>';
    });

    html += '</div>';
    contentEl.innerHTML = html;
  }

  function loadManualEmbeds(account, feedEl) {
    var contentEl = feedEl.querySelector('.instagram-feed-content');

    var html = '<div style="padding:1rem;">';
    account.featuredPosts.forEach(function (postUrl) {
      html += '<blockquote class="instagram-media" data-instgrm-permalink="' + escapeAttr(postUrl) + '" ' +
        'data-instgrm-version="14" style="max-width:100%;margin:0 0 1rem;"></blockquote>';
    });
    html += '</div>';

    contentEl.innerHTML = html;

    // Load Instagram embed script
    if (!document.querySelector('script[src*="instagram.com/embed.js"]')) {
      var script = document.createElement('script');
      script.async = true;
      script.src = 'https://www.instagram.com/embed.js';
      document.body.appendChild(script);
    } else if (window.instgrm) {
      window.instgrm.Embeds.process();
    }
  }

  // =============================================
  // HELPERS
  // =============================================
  function formatDate(dateStr) {
    var d = new Date(dateStr + 'T00:00:00');
    var months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
    return months[d.getMonth()] + ' ' + d.getFullYear();
  }

  function getCategoryIcon(category) {
    switch (category) {
      case 'Wine': return '\uD83C\uDF77';
      case 'Science': return '\uD83E\uDDEC';
      case 'Personal': return '\uD83C\uDF93';
      default: return '\uD83D\uDCF0';
    }
  }

  function escapeHtml(str) {
    var div = document.createElement('div');
    div.appendChild(document.createTextNode(str || ''));
    return div.innerHTML;
  }

  function escapeAttr(str) {
    return (str || '').replace(/&/g, '&amp;').replace(/"/g, '&quot;').replace(/'/g, '&#39;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  }
})();
