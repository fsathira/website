/* ============================================
   Bio Generator — AI-powered bio with pre-written fallbacks
   ============================================ */

(function () {
  'use strict';

  // --- Configuration ---
  // To enable live AI generation, set your Anthropic API endpoint and key.
  // For security, use a proxy/serverless function rather than exposing your key client-side.
  var BIO_CONFIG = {
    apiEndpoint: '', // e.g., 'https://your-proxy.example.com/api/generate-bio'
    apiKey: '',      // Only if using a direct API call (not recommended client-side)
    model: 'claude-sonnet-4-5-20250929'
  };

  // --- Bio Data (used for system prompt and fallbacks) ---
  var BIO_DATA = {
    name: 'Fah Sathirapongsasuti',
    fullName: 'J. Fah Sathirapongsasuti, PhD',
    education: [
      'PhD in Biostatistics, Harvard University (2013)',
      'MS in Computer Science, Stanford University (2009)',
      'BS (Honors) in Mathematical & Computational Science, Stanford University (2009)'
    ],
    roles: [
      'Consultant at Novellia (current)',
      'Director of Bioinformatics & Computational Biology, Alloy Therapeutics',
      'Senior Scientist / Head of Bioinformatics, BigHat Biosciences',
      'Director & Head of Bioinformatics, MedGenome Genomic Medicine Division',
      'Senior Scientist, Computational Biology, 23andMe Research (founding member, Therapeutics Division)'
    ],
    expertise: [
      'AI/ML-driven drug discovery',
      'Computational biology & bioinformatics',
      'Genomics & precision medicine',
      'Phenome-wide association studies (PheWAS)',
      'Antibody discovery & engineering'
    ],
    publications: '20+ peer-reviewed papers in Nature, PLoS Computational Biology, Bioinformatics, Nucleic Acids Research',
    wine: 'Co-owner & winemaker at Sunset Cellars, Suisun Valley, CA. AAPI wine advocacy leader.',
    fun: 'Thai-born, Stanford & Harvard educated, turns genomic data into drug targets by day and grapes into wine by weekend.'
  };

  // --- Pre-written Bios (fallback when no API configured) ---
  var FALLBACK_BIOS = {
    'oneliner-professional': 'Fah Sathirapongsasuti, PhD, is an AI and computational biology leader with experience spanning 23andMe, MedGenome, Alloy Therapeutics, and BigHat Bio, leveraging machine learning to accelerate drug discovery and therapeutic development.',
    'oneliner-casual': 'Scientist by day, winemaker by weekend \u2014 Fah turns genomic data into drug targets and grapes into award-worthy wines at Sunset Cellars.',
    'oneliner-academic': 'Dr. J. Fah Sathirapongsasuti (Harvard PhD, Stanford BS/MS) specializes in integrative genomics, phenome-wide association studies, and machine learning applications in therapeutic target identification.',
    'oneliner-haiku': 'Genes and grapes entwine\nAlgorithms find the cure\nData pours like wine',
    'oneliner-limerick': 'A scientist known as Fah,\nWhose code could see near and far,\nHe\'d sequence your genes,\nKnow just what it means,\nThen toast with a Suisun reserva!',
    'oneliner-eli5': 'Fah is a super-smart person who uses computers to figure out how our bodies work, and then he also makes yummy grape juice for grown-ups!',

    'short-professional': 'Dr. Fah Sathirapongsasuti is a seasoned bioinformatics and AI leader who has spent over 15 years at the forefront of computational drug discovery. With a PhD from Harvard and dual degrees from Stanford, he has held leadership roles at 23andMe (founding the Therapeutics division), MedGenome, Alloy Therapeutics, and BigHat Biosciences. His work integrates machine learning, genomics, and clinical data to accelerate therapeutic development. He is currently consulting at Novellia. When he\'s not decoding genomes, he crafts wines at Sunset Cellars in Suisun Valley, California, and advocates for AAPI representation in the wine industry.',
    'short-casual': 'Fah Sathirapongsasuti is one of those rare people who\'s equally at home in a bioinformatics lab and a wine cellar. Armed with a Harvard PhD and Stanford degrees in CS and math, he\'s spent his career at the intersection of AI and biology \u2014 from building 23andMe\'s Therapeutics division from scratch to leading computational teams at biotech startups. These days, he\'s consulting at Novellia while running Sunset Cellars, his winery in Suisun Valley where he and his husband craft wines that celebrate Japanese-Thai heritage. Because why choose between genomes and grapes?',
    'short-academic': 'J. Fah Sathirapongsasuti, PhD, is a computational biologist and bioinformatician specializing in integrative genomics, phenome-wide association studies, and machine learning applications in drug development. He received his PhD in Biostatistics from Harvard University and BS/MS degrees from Stanford University. His publication record includes contributions to landmark studies such as the Human Microbiome Project (Nature, 2012) and the development of ExomeCNV for copy-number variation detection. He has held principal investigator and leadership roles at 23andMe, MedGenome, Alloy Therapeutics, and BigHat Biosciences.',
    'short-haiku': 'Code reads the genome,\nMachine learning lights the path\u2014\nVineyard sunset glows.\n\nFrom Stanford\'s bright halls\nTo Harvard\'s deep statistics,\nTruth lives in the data.\n\nGrapes ferment like thoughts,\nAlgorithms and terroir blend\u2014\nScience meets the soil.',
    'short-limerick': 'A bioinformatician named Fah,\nStudied genomes both near and far,\nFrom Stanford to Harvard,\nHis datasets he harvested,\nThen made Pinot Noir from his terroir!\n\nWith PheWAS he\'d find a drug lead,\nAnd publish in Nature with speed,\nBut come Friday night,\nHe\'d barrel with might,\nBecause science and wine are his creed!',
    'short-eli5': 'Imagine you have a really, really big instruction book inside your body called DNA. Fah\'s job is to use super-powerful computers to read these instruction books and figure out how to make medicines that help people feel better. He went to two really famous schools \u2014 Stanford and Harvard \u2014 to learn how to do this. And when he\'s done working on computers, he goes to his special garden where he grows grapes and turns them into wine! He\'s been doing the science part for over 15 years and has helped write more than 20 important science papers.',

    'full-professional': 'Dr. J. Fah Sathirapongsasuti is a distinguished leader in AI-driven computational biology with over 15 years of experience applying machine learning and bioinformatics to accelerate drug discovery and precision medicine.\n\nFah earned his PhD in Biostatistics from Harvard University (2013) and holds both a Master\'s in Computer Science and a Bachelor\'s (Honors) in Mathematical & Computational Science from Stanford University (2009). His doctoral research focused on integrative genomics, developing novel statistical methods for analyzing next-generation sequencing data.\n\nHis career has spanned multiple pioneering roles in biotechnology. At 23andMe, he was a founding member of the Therapeutics Division, where he pioneered the application of phenome-wide association studies (PheWAS) in drug development \u2014 presenting the first statistical evidence that consumer genomics data could predict drug efficacy and side effects. He subsequently served as Director and Head of Bioinformatics at MedGenome\'s Genomic Medicine Division, Director of Bioinformatics and Computational Biology at Alloy Therapeutics, and Senior Scientist / Head of Bioinformatics at BigHat Biosciences, where AI-driven antibody engineering was a core focus.\n\nHe is currently a consultant at Novellia, continuing to bridge AI/ML and therapeutic development.\n\nFah has authored over 20 peer-reviewed publications in leading journals including Nature, PLoS Computational Biology, Bioinformatics, and Nucleic Acids Research. His contributions include co-first authorship on the Human Microbiome Project\'s microbial co-occurrence network paper, the development of ExomeCNV for copy-number variation detection from exome sequencing, and foundational work that evolved into GREAT (Genomic Regions Enrichment of Annotations Tool).\n\nBeyond science, Fah is co-owner and winemaker at Sunset Cellars in Suisun Valley, California, where he and his husband craft wines celebrating their Japanese-Thai heritage. He is an active advocate for AAPI representation in the wine industry, co-creating the #drinkAAPIwine campaign.',
    'full-casual': 'Meet Fah Sathirapongsasuti \u2014 a computational biologist turned AI drug discovery leader turned winemaker, which is exactly the kind of career trajectory that happens when you\'re endlessly curious and refuse to pick just one passion.\n\nFah grew up in Thailand, then headed to Stanford where he picked up degrees in both Mathematical & Computational Science and Computer Science (because one wasn\'t enough). From there, he went to Harvard for a PhD in Biostatistics, because if you\'re going to decode the human genome, you\'d better know your statistics.\n\nHis career reads like a tour of biotech\'s greatest hits. At 23andMe, he was one of the founding members of their Therapeutics Division, figuring out how to use data from millions of spit kits to actually develop drugs. He proved that consumer genomics data could predict how well drugs would work \u2014 a pretty big deal. After that, he led bioinformatics teams at MedGenome, Alloy Therapeutics, and BigHat Biosciences, always at the intersection of AI, machine learning, and biology.\n\nHe\'s published 20+ papers in journals like Nature and Bioinformatics, contributed to the Human Microbiome Project, built tools that are still used in labs around the world, and helped develop technology that went from a Stanford research project into Nature Biotechnology.\n\nBut here\'s where it gets fun: Fah is also the co-owner and winemaker at Sunset Cellars in Suisun Valley, California. He and his husband craft wines that reflect their Japanese-Thai heritage, and Fah has become a vocal advocate for AAPI representation in the wine world, co-creating the #drinkAAPIwine campaign featured in Wine Enthusiast.\n\nCurrently consulting at Novellia, Fah continues to find new ways to make AI work for human health \u2014 and continues to prove that the best scientists also know their way around a barrel room.',
    'full-academic': 'J. Fah Sathirapongsasuti, PhD\n\nDr. Sathirapongsasuti is a computational biologist and bioinformatician with extensive expertise in integrative genomics, statistical genetics, and the application of machine learning to therapeutic development. He received his Doctor of Philosophy in Biostatistics from Harvard University (2013) and his Master of Science in Computer Science and Bachelor of Science (Honors) in Mathematical & Computational Science from Stanford University (2009).\n\nDr. Sathirapongsasuti\'s research career began at Stanford University, where he worked in the Bejerano Lab developing statistical tools for cis-regulatory element analysis\u2014work that later evolved into GREAT (Genomic Regions Enrichment of Annotations Tool, Nature Biotechnology, 2010). At Harvard, his doctoral work focused on eQTL network analysis, microbial community ecology, and copy-number variation detection. He was co-first author on a key Human Microbiome Project publication characterizing microbial co-occurrence networks across body sites (PLoS Computational Biology, 2012) and developed ExomeCNV, a widely-used R package for exome sequencing-based CNV detection (Bioinformatics, 2011).\n\nProfessionally, Dr. Sathirapongsasuti was a founding member of the 23andMe Therapeutics Division, where he pioneered phenome-wide association studies (PheWAS) to predict drug efficacy and adverse effects. He subsequently held leadership positions as Director and Head of Bioinformatics at MedGenome\'s Genomic Medicine Division, Director of Bioinformatics and Computational Biology at Alloy Therapeutics, and Head of Bioinformatics at BigHat Biosciences.\n\nHis publication record includes over 20 peer-reviewed articles in journals including Nature, Nature Communications, PLoS Computational Biology, Bioinformatics, Nucleic Acids Research, and BMC Proceedings. His h-index reflects significant citation impact across the fields of human genetics, metagenomics, and computational genomics.\n\nDr. Sathirapongsasuti is currently a consultant at Novellia.',
    'full-haiku': 'I. Origin\nBangkok sun rises\nA boy dreams in binary\nThailand to Stanford\n\nII. Education\nCardinal red campus\nMath and code weave double helix\nKnowledge takes its shape\n\nCrimson leaves now fall\nHarvard\'s biostatistics\nPhD earned in full\n\nIII. Discovery\nFour million genomes\nTwenty-three and me reveals\nDrug targets emerge\n\nMicrobiome\'s web\nNature publishes the map\nBacteria speak\n\nIV. Leadership\nMedGenome\'s data\nAlloy forges antibodies\nBigHat thinks in AI\n\nV. The Vineyard\nSuisun Valley fog\nPinot noir in oak barrels\nScience meets the soil\n\nVI. Now\nNovellia calls\nAlgorithms heal the world\nSunset, pour the wine',
    'full-limerick': 'From Bangkok to Stanford he came,\nWith math and with code as his aim,\nHe mastered CS,\nGot degrees\u2014nothing less,\nThen Harvard secured his great name!\n\nAt 23andMe, he was first,\nIn Therapeutics, well-versed,\nWith PheWAS so keen,\nThe best ever seen,\nHe predicted which drugs were the worst!\n\nTo MedGenome next he would go,\nWhere genomes were studied just so,\nThen Alloy\'s great mission,\nWith BigHat\'s precision,\nMade antibodies learn, think, and grow!\n\nBut science was only half fun,\nWhen the day\'s computation was done,\nHe\'d head to his vines,\nCraft remarkable wines,\nAt Sunset Cellars, second to none!\n\nNow consulting, he bridges the gap,\nBetween AI and biology\'s map,\nWith twenty-plus papers,\nAnd vineyards with capers,\nThis scientist-winemaker\'s no cap!',
    'full-eli5': 'Okay, so imagine your body is like a really, really, REALLY big LEGO set. Like, billions and billions of tiny pieces. And all those pieces have a special instruction book called DNA that tells them how to fit together.\n\nFah\'s job is to use super-smart computers to read those instruction books. He went to two of the most famous schools in the world \u2014 Stanford (that\'s in sunny California!) and Harvard (that\'s near Boston where it snows a lot!) \u2014 to learn how to do this.\n\nAfter school, he worked at a company called 23andMe. You know how some people spit into a tube to learn about where their family comes from? Well, Fah helped use ALL that information \u2014 from MILLIONS of people! \u2014 to figure out which medicines might work best for different people. Pretty cool, right?\n\nThen he worked at a bunch of other science companies, always using computers and something called "artificial intelligence" (that\'s when we teach computers to think and learn by themselves!) to find new ways to make medicine.\n\nHe\'s written more than 20 important science papers \u2014 that\'s like writing 20 really important book reports, except scientists all over the world read them!\n\nBut here\'s the really fun part: when Fah isn\'t doing science, he makes WINE! He has his own winery called Sunset Cellars where he grows grapes and turns them into wine. And he and his husband make wines that celebrate their family\'s Japanese and Thai heritage.\n\nSo basically, Fah is a super-smart scientist who teaches computers to find new medicines, AND he makes wine. He\'s like a real-life superhero, except his superpower is reading DNA and making really good grape juice for grown-ups!'
  };

  // --- DOM Elements ---
  var lengthSelect = document.getElementById('bioLength');
  var toneSelect = document.getElementById('bioTone');
  var generateBtn = document.getElementById('bioGenerateBtn');
  var outputText = document.getElementById('bioOutputText');
  var outputActions = document.getElementById('bioOutputActions');
  var copyBtn = document.getElementById('bioCopyBtn');
  var regenerateBtn = document.getElementById('bioRegenerateBtn');
  var configLink = document.getElementById('bioConfigLink');

  if (!generateBtn) return;

  var currentBio = '';
  var isGenerating = false;

  // --- Typing Animation ---
  function typeText(text, element, speed, callback) {
    element.textContent = '';
    var i = 0;
    var cursor = document.createElement('span');
    cursor.className = 'cursor';
    element.appendChild(cursor);

    function typeChar() {
      if (i < text.length) {
        var charNode = document.createTextNode(text.charAt(i));
        element.insertBefore(charNode, cursor);
        i++;
        setTimeout(typeChar, speed);
      } else {
        setTimeout(function () {
          if (cursor.parentNode) cursor.parentNode.removeChild(cursor);
          if (callback) callback();
        }, 600);
      }
    }

    typeChar();
  }

  // --- Generate Bio ---
  function generateBio() {
    if (isGenerating) return;
    isGenerating = true;

    var length = lengthSelect.value;
    var tone = toneSelect.value;

    generateBtn.textContent = 'Generating...';
    generateBtn.disabled = true;
    outputActions.style.display = 'none';

    if (BIO_CONFIG.apiEndpoint) {
      generateWithAPI(length, tone);
    } else {
      generateFallback(length, tone);
    }
  }

  function generateFallback(length, tone) {
    var key = length + '-' + tone;
    var bio = FALLBACK_BIOS[key] || FALLBACK_BIOS['short-professional'];

    var speed = bio.length > 500 ? 8 : bio.length > 200 ? 15 : 25;

    typeText(bio, outputText, speed, function () {
      currentBio = bio;
      outputActions.style.display = 'flex';
      generateBtn.textContent = 'Generate \u2192';
      generateBtn.disabled = false;
      isGenerating = false;
    });
  }

  function generateWithAPI(length, tone) {
    var lengthInstruction = {
      'oneliner': 'Write exactly one sentence.',
      'short': 'Write a short paragraph (3-5 sentences).',
      'full': 'Write a comprehensive bio (4-6 paragraphs).'
    };

    var toneInstruction = {
      'professional': 'Use a professional, polished tone suitable for a conference program or company bio.',
      'casual': 'Use a casual, engaging, conversational tone with personality and humor.',
      'academic': 'Use a formal academic tone suitable for a journal or university profile.',
      'haiku': 'Write the bio entirely as a series of haiku (5-7-5 syllable pattern).',
      'limerick': 'Write the bio entirely as a series of limericks (AABBA rhyme scheme).',
      'eli5': 'Explain this person\'s career as if talking to a curious 5-year-old. Use simple words and fun analogies.'
    };

    var systemPrompt = 'You are writing a biographical description of ' + BIO_DATA.fullName + '. Here are the facts:\n\n' +
      'Education: ' + BIO_DATA.education.join('; ') + '\n' +
      'Career: ' + BIO_DATA.roles.join('; ') + '\n' +
      'Expertise: ' + BIO_DATA.expertise.join(', ') + '\n' +
      'Publications: ' + BIO_DATA.publications + '\n' +
      'Wine: ' + BIO_DATA.wine + '\n\n' +
      'Only use facts provided above. Do not invent details.';

    var userPrompt = (lengthInstruction[length] || '') + ' ' + (toneInstruction[tone] || '');

    fetch(BIO_CONFIG.apiEndpoint, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': BIO_CONFIG.apiKey ? ('Bearer ' + BIO_CONFIG.apiKey) : ''
      },
      body: JSON.stringify({
        model: BIO_CONFIG.model,
        system: systemPrompt,
        prompt: userPrompt,
        max_tokens: 1024
      })
    })
    .then(function (res) { return res.json(); })
    .then(function (data) {
      var bio = data.text || data.content || data.completion || '';
      if (!bio) throw new Error('Empty response');

      var speed = bio.length > 500 ? 8 : bio.length > 200 ? 15 : 25;
      typeText(bio, outputText, speed, function () {
        currentBio = bio;
        outputActions.style.display = 'flex';
        generateBtn.textContent = 'Generate \u2192';
        generateBtn.disabled = false;
        isGenerating = false;
      });
    })
    .catch(function () {
      // Fall back to pre-written bios on API error
      generateFallback(length, tone);
    });
  }

  // --- Event Listeners ---
  generateBtn.addEventListener('click', generateBio);

  if (regenerateBtn) {
    regenerateBtn.addEventListener('click', generateBio);
  }

  if (copyBtn) {
    copyBtn.addEventListener('click', function () {
      if (currentBio && navigator.clipboard) {
        navigator.clipboard.writeText(currentBio).then(function () {
          copyBtn.textContent = 'Copied!';
          setTimeout(function () { copyBtn.textContent = 'Copy'; }, 2000);
        });
      }
    });
  }

  if (configLink) {
    configLink.addEventListener('click', function (e) {
      e.preventDefault();
      var msg = 'To enable live AI generation, edit js/bio-generator.js and set:\n\n' +
        '  BIO_CONFIG.apiEndpoint = "your-api-url"\n' +
        '  BIO_CONFIG.apiKey = "your-key"\n\n' +
        'Recommended: Use a serverless proxy (Vercel/Netlify function) to keep your API key secure.\n\n' +
        'Currently using pre-written bios (which are still fun!).';
      alert(msg);
    });
  }
})();
