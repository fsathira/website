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
      'PhD in Biostatistics, Harvard University (2013). Advisors: John Quackenbush, Curtis Huttenhower, Dawn DeMeo, Christoph Lange',
      'MS in Computer Science, Stanford University (2009)',
      'BS (Honors) in Mathematical & Computational Science, Stanford University (2009). Advisor: Brad Efron'
    ],
    roles: [
      'Consultant, EpiCurious Innovations (2025-present) \u2014 Decade Bio: AI-driven personalized antibody-drug conjugate design; Novella Health: academic collaboration for EHR data integration',
      'Senior Director, Data Science, BigHat Biosciences (2024-2025)',
      'Director, Bioinformatics & Computational Biology, Alloy Therapeutics / 82VS Venture Studio (2021-2024). Contributed to the creation of 9 companies. Co-founded an immuno-oncology company. 3 provisional patents.',
      'Director & Head of Bioinformatics, Genomic Medicine, MedGenome (2019-2021). Led South Asian genomic medicine initiative. Designed SARGAM SNP array with Thermo Fisher.',
      'Senior Scientist, Computational Biology, 23andMe (2013-2019). Founding scientist of the Therapeutics Division. >60 co-authored publications.',
      'Co-Founder & CTO, SQ Technologies (2013-2017). Built HealthIQ App. Harvard Innovation Lab resident. Harvard Deans\u2019 Challenge 2nd Runner-up.'
    ],
    expertise: [
      'AI/ML-driven drug discovery & antibody engineering',
      'Computational biology & bioinformatics',
      'Genomics, population genetics & precision medicine',
      'Phenome-wide association studies (PheWAS)',
      'Antisense oligonucleotide (ASO) design',
      'Target discovery & evaluation'
    ],
    publications: '86 peer-reviewed papers in Nature, Science, Lancet Neurology, PLoS Computational Biology, Bioinformatics, Nucleic Acids Research, and more',
    awards: 'One of 20 LGBTQ+ Biopharma Leaders 2024 (Endpoints News). King\u2019s Scholarship, Thailand.',
    patent: 'Methods and compositions for targeting efemp1 (WO2023154964A1)',
    wine: 'Co-owner & winemaker at Sunset Cellars, Suisun Valley, CA. Certified Sommelier. AAPI wine advocacy leader.',
    fun: 'Thai-born, Stanford & Harvard educated, turns genomic data into drug targets by day and grapes into wine by weekend. Published two books. Ordained Buddhist monk. Black belt in Kyudo (Japanese Archery).'
  };

  // --- Pre-written Bios (fallback when no API configured) ---
  var FALLBACK_BIOS = {
    'oneliner-professional': 'Fah Sathirapongsasuti, PhD, is an AI and computational biology leader with 86 publications and experience spanning 23andMe (founding the Therapeutics Division), MedGenome, Alloy Therapeutics, and BigHat Biosciences, currently consulting on AI-driven antibody-drug conjugate design at Decade Bio.',
    'oneliner-casual': 'Scientist by day, winemaker by weekend \u2014 Fah turns genomic data into drug targets and grapes into award-worthy wines at Sunset Cellars.',
    'oneliner-academic': 'Dr. J. Fah Sathirapongsasuti (Harvard PhD, Stanford BS/MS) is a computational biologist specializing in integrative genomics, phenome-wide association studies, and AI-driven drug discovery, with 86 publications in Nature, Science, and Lancet Neurology.',
    'oneliner-haiku': 'Genes and grapes entwine\nAlgorithms find the cure\nData pours like wine',
    'oneliner-limerick': 'A scientist known as Fah,\nWhose code could see near and far,\nHe\'d sequence your genes,\nKnow just what it means,\nThen toast with a Suisun reserva!',
    'oneliner-eli5': 'Fah is a super-smart person who uses computers to figure out how our bodies work, and then he also makes yummy grape juice for grown-ups!',

    'short-professional': 'Dr. Fah Sathirapongsasuti is a bioinformatics and AI leader who has spent over 15 years at the forefront of computational drug discovery. With a PhD from Harvard and dual degrees from Stanford, he has held leadership roles at 23andMe (founding the Therapeutics Division), MedGenome, Alloy Therapeutics (contributing to the creation of 9 companies), and BigHat Biosciences. He is currently consulting through EpiCurious Innovations, designing AI-driven personalized antibody-drug conjugates at Decade Bio and driving academic collaboration for Novella Health. With 86 peer-reviewed publications in Nature, Science, and Lancet Neurology, and a patent in antisense oligonucleotide therapeutics, his work bridges computational genomics with real-world drug development. When not decoding genomes, he crafts wines at Sunset Cellars in Suisun Valley, California.',
    'short-casual': 'Fah Sathirapongsasuti is one of those rare people who\'s equally at home in a bioinformatics lab and a wine cellar. Armed with a Harvard PhD and Stanford degrees in CS and math, he\'s spent his career at the intersection of AI and biology \u2014 from building 23andMe\'s Therapeutics Division from scratch (86 publications later!) to leading computational teams at MedGenome, Alloy Therapeutics, and BigHat Biosciences. These days, he\'s consulting through EpiCurious Innovations, designing AI-driven antibody-drug conjugates at Decade Bio, while running Sunset Cellars, his winery in Suisun Valley. Because why choose between genomes and grapes?',
    'short-academic': 'J. Fah Sathirapongsasuti, PhD, is a computational biologist specializing in integrative genomics, phenome-wide association studies, and AI-driven drug discovery. He received his PhD in Biostatistics from Harvard University (advisors: Quackenbush, Huttenhower, DeMeo, Lange) and BS/MS degrees from Stanford University (advisor: Efron). His 86 publications include landmark studies in Nature, Science, and Lancet Neurology, and he is the developer of ExomeCNV. He has held leadership roles at 23andMe (founding scientist, Therapeutics Division), MedGenome, Alloy Therapeutics, and BigHat Biosciences. He currently consults on AI-driven antibody-drug conjugate design at Decade Bio through EpiCurious Innovations.',
    'short-haiku': 'Code reads the genome,\nMachine learning lights the path\u2014\nVineyard sunset glows.\n\nFrom Stanford\'s bright halls\nTo Harvard\'s deep statistics,\nTruth lives in the data.\n\nGrapes ferment like thoughts,\nAlgorithms and terroir blend\u2014\nScience meets the soil.',
    'short-limerick': 'A bioinformatician named Fah,\nStudied genomes both near and far,\nFrom Stanford to Harvard,\nHis datasets he harvested,\nThen made Pinot Noir from his terroir!\n\nWith PheWAS he\'d find a drug lead,\nAnd publish in Nature with speed,\nBut come Friday night,\nHe\'d barrel with might,\nBecause science and wine are his creed!',
    'short-eli5': 'Imagine you have a really, really big instruction book inside your body called DNA. Fah\'s job is to use super-powerful computers to read these instruction books and figure out how to make medicines that help people feel better. He went to two really famous schools \u2014 Stanford and Harvard \u2014 to learn how to do this. And when he\'s done working on computers, he goes to his special garden where he grows grapes and turns them into wine! He\'s been doing the science part for over 15 years and has helped write more than 80 important science papers.',

    'full-professional': 'Dr. J. Fah Sathirapongsasuti is a leader in AI-driven computational biology with over 15 years of experience applying machine learning and bioinformatics to accelerate drug discovery and precision medicine.\n\nFah earned his PhD in Biostatistics from Harvard University (2013), with advisors John Quackenbush, Curtis Huttenhower, Dawn DeMeo, and Christoph Lange, and coursework in Medical Ventures (MIT Sloan) and Entrepreneurship (Harvard Business School). He holds both a Master\'s in Computer Science and a Bachelor\'s (Honors) in Mathematical & Computational Science from Stanford University (2009), where he studied under Brad Efron.\n\nHis career has spanned multiple pioneering roles in biotechnology. At 23andMe, he was a founding scientist of the Therapeutics Division, pioneering phenome-wide association studies (PheWAS), fine-mapping methodologies, microbiome, and CNV detection initiatives. His drug target validation results helped launch the internal therapeutics division and close deals with key partners. He led academic and industry collaborations resulting in over 60 co-authored publications, including a seminal Science paper on the genetics of sexual behavior.\n\nHe subsequently served as Director and Head of Bioinformatics at MedGenome\'s Genomic Medicine Division, where he led the South Asian genomic medicine initiative and collaborated with Thermo Fisher to design the first South Asian-specific SNP array "SARGAM." At Alloy Therapeutics / 82VS Venture Studio, he led Target Discovery and Evaluation, contributing to the creation of nine companies and co-founding an immuno-oncology company. He authored three provisional patent applications and designed antisense oligonucleotides now in pre-clinical development. At BigHat Biosciences, as Senior Director of Data Science, he led a seven-person cross-disciplinary team building platform capabilities across 5+ active pipelines.\n\nHe currently consults through EpiCurious Innovations, designing AI-driven personalized antibody-drug conjugates at Decade Bio and driving academic collaboration for Novella Health.\n\nFah has authored 86 peer-reviewed publications in leading journals including Nature, Science, Lancet Neurology, PLoS Computational Biology, and Bioinformatics. He was recognized as one of 20 LGBTQ+ Biopharma Leaders in 2024 by Endpoints News and holds a patent for antisense oligonucleotide therapeutics targeting EFEMP1.\n\nBeyond science, Fah is co-owner and winemaker at Sunset Cellars in Suisun Valley, California. A Certified Sommelier, he is an active advocate for AAPI representation in the wine industry.',
    'full-casual': 'Meet Fah Sathirapongsasuti \u2014 a computational biologist turned AI drug discovery leader turned winemaker, which is exactly the kind of career trajectory that happens when you\'re endlessly curious and refuse to pick just one passion.\n\nFah grew up in Thailand on a King\'s Scholarship, then headed to Stanford where he picked up degrees in both Mathematical & Computational Science and Computer Science (because one wasn\'t enough, and studying under Brad Efron doesn\'t hurt). From there, he went to Harvard for a PhD in Biostatistics, because if you\'re going to decode the human genome, you\'d better know your statistics.\n\nHis career reads like a tour of biotech\'s greatest hits. At 23andMe, he was a founding scientist of their Therapeutics Division, figuring out how to use data from millions of spit kits to actually develop drugs. He proved that consumer genomics data could predict how well drugs would work, co-authored 60+ publications including a seminal Science paper on the genetics of sexual behavior, and helped close deals with key pharma partners. After that, he led genomic medicine at MedGenome (designing the first South Asian-specific SNP array), ran target discovery at Alloy Therapeutics (contributing to creating 9 companies!), and built data science at BigHat Biosciences.\n\nWith 86 publications in journals like Nature, Science, and Lancet Neurology, plus a patent in antisense oligonucleotide therapeutics, he\'s also co-founded a health tech startup and was named one of 20 LGBTQ+ Biopharma Leaders by Endpoints News in 2024.\n\nBut here\'s where it gets fun: Fah is also the co-owner and winemaker at Sunset Cellars in Suisun Valley, California. A Certified Sommelier, he\'s become a vocal advocate for AAPI representation in the wine world.\n\nCurrently consulting through EpiCurious Innovations \u2014 designing AI-driven antibody-drug conjugates at Decade Bio \u2014 Fah continues to find new ways to make AI work for human health, and continues to prove that the best scientists also know their way around a barrel room.',
    'full-academic': 'J. Fah Sathirapongsasuti, PhD\n\nDr. Sathirapongsasuti is a computational biologist with extensive expertise in integrative genomics, statistical genetics, and the application of machine learning to therapeutic development. He received his Doctor of Philosophy in Biostatistics from Harvard University (2013; advisors: Quackenbush, Huttenhower, DeMeo, Lange) and his Master of Science in Computer Science and Bachelor of Science (Honors) in Mathematical & Computational Science from Stanford University (2009; advisor: Efron).\n\nDr. Sathirapongsasuti\'s research career began at Stanford University, where he worked in the Bejerano Lab developing statistical tools for cis-regulatory element analysis\u2014work that later evolved into GREAT (Genomic Regions Enrichment of Annotations Tool, Nature Biotechnology, 2010). At Harvard, his doctoral work focused on eQTL network analysis, microbial community ecology, and copy-number variation detection. He was co-first author on a key Human Microbiome Project publication characterizing microbial co-occurrence networks across body sites (PLoS Computational Biology, 2012; 956 citations) and developed ExomeCNV, a widely-used R package for exome sequencing-based CNV detection (Bioinformatics, 2011; 403 citations).\n\nProfessionally, Dr. Sathirapongsasuti was a founding scientist of the 23andMe Therapeutics Division, where he pioneered phenome-wide association studies (PheWAS), fine-mapping methodologies, and microbiome analysis. He was lead analyst on a seminal Science paper on the genetics of same-sex sexual behavior (2019). He subsequently led genomic medicine at MedGenome (where he co-designed the SARGAM South Asian SNP array with Thermo Fisher), directed bioinformatics at Alloy Therapeutics / 82VS Venture Studio (contributing to the creation of 9 companies), and led data science at BigHat Biosciences. He holds a patent for antisense oligonucleotide therapeutics targeting EFEMP1 (WO2023154964A1).\n\nHis publication record includes 86 peer-reviewed articles in journals including Nature, Science, Lancet Neurology, Nature Communications, PLoS Computational Biology, and Bioinformatics. He serves as a referee for Nature Methods, Nature Communications, Genome Research, and other leading journals, and served on Illumina\'s Global Screening Array design working group. He was an external examiner for the Dissertation Advisory Committee at Mahidol University (2022-2024).\n\nDr. Sathirapongsasuti currently consults through EpiCurious Innovations on AI-driven personalized antibody-drug conjugate design at Decade Bio and academic collaboration for Novella Health.',
    'full-haiku': 'I. Origin\nBangkok sun rises\nA boy dreams in binary\nThailand to Stanford\n\nII. Education\nCardinal red campus\nMath and code weave double helix\nKnowledge takes its shape\n\nCrimson leaves now fall\nHarvard\'s biostatistics\nPhD earned in full\n\nIII. Discovery\nFour million genomes\nTwenty-three and me reveals\nDrug targets emerge\n\nMicrobiome\'s web\nNature publishes the map\nBacteria speak\n\nIV. Leadership\nMedGenome\'s data\nAlloy forges antibodies\nBigHat thinks in AI\n\nV. The Vineyard\nSuisun Valley fog\nPinot noir in oak barrels\nScience meets the soil\n\nVI. Now\nNovellia calls\nAlgorithms heal the world\nSunset, pour the wine',
    'full-limerick': 'From Bangkok to Stanford he came,\nWith math and with code as his aim,\nHe mastered CS,\nGot degrees\u2014nothing less,\nThen Harvard secured his great name!\n\nAt 23andMe, he was first,\nIn Therapeutics, well-versed,\nWith PheWAS so keen,\nThe best ever seen,\nHe predicted which drugs were the worst!\n\nTo MedGenome next he would go,\nWhere genomes were studied just so,\nThen Alloy\'s great mission,\nWith BigHat\'s precision,\nMade antibodies learn, think, and grow!\n\nBut science was only half fun,\nWhen the day\'s computation was done,\nHe\'d head to his vines,\nCraft remarkable wines,\nAt Sunset Cellars, second to none!\n\nNow consulting, he bridges the gap,\nBetween AI and biology\'s map,\nWith eighty-six papers,\nAnd vineyards with capers,\nThis scientist-winemaker\'s no cap!',
    'full-eli5': 'Okay, so imagine your body is like a really, really, REALLY big LEGO set. Like, billions and billions of tiny pieces. And all those pieces have a special instruction book called DNA that tells them how to fit together.\n\nFah\'s job is to use super-smart computers to read those instruction books. He went to two of the most famous schools in the world \u2014 Stanford (that\'s in sunny California!) and Harvard (that\'s near Boston where it snows a lot!) \u2014 to learn how to do this.\n\nAfter school, he worked at a company called 23andMe. You know how some people spit into a tube to learn about where their family comes from? Well, Fah helped use ALL that information \u2014 from MILLIONS of people! \u2014 to figure out which medicines might work best for different people. Pretty cool, right?\n\nThen he worked at a bunch of other science companies, always using computers and something called "artificial intelligence" (that\'s when we teach computers to think and learn by themselves!) to find new ways to make medicine.\n\nHe\'s written more than 80 important science papers \u2014 that\'s like writing 20 really important book reports, except scientists all over the world read them!\n\nBut here\'s the really fun part: when Fah isn\'t doing science, he makes WINE! He has his own winery called Sunset Cellars where he grows grapes and turns them into wine. And he and his husband make wines that celebrate their family\'s Japanese and Thai heritage.\n\nSo basically, Fah is a super-smart scientist who teaches computers to find new medicines, AND he makes wine. He\'s like a real-life superhero, except his superpower is reading DNA and making really good grape juice for grown-ups!'
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
