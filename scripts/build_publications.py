#!/usr/bin/env python3
"""Merge PubMed and Google Scholar CSVs into publications.json."""

import csv
import json
import re
import sys
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
PUBMED_CSV = DATA_DIR / "pubmed-20260218.csv"
GSCHOLAR_CSV = DATA_DIR / "googlescholar-20260218.csv"
EXISTING_JSON = DATA_DIR / "publications.json"
OUTPUT_JSON = DATA_DIR / "publications.json"

AUTHOR_NAME_PATTERNS = [
    "Sathirapongsasuti",
]


def is_correction(title):
    """Check if a paper is a correction/erratum."""
    lower = title.lower()
    return any(w in lower for w in [
        "author correction",
        "publisher correction",
        "correction to:",
        "correction:",
        "corrigendum",
        "erratum",
    ])


def normalize_title(title):
    """Normalize title for comparison (handle en-dash, quotes, etc.)."""
    t = title.lower()
    # Normalize dashes
    t = t.replace('\u2013', '-').replace('\u2014', '-').replace('\u2012', '-')
    # Normalize quotes
    t = t.replace('\u201c', '"').replace('\u201d', '"')
    t = t.replace('\u2018', "'").replace('\u2019', "'")
    # Remove extra whitespace
    t = re.sub(r'\s+', ' ', t).strip()
    return t


def author_is_named(authors_str):
    """Check if Sathirapongsasuti appears directly in the author list."""
    return "Sathirapongsasuti" in authors_str


def determine_role(authors_str, title):
    """Determine the author's role in the paper."""
    if not author_is_named(authors_str):
        return "consortium"

    # Check if first author
    authors_lower = authors_str.lower()
    if authors_lower.startswith("sathirapongsasuti"):
        return "first_author"

    # Check if co-first (indicated by * in Google Scholar data or known)
    # For the microbial co-occurrence paper
    if "co-occurrence" in title.lower() and "Sathirapongsasuti" in authors_str:
        return "co_first_author"

    return "named_author"


def slugify(text, year):
    """Create a URL-friendly slug from title and year."""
    # Take first few meaningful words
    words = re.sub(r'[^a-z0-9\s]', '', text.lower()).split()
    # Remove common words
    stop = {'a', 'an', 'the', 'of', 'in', 'for', 'and', 'or', 'to', 'with', 'from', 'by', 'on', 'is', 'are', 'was', 'were', 'that', 'this', 'its'}
    words = [w for w in words if w not in stop][:5]
    slug = '-'.join(words)
    return f"{slug}-{year}"


def abbreviate_authors_pubmed(authors_str, max_authors=6):
    """Abbreviate a PubMed-style author list."""
    # PubMed uses: "LastName1 AB, LastName2 CD, ...; ConsortiumName; ..."
    # Split on semicolons to separate consortia
    parts = authors_str.split(";")
    main_authors = parts[0].strip() if parts else authors_str

    # Split individual authors by comma
    author_list = [a.strip() for a in main_authors.split(",") if a.strip()]

    # Highlight Sathirapongsasuti
    result = []
    for a in author_list[:max_authors]:
        if "Sathirapongsasuti" in a:
            result.append(a)
        else:
            result.append(a)

    if len(author_list) > max_authors:
        # Check if Sathirapongsasuti is in the truncated part
        remaining = author_list[max_authors:]
        has_sathi = any("Sathirapongsasuti" in a for a in remaining)
        if has_sathi:
            sathi = next(a for a in remaining if "Sathirapongsasuti" in a)
            result.append("...")
            result.append(sathi)
            result.append("et al.")
        else:
            result.append("et al.")

    # Add consortium info
    consortia = [p.strip() for p in parts[1:] if p.strip() and not p.strip().endswith(".")]
    if consortia:
        consortium_str = "; ".join(consortia[:2])
        if consortium_str:
            result.append(f"({consortium_str})")

    return ", ".join(result)


def abbreviate_authors_gscholar(authors_str):
    """Abbreviate a Google Scholar-style author list."""
    if not authors_str:
        return ""
    # GScholar uses: "Last, First; Last, First; ..."
    authors = [a.strip().rstrip(", ") for a in authors_str.split(";") if a.strip()]

    result = []
    for a in authors[:6]:
        if "Sathirapongsasuti" in a:
            result.append(a.strip())
        else:
            result.append(a.strip())

    if len(authors) > 6:
        remaining = authors[6:]
        has_sathi = any("Sathirapongsasuti" in a for a in remaining)
        if has_sathi:
            sathi = next(a for a in remaining if "Sathirapongsasuti" in a)
            result.append("...")
            result.append(sathi.strip())
            result.append("et al.")
        else:
            result.append("et al.")

    return ", ".join(result)


def parse_pubmed_csv(filepath):
    """Parse the PubMed CSV export."""
    pubs = {}
    with open(filepath, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            doi = row.get("DOI", "").strip()
            pmid = row.get("PMID", "").strip()
            title = row.get("Title", "").strip()

            if is_correction(title):
                continue

            # Skip known correction DOIs that have same title as originals
            skip_dois = {
                "10.1016/j.neuron.2018.08.029",  # Correction to migraine paper
                "10.1172/JCI143863",  # Republication/correction of HDAC6 COPD paper
            }
            if doi in skip_dois:
                continue

            authors_raw = row.get("Authors", "").strip()
            role = determine_role(authors_raw, title)
            authors = abbreviate_authors_pubmed(authors_raw)

            journal = row.get("Journal/Book", "").strip()
            year = int(row.get("Publication Year", "0") or "0")
            citation = row.get("Citation", "").strip()
            pmcid = row.get("PMCID", "").strip()

            # Extract volume info from citation
            volume = ""
            vol_match = re.search(r'(\d{4})\s*;?\s*(\d+.*?)\.?\s*doi:', citation) if citation else None
            if vol_match:
                volume = vol_match.group(2).strip().rstrip(".")

            key = doi if doi else f"pmid-{pmid}"
            pubs[key] = {
                "title": title,
                "authors": authors,
                "authors_raw": authors_raw,
                "journal": journal,
                "year": year,
                "volume": volume,
                "doi": doi,
                "pubmed_id": pmid,
                "pmcid": pmcid,
                "role": role,
            }
    return pubs


def parse_gscholar_csv(filepath):
    """Parse the Google Scholar CSV export."""
    pubs = {}
    with open(filepath, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            title = row.get("Title", "").strip()
            if not title or is_correction(title):
                continue

            authors_raw = row.get("Authors", "").strip()
            authors = abbreviate_authors_gscholar(authors_raw)

            journal = row.get("Publication", "").strip()
            year = int(row.get("Year", "0") or "0")
            volume_parts = []
            if row.get("Volume"):
                volume_parts.append(row["Volume"].strip())
            if row.get("Number"):
                volume_parts.append(f"({row['Number'].strip()})")
            if row.get("Pages"):
                volume_parts.append(f":{row['Pages'].strip()}")
            volume = "".join(volume_parts)

            role = determine_role(authors_raw, title)

            # Use title as key since GScholar doesn't have DOI
            pubs[title.lower()] = {
                "title": title,
                "authors": authors,
                "authors_raw": authors_raw,
                "journal": journal,
                "year": year,
                "volume": volume,
                "role": role,
                "publisher": row.get("Publisher", "").strip(),
            }
    return pubs


def load_existing_json(filepath):
    """Load existing publications.json."""
    with open(filepath, "r", encoding="utf-8") as f:
        return json.load(f)


def auto_tags(title, journal, role):
    """Generate automatic tags based on title and journal."""
    tags = []
    lower_title = title.lower()
    lower_journal = journal.lower()

    if "gwas" in lower_title or "genome-wide association" in lower_title:
        tags.append("GWAS")
    if "microbiome" in lower_title or "microbial" in lower_title:
        tags.append("microbiome")
    if "parkinson" in lower_title:
        tags.append("neurogenetics")
    if "cancer" in lower_title or "melanoma" in lower_title or "oncol" in lower_journal:
        tags.append("cancer")
    if "psychiatr" in lower_title or "depress" in lower_title or "adhd" in lower_title or "attention" in lower_title:
        tags.append("psychiatric-genetics")
    if "covid" in lower_title or "sars" in lower_title or "ace2" in lower_title:
        tags.append("COVID-19")
    if "copy-number" in lower_title or "cnv" in lower_title or "exome" in lower_title:
        tags.append("bioinformatics")
    if "insomnia" in lower_title or "sleep" in lower_title or "nap" in lower_title:
        tags.append("sleep-genetics")
    if "prion" in lower_title:
        tags.append("neurogenetics")

    if role == "consortium":
        tags.append("consortium")
    if "23andme" in lower_title.lower() or role in ("named_author", "first_author", "co_first_author"):
        if "human-genetics" not in tags:
            tags.append("human-genetics")

    if not tags:
        tags.append("human-genetics")

    return tags


def merge_publications(pubmed, gscholar, existing):
    """Merge all publication sources."""
    # Index existing by DOI and PMID
    existing_by_doi = {}
    existing_by_pmid = {}
    existing_custom = []  # entries without DOI/PMID (conference, unpublished)

    for entry in existing:
        doi = entry.get("doi", "")
        pmid = entry.get("pubmed_id", "")
        if doi:
            existing_by_doi[doi] = entry
        if pmid:
            existing_by_pmid[pmid] = entry
        if not doi and not pmid:
            existing_custom.append(entry)

    # Build final list
    merged = {}

    # 1. Add all PubMed entries
    for key, pub in pubmed.items():
        doi = pub["doi"]
        pmid = pub["pubmed_id"]

        # Check if we have existing data to merge
        existing_entry = None
        if doi and doi in existing_by_doi:
            existing_entry = existing_by_doi[doi]
        elif pmid and pmid in existing_by_pmid:
            existing_entry = existing_by_pmid[pmid]

        entry = {
            "title": pub["title"],
            "authors": pub["authors"],
            "journal": pub["journal"],
            "year": pub["year"],
            "volume": pub.get("volume", ""),
            "doi": doi,
            "pubmed_id": pmid,
            "pmcid": pub.get("pmcid", ""),
            "role": pub["role"],
        }

        if existing_entry:
            # Preserve custom abstract and tags
            entry["abstract"] = existing_entry.get("abstract", "")
            entry["tags"] = existing_entry.get("tags", auto_tags(pub["title"], pub["journal"], pub["role"]))
            entry["id"] = existing_entry.get("id", slugify(pub["title"], pub["year"]))
            # Keep existing authors format if it has useful info (like co-first notation)
            if "co-first" in existing_entry.get("authors", "").lower() or "*" in existing_entry.get("authors", ""):
                entry["authors"] = existing_entry["authors"]
        else:
            entry["abstract"] = ""
            entry["tags"] = auto_tags(pub["title"], pub["journal"], pub["role"])
            entry["id"] = slugify(pub["title"], pub["year"])

        merged[doi if doi else f"pmid-{pmid}"] = entry

    # 2. Add Google Scholar entries not already in PubMed (by title matching)
    for key, pub in gscholar.items():
        norm_title = normalize_title(pub["title"])
        # Check if already in merged (by normalized title similarity)
        already_exists = False
        for existing_key, existing_pub in merged.items():
            if norm_title == normalize_title(existing_pub["title"]):
                already_exists = True
                # Update with Google Scholar info if PubMed was missing volume
                if not existing_pub.get("volume") and pub.get("volume"):
                    existing_pub["volume"] = pub["volume"]
                break

        if not already_exists:
            entry = {
                "id": slugify(pub["title"], pub["year"]),
                "title": pub["title"],
                "authors": pub["authors"],
                "journal": pub["journal"],
                "year": pub["year"],
                "volume": pub.get("volume", ""),
                "doi": "",
                "pubmed_id": "",
                "pmcid": "",
                "role": pub["role"],
                "abstract": "",
                "tags": auto_tags(pub["title"], pub["journal"], pub["role"]),
            }
            merged[norm_title] = entry

    # 3. Add custom entries (conference posters, unpublished work)
    for entry in existing_custom:
        entry_id = entry.get("id", "custom-" + str(len(merged)))
        if "role" not in entry:
            entry["role"] = "first_author"
        if "pmcid" not in entry:
            entry["pmcid"] = ""
        merged[entry_id] = entry

    return merged


def format_output(merged):
    """Format the merged publications as a sorted list."""
    pubs = list(merged.values())

    # Sort: by year (desc), then by role priority, then by title
    role_priority = {
        "first_author": 0,
        "co_first_author": 1,
        "named_author": 2,
        "consortium": 3,
        "other": 4,
    }

    pubs.sort(key=lambda p: (
        -p.get("year", 0),
        role_priority.get(p.get("role", "other"), 4),
        p.get("title", ""),
    ))

    # Clean up: ensure all fields present
    output = []
    for p in pubs:
        entry = {
            "id": p.get("id", ""),
            "title": p.get("title", ""),
            "authors": p.get("authors", ""),
            "journal": p.get("journal", ""),
            "year": p.get("year", 0),
            "volume": p.get("volume", ""),
            "doi": p.get("doi", ""),
            "pubmed_id": p.get("pubmed_id", ""),
            "pmcid": p.get("pmcid", ""),
            "role": p.get("role", "other"),
            "abstract": p.get("abstract", ""),
            "tags": p.get("tags", []),
        }
        output.append(entry)

    return output


def main():
    pubmed = parse_pubmed_csv(PUBMED_CSV)
    print(f"PubMed entries (excluding corrections): {len(pubmed)}")

    gscholar = parse_gscholar_csv(GSCHOLAR_CSV)
    print(f"Google Scholar entries: {len(gscholar)}")

    existing = load_existing_json(EXISTING_JSON)
    print(f"Existing JSON entries: {len(existing)}")

    merged = merge_publications(pubmed, gscholar, existing)
    print(f"Merged total: {len(merged)}")

    output = format_output(merged)

    # Count by role
    roles = {}
    for p in output:
        role = p.get("role", "other")
        roles[role] = roles.get(role, 0) + 1
    print(f"By role: {roles}")

    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    print(f"Written {len(output)} publications to {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
