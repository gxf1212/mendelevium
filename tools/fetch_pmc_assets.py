#!/usr/bin/env python3
"""
Fetch open-access fulltext assets from NCBI PubMed Central (PMC).

What this tool does well:
- Given a PMCID (e.g., "PMC1234567") or a PMC article URL, finds the best PDF link
  and *attempts* to download it.
- Also collects supplementary file links listed on the PMC page.
- Writes a small metadata JSON alongside downloads for traceability/reuse.

Important limitation (as of 2026-01):
- Some PMC downloads are protected by an anti-bot Proof-of-Work (PoW) page.
  In that case, direct viewer links may fail. This tool will try the official
  NCBI OA service (often provides FTP links) as a fallback; if still blocked,
  it will tell you to download via a normal browser.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Optional


PMC_BASE = "https://pmc.ncbi.nlm.nih.gov"
OA_FCGI = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi"


@dataclass
class FetchResult:
    input: str
    resolved_pmcid: Optional[str]
    pmc_url: Optional[str]
    pdf_url: Optional[str]
    pdf_path: Optional[str]
    supplementary_urls: list[str]
    blocked_by_pow: bool
    notes: list[str]


def _http_get(url: str, *, timeout_s: int = 40) -> bytes:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) fetch_pmc_assets/1.0",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
        },
    )
    with urllib.request.urlopen(req, timeout=timeout_s) as resp:
        return resp.read()


def _http_download(url: str, out_path: Path, *, timeout_s: int = 120) -> tuple[bool, str]:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) fetch_pmc_assets/1.0",
            "Accept": "*/*",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            content_type = (resp.headers.get("content-type") or "").lower()
            data = resp.read()
    except urllib.error.HTTPError as e:
        return False, f"HTTPError {e.code}: {e.reason}"
    except Exception as e:  # noqa: BLE001
        return False, f"{type(e).__name__}: {e}"

    # Detect PMC PoW / anti-bot HTML.
    if "text/html" in content_type:
        text = data.decode("utf-8", "ignore")
        if _looks_like_pmc_pow(text):
            return False, "Blocked by PMC anti-bot PoW (HTML challenge page)."

    out_path.write_bytes(data)
    return True, "ok"


def _looks_like_pmc_pow(html: str) -> bool:
    # Observed patterns on PMC when it blocks automated downloads:
    # - Title: "Preparing to download ..."
    # - JS includes a PoW initializer and cookie name "cloudpmc-viewer-pow"
    h = html.lower()
    return (
        "preparing to download" in h
        or "cloudpmc-viewer-pow" in h
        or "pow_challenge" in h
        or "/static/assets/pow-" in h
    )


def _normalize_pmcid(s: str) -> Optional[str]:
    s = s.strip()
    m = re.search(r"\bPMC\d+\b", s, flags=re.IGNORECASE)
    if m:
        return m.group(0).upper()
    return None


def _resolve_pmcid_from_pmid(pmid: str) -> Optional[str]:
    # Use NCBI eutils to map PMID -> PMCID, if available.
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?"
        + urllib.parse.urlencode(
            {
                "dbfrom": "pubmed",
                "db": "pmc",
                "linkname": "pubmed_pmc",
                "id": pmid,
            }
        )
    )
    try:
        xml_bytes = _http_get(url, timeout_s=40)
    except Exception:
        return None

    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError:
        return None

    # elink returns numeric PMC IDs (e.g., 12519464). We can convert to PMCID
    # by prefixing "PMC" if we can find the canonical PMCID elsewhere.
    # Unfortunately, PMCIDs are not always 1:1 with PMCID strings in this XML.
    # So here we only return None and ask user to supply PMCID/URL directly if needed.
    # (Kept for future expansion; avoids silently returning a wrong ID.)
    _ = root
    return None


def _pick_pdf_url(pmc_url: str, pmc_html: str) -> Optional[str]:
    # Prefer relative "pdf/<name>.pdf" links on the PMC article page.
    # Example: href="pdf/id5c00517.pdf"
    candidates = re.findall(r'href="(pdf/[^"]+?\.pdf)"', pmc_html, flags=re.IGNORECASE)
    if candidates:
        return urllib.parse.urljoin(pmc_url, candidates[0])

    # Fallback: any .pdf under pmc domain
    any_pdf = re.findall(r'href="([^"]+?\.pdf)"', pmc_html, flags=re.IGNORECASE)
    for href in any_pdf:
        abs_url = urllib.parse.urljoin(pmc_url, href)
        if abs_url.startswith(PMC_BASE) and "/pdf/" in abs_url:
            return abs_url
    return None


def _collect_supp_urls(pmc_html: str) -> list[str]:
    # Common supplementary patterns on PMC:
    # - /articles/instance/<num>/bin/<file>
    # - /articles/PMCxxxxxxx/bin/<file>
    urls: set[str] = set()
    for href in re.findall(r'href="([^"]+)"', pmc_html, flags=re.IGNORECASE):
        if not href:
            continue
        if "/bin/" not in href:
            continue
        abs_url = urllib.parse.urljoin(PMC_BASE + "/", href)
        if abs_url.startswith(PMC_BASE):
            urls.add(abs_url)
    return sorted(urls)


def _get_oa_links(pmcid: str) -> list[tuple[str, str]]:
    """
    Return (format, href) from NCBI OA service.

    This often provides FTP links that are not behind the PMC viewer PoW wall,
    e.g.:
      - ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_pdf/.../main.PMCxxxx.pdf
      - ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_package/.../PMCxxxx.tar.gz
    """
    url = OA_FCGI + "?" + urllib.parse.urlencode({"id": pmcid})
    try:
        xml_bytes = _http_get(url, timeout_s=40)
    except Exception:
        return []
    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError:
        return []

    out: list[tuple[str, str]] = []
    for link in root.findall(".//link"):
        fmt = link.attrib.get("format")
        href = link.attrib.get("href")
        if fmt and href:
            out.append((fmt, href))
    return out


def _extract_pdf_from_tgz(tgz_path: Path, out_dir: Path) -> tuple[Optional[Path], list[Path], list[str]]:
    """
    Extract a PMC OA package .tar.gz and return:
      - best_pdf (largest .pdf), extracted_pdfs, notes
    """
    notes: list[str] = []
    out_dir.mkdir(parents=True, exist_ok=True)
    pdfs: list[Path] = []

    try:
        with tarfile.open(tgz_path, "r:gz") as tf:
            members = tf.getmembers()
            tf.extractall(out_dir)
    except Exception as e:  # noqa: BLE001
        return None, [], [f"Failed to extract OA package: {type(e).__name__}: {e}"]

    for m in members:
        if not m.name.lower().endswith(".pdf"):
            continue
        candidate = out_dir / m.name
        if candidate.exists():
            pdfs.append(candidate)

    if not pdfs:
        notes.append("OA package extracted, but no PDF found inside.")
        return None, [], notes

    best = max(pdfs, key=lambda p: p.stat().st_size if p.exists() else 0)
    if len(pdfs) > 1:
        notes.append(f"Multiple PDFs found in OA package; picked largest: {best.name}")
    return best, pdfs, notes


def _safe_slug(s: str) -> str:
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s.strip())
    return s[:200] if s else "pmc"


def fetch_one(inp: str, out_dir: Path, *, sleep_s: float = 0.0) -> FetchResult:
    notes: list[str] = []
    blocked = False

    pmcid = _normalize_pmcid(inp)
    if pmcid is None and inp.strip().isdigit():
        maybe = _resolve_pmcid_from_pmid(inp.strip())
        if maybe:
            pmcid = maybe
        else:
            notes.append("Input looks like a PMID, but PMID->PMCID auto-resolve is not implemented reliably; please pass a PMCID or PMC URL.")

    if pmcid is None:
        # Try treat as URL and extract PMCID.
        pmcid = _normalize_pmcid(inp)

    if pmcid is None:
        return FetchResult(
            input=inp,
            resolved_pmcid=None,
            pmc_url=None,
            pdf_url=None,
            pdf_path=None,
            supplementary_urls=[],
            blocked_by_pow=False,
            notes=notes or ["Could not resolve PMCID from input."],
        )

    pmc_url = f"{PMC_BASE}/articles/{pmcid}/"

    try:
        html = _http_get(pmc_url, timeout_s=40).decode("utf-8", "ignore")
    except Exception as e:  # noqa: BLE001
        return FetchResult(
            input=inp,
            resolved_pmcid=pmcid,
            pmc_url=pmc_url,
            pdf_url=None,
            pdf_path=None,
            supplementary_urls=[],
            blocked_by_pow=False,
            notes=[f"Failed to fetch PMC page: {type(e).__name__}: {e}"],
        )

    pdf_url = _pick_pdf_url(pmc_url, html)
    supp_urls = _collect_supp_urls(html)

    pmc_out = out_dir / pmcid
    pmc_out.mkdir(parents=True, exist_ok=True)

    # Save the HTML snapshot for debugging/traceability.
    (pmc_out / "pmc_page.html").write_text(html, encoding="utf-8")

    pdf_path: Optional[str] = None
    if pdf_url:
        pdf_name = _safe_slug(urllib.parse.urlparse(pdf_url).path.split("/")[-1] or f"{pmcid}.pdf")
        target = pmc_out / pdf_name
        ok, msg = _http_download(pdf_url, target)
        if ok:
            pdf_path = str(target)
        else:
            notes.append(f"PDF download failed: {msg}")
            if "PoW" in msg:
                blocked = True
    else:
        notes.append("No obvious PDF link found on PMC page (pdf/*.pdf).")

    # If blocked by PMC PoW (or PDF missing), try NCBI OA service links (often FTP).
    oa_links = _get_oa_links(pmcid)
    oa_pdf = [href for fmt, href in oa_links if fmt == "pdf"]
    oa_tgz = [href for fmt, href in oa_links if fmt == "tgz"]

    if (blocked or pdf_path is None) and (oa_pdf or oa_tgz):
        if oa_pdf:
            oa_pdf_url = oa_pdf[0]
            pdf_name = _safe_slug(urllib.parse.urlparse(oa_pdf_url).path.split("/")[-1] or f"main.{pmcid}.pdf")
            target = pmc_out / pdf_name
            ok, msg = _http_download(oa_pdf_url, target)
            if ok:
                pdf_path = str(target)
                pdf_url = oa_pdf_url
                blocked = False
                notes.append("PDF downloaded via NCBI OA service (FTP), bypassing PMC viewer.")
            else:
                notes.append(f"OA PDF download failed: {msg}")

        # If still no PDF, try OA package tarball and extract PDFs from it.
        if pdf_path is None and oa_tgz:
            oa_tgz_url = oa_tgz[0]
            tgz_name = _safe_slug(urllib.parse.urlparse(oa_tgz_url).path.split("/")[-1] or f"{pmcid}.tar.gz")
            tgz_path = pmc_out / tgz_name
            ok, msg = _http_download(oa_tgz_url, tgz_path)
            if ok:
                extracted_dir = pmc_out / "oa_package_extracted"
                best_pdf, _pdfs, extract_notes = _extract_pdf_from_tgz(tgz_path, extracted_dir)
                notes.extend(extract_notes)
                if best_pdf:
                    pdf_path = str(best_pdf)
                    pdf_url = oa_tgz_url
                    blocked = False
                    notes.append("PDF recovered from OA package tarball (FTP).")
                else:
                    notes.append("OA package downloaded, but no PDF was recovered.")
            else:
                notes.append(f"OA package download failed: {msg}")

    # Attempt supplementary downloads (best-effort; often PoW-protected).
    downloaded = 0
    for u in supp_urls:
        if sleep_s:
            time.sleep(sleep_s)
        name = _safe_slug(urllib.parse.urlparse(u).path.split("/")[-1] or "supp")
        target = pmc_out / "supplementary" / name
        ok, msg = _http_download(u, target)
        if ok:
            downloaded += 1
        else:
            if "PoW" in msg:
                blocked = True
            notes.append(f"Supp download failed ({name}): {msg}")

    if supp_urls and downloaded == 0:
        notes.append("Supplementary links were found, but none could be downloaded automatically (likely PoW-protected).")

    return FetchResult(
        input=inp,
        resolved_pmcid=pmcid,
        pmc_url=pmc_url,
        pdf_url=pdf_url,
        pdf_path=pdf_path,
        supplementary_urls=supp_urls,
        blocked_by_pow=blocked,
        notes=notes,
    )


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Fetch PMC PDFs/supplementary assets (best-effort; detects PoW blocks).")
    parser.add_argument("inputs", nargs="+", help="PMCID (PMCxxxx), PMC URL, or PMID (limited support).")
    parser.add_argument("--out", default="downloads/pmc", help="Output directory (default: downloads/pmc).")
    parser.add_argument("--sleep", type=float, default=0.0, help="Sleep seconds between supplementary downloads (default: 0).")
    args = parser.parse_args(argv)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    results: list[FetchResult] = []
    for inp in args.inputs:
        results.append(fetch_one(inp, out_dir, sleep_s=args.sleep))

    # Write an aggregated index file.
    index_path = out_dir / "index.json"
    index_path.write_text(json.dumps([asdict(r) for r in results], ensure_ascii=False, indent=2), encoding="utf-8")

    # Print a human-readable summary.
    for r in results:
        print(f"- input: {r.input}")
        print(f"  pmcid: {r.resolved_pmcid}")
        print(f"  pmc_url: {r.pmc_url}")
        print(f"  pdf_url: {r.pdf_url}")
        print(f"  pdf_path: {r.pdf_path}")
        print(f"  supp_urls: {len(r.supplementary_urls)}")
        print(f"  blocked_by_pow: {r.blocked_by_pow}")
        if r.notes:
            for n in r.notes[:20]:
                print(f"  note: {n}")
        print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
