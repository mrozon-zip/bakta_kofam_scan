"""
gene_presence_absence.py
========================
Wyszukuje obecność określonych genów w adnotowanych genomach bakteryjnych.
Źródła danych: pliki .tsv z bakty + wyniki KofamScan (.tsv)
Wynik: macierz obecności/nieobecności (0/1) w formacie TSV

Struktura katalogów (domyślna, konfigurowalna):
    genomes/
    ├── genome1/
    │   ├── genome1.tsv          # adnotacja bakty
    │   └── genome1_kofam.tsv    # wyniki KofamScan
    ├── genome2/
    │   ├── genome2.tsv
    │   └── genome2_kofam.tsv
    └── ...

Użycie:
    python gene_presence_absence.py \
        --genes genes.tsv \
        --genomes_dir genomes/ \
        --output matrix.tsv

    # Jeśli pliki są w płaskiej strukturze (nie w podkatalogach):
    python gene_presence_absence.py \
        --genes genes.tsv \
        --bakta_dir bakta_results/ \
        --kofam_dir kofam_results/ \
        --output matrix.tsv
"""

import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Konfiguracja kolumn — dostosuj do swoich plików jeśli nazwy się różnią
# ---------------------------------------------------------------------------

# Kolumny w pliku TSV z bakty (standardowy eksport bakty)
BAKTA_COL_GENE    = "gene"        # nazwa genu (np. gyrA)
BAKTA_COL_PRODUCT = "product"     # opis produktu
BAKTA_COL_DBXREF  = "db_xrefs"   # referencje do baz (zawiera K-numery)

# Kolumny w pliku TSV KofamScan
# Standardowy output KofamScan: znacznik, nazwa sekwencji, K-numer, próg, wynik, E-value, opis
KOFAM_COL_SIGNIFICANT = 0   # indeks kolumny ze znacznikiem "*" (znaczy: przeszło próg)
KOFAM_COL_GENE_ID     = 1   # indeks: nazwa sekwencji w genomie
KOFAM_COL_KO          = 2   # indeks: K-numer (np. K00001)

# Próg dla KofamScan: czy akceptować tylko wiersze oznaczone "*" (znacznik przejścia progu)
KOFAM_SIGNIFICANT_ONLY = True


# ---------------------------------------------------------------------------
# Parsowanie pliku z szukanymi genami
# ---------------------------------------------------------------------------

def load_gene_list(path: str) -> pd.DataFrame:
    """
    Wczytuje tabelę szukanych genów.
    Oczekiwane kolumny (TSV): 'gene' i opcjonalnie 'ko' (K-numer, np. K00789)
    Jeśli nie ma kolumny 'ko', wyszukiwanie odbywa się tylko po nazwie.

    Przykład pliku genes.tsv:
        gene    ko          function
        gyrA    K02469      DNA gyrase subunit A
        cbh     K01442      choloylglycine hydrolase
        luxS    K07173      S-ribosylhomocysteine lyase
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = [c.strip().lower() for c in df.columns]

    if "gene" not in df.columns:
        sys.exit(f"[BŁĄD] Plik '{path}' musi zawierać kolumnę 'gene'.")

    df["gene"] = df["gene"].str.strip().str.lower()

    if "ko" not in df.columns:
        print("[INFO] Brak kolumny 'ko' w pliku genów — wyszukiwanie tylko po nazwach.")
        df["ko"] = pd.NA
    else:
        # Normalizacja K-numerów: usuń białe znaki, upewnij się że zaczyna się od K
        df["ko"] = df["ko"].str.strip().str.upper()
        df.loc[~df["ko"].str.match(r"^K\d{5}$", na=False), "ko"] = pd.NA

    return df[["gene", "ko"]].drop_duplicates(subset="gene").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Parsowanie pliku TSV bakty
# ---------------------------------------------------------------------------

def load_bakta_tsv(path: str) -> pd.DataFrame:
    """
    Wczytuje plik TSV z bakty. Standardowy nagłówek bakty:
    #Sequence Id | Type | Start | Stop | Strand | Locus Tag | Gene | Product | DbXrefs

    Zwraca DataFrame z kolumnami: gene_name (małe litery), ko_numbers (set), raw_product
    """
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, comment=None)
    except Exception as e:
        print(f"  [OSTRZEŻENIE] Nie można wczytać {path}: {e}")
        return pd.DataFrame(columns=["gene_name", "ko_numbers", "raw_product"])

    # Normalizacja nagłówków — bakta używa '#Sequence Id' z hashem
    df.columns = [c.strip().lstrip("#").strip().lower().replace(" ", "_") for c in df.columns]

    # Mapowanie możliwych nazw kolumn
    col_map = {
        "gene": BAKTA_COL_GENE.lower(),
        "product": BAKTA_COL_PRODUCT.lower(),
        "db_xrefs": BAKTA_COL_DBXREF.lower(),
    }

    # Znajdź rzeczywiste nazwy kolumn w pliku
    actual_gene_col    = next((c for c in df.columns if "gene" in c), None)
    actual_product_col = next((c for c in df.columns if "product" in c), None)
    actual_dbxref_col  = next((c for c in df.columns if "xref" in c or "db_" in c), None)

    records = []
    for _, row in df.iterrows():
        gene_name = ""
        if actual_gene_col and pd.notna(row.get(actual_gene_col)):
            gene_name = str(row[actual_gene_col]).strip().lower()

        product = ""
        if actual_product_col and pd.notna(row.get(actual_product_col)):
            product = str(row[actual_product_col]).strip().lower()

        # Wyciągnij K-numery z db_xrefs (np. "KEGG:K02469, COG:COG0188")
        ko_numbers = set()
        if actual_dbxref_col and pd.notna(row.get(actual_dbxref_col)):
            xrefs = str(row[actual_dbxref_col])
            found = re.findall(r"K\d{5}", xrefs, re.IGNORECASE)
            ko_numbers = {k.upper() for k in found}

        if gene_name or ko_numbers:
            records.append({
                "gene_name": gene_name,
                "ko_numbers": ko_numbers,
                "raw_product": product,
            })

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Parsowanie wyników KofamScan
# ---------------------------------------------------------------------------

def load_kofamscan_tsv(path: str) -> set:
    """
    Wczytuje wyniki KofamScan i zwraca zbiór K-numerów obecnych w genomie.
    Standardowy format KofamScan (z nagłówkiem zaczynającym się od '#'):

        # * means the K-number assigned to the gene with highest score
        *  seq_id  K00001  threshold  score  evalue  description

    Jeśli KOFAM_SIGNIFICANT_ONLY=True, bierze tylko wiersze z '*' w pierwszej kolumnie.
    """
    ko_numbers = set()

    try:
        with open(path, "r") as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("##") or not line.strip():
                    continue
                # Linie z wynikami mogą zaczynać się od '*' lub ' '
                if line.startswith("#"):
                    continue

                parts = line.split()
                if len(parts) < 3:
                    continue

                significant = parts[KOFAM_COL_SIGNIFICANT]
                ko = parts[KOFAM_COL_KO]

                if not re.match(r"^K\d{5}$", ko):
                    continue

                if KOFAM_SIGNIFICANT_ONLY and significant != "*":
                    continue

                ko_numbers.add(ko.upper())

    except Exception as e:
        print(f"  [OSTRZEŻENIE] Nie można wczytać KofamScan {path}: {e}")

    return ko_numbers


# ---------------------------------------------------------------------------
# Wyszukiwanie genów w jednym genomie
# ---------------------------------------------------------------------------

def search_genes_in_genome(
    gene_list: pd.DataFrame,
    bakta_df: pd.DataFrame,
    kofam_ko_set: set,
) -> dict:
    """
    Dla jednego genomu sprawdza obecność każdego szukanego genu.
    Dopasowanie uznawane za sukces jeśli:
      (A) nazwa genu pasuje do nazwy w adnotacji bakty (dokładne lub jako prefiks)
      LUB
      (B) K-numer genu jest w adnotacji bakty (db_xrefs)
      LUB
      (C) K-numer genu jest w wynikach KofamScan

    Zwraca słownik: {nazwa_genu: 0 lub 1}
    """
    results = {}

    # Zbuduj pomocnicze struktury dla szybkiego wyszukiwania
    bakta_gene_names = set(bakta_df["gene_name"].dropna())
    bakta_ko_all = set()
    for ko_set in bakta_df["ko_numbers"]:
        bakta_ko_all.update(ko_set)

    all_ko = bakta_ko_all | kofam_ko_set  # suma K-numerów z obu źródeł

    for _, row in gene_list.iterrows():
        query_gene = row["gene"].lower()  # np. "gyra"
        query_ko   = row["ko"]            # np. "K02469" lub pd.NA

        found = False

        # --- Kanał 1: dopasowanie nazwy genu ---
        # Szukamy: dokładne dopasowanie LUB query jest prefiksem nazwy w genomie
        # Przykład: "clpp" znajdzie "clpp", ale też "clppa" jeśli taka forma istnieje
        if query_gene:
            for bname in bakta_gene_names:
                if bname == query_gene or bname.startswith(query_gene):
                    found = True
                    break

        # --- Kanał 2: dopasowanie K-numeru ---
        if not found and pd.notna(query_ko) and str(query_ko) in all_ko:
            found = True

        results[row["gene"]] = 1 if found else 0

    return results


# ---------------------------------------------------------------------------
# Odkrywanie plików genomów
# ---------------------------------------------------------------------------

def find_genome_files(
    genomes_dir: str | None,
    bakta_dir: str | None,
    kofam_dir: str | None,
) -> list[dict]:
    """
    Zwraca listę słowników z polami:
      - genome_name: nazwa genomu (do macierzy)
      - bakta_path:  ścieżka do pliku TSV bakty (lub None)
      - kofam_path:  ścieżka do pliku TSV KofamScan (lub None)

    Obsługuje dwa tryby:
      1. genomes_dir: struktura podkatalogów (każdy podkatalog = jeden genom)
      2. bakta_dir + kofam_dir: płaska struktura, pliki nazwane wg genomu
    """
    genomes = []

    if genomes_dir:
        base = Path(genomes_dir)
        for genome_dir in sorted(base.iterdir()):
            if not genome_dir.is_dir():
                continue
            name = genome_dir.name

            # Szukaj pliku TSV bakty (dowolna nazwa z .tsv, nie zawierająca 'kofam')
            bakta_files = [
                f for f in genome_dir.glob("*.tsv")
                if "bakta" in f.name.lower()
            ]
            kofam_files = [
                f for f in genome_dir.glob("*.tsv")
                if "kofam" in f.name.lower()
            ]
            # Alternatywnie KofamScan może mieć rozszerzenie .txt
            if not kofam_files:
                kofam_files = list(genome_dir.glob("*kofam*.txt"))

            genomes.append({
                "genome_name": name,
                "bakta_path":  str(bakta_files[0]) if bakta_files else None,
                "kofam_path":  str(kofam_files[0]) if kofam_files else None,
            })

    elif bakta_dir and kofam_dir:
        bakta_base = Path(bakta_dir)
        kofam_base = Path(kofam_dir)

        bakta_files = {f.stem: f for f in bakta_base.glob("*.tsv")}
        kofam_files = {f.stem: f for f in kofam_base.glob("*.tsv")}
        kofam_files.update({f.stem: f for f in kofam_base.glob("*.txt")})

        # Zbierz nazwy genomów z obu katalogów
        all_names = sorted(set(bakta_files.keys()) | set(kofam_files.keys()))
        for name in all_names:
            genomes.append({
                "genome_name": name,
                "bakta_path":  str(bakta_files[name]) if name in bakta_files else None,
                "kofam_path":  str(kofam_files[name]) if name in kofam_files else None,
            })

    return genomes


# ---------------------------------------------------------------------------
# Główna funkcja
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Macierz obecności/nieobecności genów w genomach bakteryjnych.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--genes", required=True,
        help="Plik TSV z listą szukanych genów. Wymagana kolumna: 'gene'. "
             "Opcjonalna: 'ko' (K-numer KEGG)."
    )
    parser.add_argument(
        "--genomes_dir", default=None,
        help="Katalog z podkatalogami dla każdego genomu."
    )
    parser.add_argument(
        "--bakta_dir", default=None,
        help="Katalog z plikami TSV bakty (tryb płaski)."
    )
    parser.add_argument(
        "--kofam_dir", default=None,
        help="Katalog z plikami TSV/TXT KofamScan (tryb płaski)."
    )
    parser.add_argument(
        "--output", default="presence_absence_matrix.tsv",
        help="Ścieżka do pliku wynikowego (TSV). Domyślnie: presence_absence_matrix.tsv"
    )
    parser.add_argument(
        "--include_details", action="store_true",
        help="Jeśli podane, obok macierzy 0/1 zapisz też plik z detalami dopasowań."
    )
    parser.add_argument(
        "--no_significant_filter", action="store_true",
        help="Jeśli podane, w KofamScan akceptuj też wyniki bez znacznika '*'."
    )
    args = parser.parse_args()

    # Walidacja argumentów
    if not args.genomes_dir and not (args.bakta_dir or args.kofam_dir):
        sys.exit(
            "[BŁĄD] Podaj --genomes_dir LUB przynajmniej jeden z --bakta_dir / --kofam_dir."
        )

    global KOFAM_SIGNIFICANT_ONLY
    if args.no_significant_filter:
        KOFAM_SIGNIFICANT_ONLY = False

    # 1. Wczytaj listę szukanych genów
    print(f"[1/4] Wczytywanie listy genów z: {args.genes}")
    gene_list = load_gene_list(args.genes)
    print(f"      Znaleziono {len(gene_list)} genów do wyszukania.")
    n_with_ko = gene_list["ko"].notna().sum()
    print(f"      Geny z K-numerem: {n_with_ko}, bez K-numeru: {len(gene_list) - n_with_ko}")

    # 2. Znajdź pliki genomów
    print(f"\n[2/4] Wykrywanie plików genomów...")
    genomes = find_genome_files(args.genomes_dir, args.bakta_dir, args.kofam_dir)
    if not genomes:
        sys.exit("[BŁĄD] Nie znaleziono żadnych genomów. Sprawdź ścieżki.")
    print(f"      Znaleziono {len(genomes)} genomów.")

    # 3. Przetwarzaj każdy genom
    print(f"\n[3/4] Wyszukiwanie genów w genomach...")
    matrix_rows = []

    for genome in genomes:
        name = genome["genome_name"]
        print(f"  → {name}")

        # Wczytaj adnotacje bakty
        if genome["bakta_path"]:
            bakta_df = load_bakta_tsv(genome["bakta_path"])
            print(f"      bakta: {len(bakta_df)} wpisów z {genome['bakta_path']}")
        else:
            print(f"      [OSTRZEŻENIE] Brak pliku bakty dla {name}")
            bakta_df = pd.DataFrame(columns=["gene_name", "ko_numbers", "raw_product"])

        # Wczytaj wyniki KofamScan
        if genome["kofam_path"]:
            kofam_ko_set = load_kofamscan_tsv(genome["kofam_path"])
            print(f"      kofam: {len(kofam_ko_set)} unikalnych K-numerów")
        else:
            print(f"      [OSTRZEŻENIE] Brak pliku KofamScan dla {name}")
            kofam_ko_set = set()

        # Wyszukaj geny
        results = search_genes_in_genome(gene_list, bakta_df, kofam_ko_set)
        results["genome"] = name
        matrix_rows.append(results)

    # 4. Zbuduj i zapisz macierz
    print(f"\n[4/4] Zapisywanie wyników do: {args.output}")
    matrix_df = pd.DataFrame(matrix_rows)

    # Przesuń kolumnę 'genome' na pierwszą pozycję
    cols = ["genome"] + [c for c in matrix_df.columns if c != "genome"]
    matrix_df = matrix_df[cols]

    matrix_df.to_csv(args.output, sep="\t", index=False)

    # Podsumowanie
    gene_cols = [c for c in matrix_df.columns if c != "genome"]
    total_cells = len(matrix_df) * len(gene_cols)
    present_cells = matrix_df[gene_cols].sum().sum()
    print(f"\n✓ Zapisano macierz: {len(matrix_df)} genomów × {len(gene_cols)} genów")
    print(f"  Obecność: {present_cells}/{total_cells} ({100*present_cells/total_cells:.1f}%)")

    # Geny nieznalezione w żadnym genomie
    absent_genes = [g for g in gene_cols if matrix_df[g].sum() == 0]
    if absent_genes:
        print(f"\n  [INFO] Geny nieznalezione w żadnym genomie ({len(absent_genes)}):")
        for g in absent_genes:
            print(f"    - {g}")

    # Opcjonalnie: szczegółowe statystyki per gen
    print(f"\n  Statystyki per gen (top 10 najczęstszych):")
    gene_counts = matrix_df[gene_cols].sum().sort_values(ascending=False)
    for gene, count in gene_counts.head(10).items():
        bar = "█" * int(count) + "░" * (len(matrix_df) - int(count))
        print(f"    {gene:<20} {bar} {int(count)}/{len(matrix_df)}")

    print(f"\nGotowe! Wynik zapisany w: {args.output}")


if __name__ == "__main__":
    main()