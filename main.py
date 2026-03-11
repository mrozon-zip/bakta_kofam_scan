"""
gene_presence_absence.py
========================
Wyszukuje obecność określonych genów w adnotowanych genomach bakteryjnych.
Źródła danych: pliki .tsv z bakty + wyniki KofamScan (.tsv)

Wynik: dla każdego genomu osobny plik TSV z listą OBECNYCH genów,
       wzbogacony o kolumny 'function' i 'trait' z tabeli zapytań.

Wymagana struktura tabeli zapytań (genes.tsv / query_table.tsv):
    gene    ko       function                            Trait
    cbh     K01442   Choloylglycine hydrolase (EC...)    Bile tolerance
    luxS    K07173   S-ribosylhomocysteine lyase...      Anti-pathogenic effect

Struktura katalogów (tryb podkatalogów):
    genomes/
    ├── genome1/
    │   ├── genome1.tsv          # adnotacja bakty
    │   └── genome1_kofam.tsv    # wyniki KofamScan
    ├── genome2/
    │   └── ...

Użycie:
    # Tryb podkatalogów:
    python gene_presence_absence.py --genes query_table.tsv --genomes_dir genomes/

    # Tryb płaski:
    python gene_presence_absence.py --genes query_table.tsv \
        --bakta_dir bakta_results/ --kofam_dir kofam_results/ --output_dir results/
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Konfiguracja kolumn bakty — dostosuj jeśli Twoje pliki mają inne nagłówki
# ---------------------------------------------------------------------------
BAKTA_COL_GENE    = "gene"
BAKTA_COL_PRODUCT = "product"
BAKTA_COL_DBXREF  = "db_xrefs"

# KofamScan: indeksy kolumn w pliku wynikowym
KOFAM_COL_SIGNIFICANT = 0   # '*' = przeszło próg
KOFAM_COL_GENE_ID     = 1   # nazwa sekwencji
KOFAM_COL_KO          = 2   # K-numer

# Czy akceptować tylko wyniki KofamScan z '*' (przejście progu HMM)?
KOFAM_SIGNIFICANT_ONLY = True


# ---------------------------------------------------------------------------
# Parsowanie tabeli zapytań (genes.tsv / query_table.tsv)
# ---------------------------------------------------------------------------

def load_gene_list(path: str) -> pd.DataFrame:
    """
    Wczytuje tabelę szukanych genów.

    Wymagane kolumny: 'gene'
    Opcjonalne:       'ko' (K-numer KEGG), 'function', 'Trait' (lub 'trait')

    Wszystkie dodatkowe kolumny są zachowywane i trafią do pliku wynikowego.
    Wewnętrznie dodawane są kolumny pomocnicze _gene_lower i _ko (usuwane z outputu).
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    col_lower = {c.lower(): c for c in df.columns}

    if "gene" not in col_lower:
        sys.exit(f"[BŁĄD] Plik '{path}' musi zawierać kolumnę 'gene'.")

    gene_col = col_lower["gene"]
    df["_gene_lower"] = df[gene_col].str.strip().str.lower()

    # Normalizacja K-numerów
    if "ko" in col_lower:
        ko_col = col_lower["ko"]
        df["_ko"] = df[ko_col].str.strip().str.upper()
        df.loc[~df["_ko"].str.match(r"^K\d{5}$", na=False), "_ko"] = pd.NA
    else:
        print("[INFO] Brak kolumny 'ko' — wyszukiwanie tylko po nazwach genów.")
        df["_ko"] = pd.NA

    df = df.drop_duplicates(subset=["_gene_lower"]).reset_index(drop=True)

    n_with_ko = df["_ko"].notna().sum()
    print(f"      Wczytano {len(df)} genów ({n_with_ko} z K-numerem, {len(df)-n_with_ko} bez).")
    return df


# ---------------------------------------------------------------------------
# Parsowanie pliku TSV bakty
# ---------------------------------------------------------------------------

def load_bakta_tsv(path: str) -> pd.DataFrame:
    """
    Wczytuje plik TSV z bakty.
    Zwraca DataFrame z kolumnami: gene_name, ko_numbers (set), raw_product.
    """
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception as e:
        print(f"  [OSTRZEŻENIE] Nie można wczytać {path}: {e}")
        return pd.DataFrame(columns=["gene_name", "ko_numbers", "raw_product"])

    df.columns = [c.strip().lstrip("#").strip().lower().replace(" ", "_") for c in df.columns]

    actual_gene_col    = next((c for c in df.columns if c == "gene"), None) or \
                         next((c for c in df.columns if "gene" in c), None)
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

        ko_numbers = set()
        if actual_dbxref_col and pd.notna(row.get(actual_dbxref_col)):
            found = re.findall(r"K\d{5}", str(row[actual_dbxref_col]), re.IGNORECASE)
            ko_numbers = {k.upper() for k in found}

        if gene_name or ko_numbers:
            records.append({
                "gene_name":   gene_name,
                "ko_numbers":  ko_numbers,
                "raw_product": product,
            })

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Parsowanie wyników KofamScan
# ---------------------------------------------------------------------------

def load_kofamscan_tsv(path: str) -> set:
    """
    Wczytuje wyniki KofamScan, zwraca zbiór K-numerów obecnych w genomie.
    """
    ko_numbers = set()
    try:
        with open(path, "r") as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("##") or line.startswith("#") or not line.strip():
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
# Wyszukiwanie genów — zwraca wiersze z query_table dla obecnych genów
# ---------------------------------------------------------------------------

def search_genes_in_genome(
    gene_list: pd.DataFrame,
    bakta_df: pd.DataFrame,
    kofam_ko_set: set,
) -> pd.DataFrame:
    """
    Zwraca podzbiór wierszy gene_list odpowiadający genom OBECNYM w genomie.

    Dopasowanie (wystarczy jeden z kanałów):
      (A) nazwa genu == nazwa w adnotacji bakty (lub bakta-nazwa zaczyna się od query)
      (B) K-numer genu jest w db_xrefs bakty
      (C) K-numer genu jest w wynikach KofamScan

    Kolumny pomocnicze _gene_lower i _ko są usuwane z outputu.
    """
    bakta_gene_names = set(bakta_df["gene_name"].dropna())
    bakta_ko_all = set()
    for ko_set in bakta_df["ko_numbers"]:
        bakta_ko_all.update(ko_set)
    all_ko = bakta_ko_all | kofam_ko_set

    present_indices = []
    for idx, row in gene_list.iterrows():
        query_gene = row["_gene_lower"]
        query_ko   = row["_ko"]
        found = False

        # Kanał 1: nazwa genu
        if query_gene:
            for bname in bakta_gene_names:
                if bname == query_gene or bname.startswith(query_gene):
                    found = True
                    break

        # Kanał 2: K-numer (bakta db_xrefs + KofamScan)
        if not found and pd.notna(query_ko) and str(query_ko) in all_ko:
            found = True

        if found:
            present_indices.append(idx)

    result = gene_list.loc[present_indices].drop(columns=["_gene_lower", "_ko"])
    return result.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Wykrywanie plików genomów
# ---------------------------------------------------------------------------

def find_genome_files(genomes_dir, bakta_dir, kofam_dir) -> list[dict]:
    genomes = []

    if genomes_dir:
        base = Path(genomes_dir)
        for genome_dir in sorted(base.iterdir()):
            if not genome_dir.is_dir():
                continue
            name = genome_dir.name
            bakta_files = [f for f in genome_dir.glob("*.tsv") if "kofam" not in f.name.lower()]
            kofam_files = [f for f in genome_dir.glob("*.tsv") if "kofam" in f.name.lower()]
            if not kofam_files:
                kofam_files = list(genome_dir.glob("*kofam*.txt"))
            genomes.append({
                "genome_name": name,
                "genome_dir":  str(genome_dir),
                "bakta_path":  str(bakta_files[0]) if bakta_files else None,
                "kofam_path":  str(kofam_files[0]) if kofam_files else None,
            })

    elif bakta_dir or kofam_dir:
        bakta_files = {f.stem: f for f in Path(bakta_dir).glob("*.tsv")} if bakta_dir else {}
        kofam_files = {}
        if kofam_dir:
            kofam_files = {f.stem: f for f in Path(kofam_dir).glob("*.tsv")}
            kofam_files.update({f.stem: f for f in Path(kofam_dir).glob("*.txt")})
        for name in sorted(set(bakta_files) | set(kofam_files)):
            genomes.append({
                "genome_name": name,
                "genome_dir":  None,
                "bakta_path":  str(bakta_files[name]) if name in bakta_files else None,
                "kofam_path":  str(kofam_files[name]) if name in kofam_files else None,
            })

    return genomes


# ---------------------------------------------------------------------------
# Główna funkcja
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Wyszukuje geny w genomach bakteryjnych, zapisuje tabelę obecnych genów per genom.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--genes", required=True,
        help="Plik TSV z listą szukanych genów (wymagana kol. 'gene', "
             "opcjonalne: 'ko', 'function', 'Trait').")
    parser.add_argument("--genomes_dir", default=None,
        help="Katalog z podkatalogami dla każdego genomu.")
    parser.add_argument("--bakta_dir", default=None,
        help="Katalog z plikami TSV bakty (tryb płaski).")
    parser.add_argument("--kofam_dir", default=None,
        help="Katalog z plikami KofamScan (tryb płaski).")
    parser.add_argument("--output_dir", default=None,
        help="Katalog wyjściowy (tryb płaski). Domyślnie: bieżący katalog.")
    parser.add_argument("--output_filename", default="genes_present.tsv",
        help="Nazwa pliku wynikowego per genom. Domyślnie: genes_present.tsv")
    parser.add_argument("--no_significant_filter", action="store_true",
        help="Akceptuj też wyniki KofamScan bez znacznika '*'.")
    args = parser.parse_args()

    if not args.genomes_dir and not (args.bakta_dir or args.kofam_dir):
        sys.exit("[BŁĄD] Podaj --genomes_dir LUB przynajmniej --bakta_dir / --kofam_dir.")

    global KOFAM_SIGNIFICANT_ONLY
    if args.no_significant_filter:
        KOFAM_SIGNIFICANT_ONLY = False

    # 1. Wczytaj tabelę zapytań
    print(f"[1/3] Wczytywanie tabeli genów z: {args.genes}")
    gene_list = load_gene_list(args.genes)

    # 2. Znajdź pliki genomów
    print(f"\n[2/3] Wykrywanie plików genomów...")
    genomes = find_genome_files(args.genomes_dir, args.bakta_dir, args.kofam_dir)
    if not genomes:
        sys.exit("[BŁĄD] Nie znaleziono żadnych genomów. Sprawdź ścieżki.")
    print(f"      Znaleziono {len(genomes)} genomów.")

    # 3. Przetwarzaj każdy genom i zapisuj wynik
    print(f"\n[3/3] Wyszukiwanie genów i zapis wyników...")
    summary = []

    for genome in genomes:
        name = genome["genome_name"]
        print(f"\n  → {name}")

        if genome["bakta_path"]:
            bakta_df = load_bakta_tsv(genome["bakta_path"])
            print(f"      bakta:  {len(bakta_df)} wpisów")
        else:
            print(f"      [OSTRZEŻENIE] Brak pliku bakty")
            bakta_df = pd.DataFrame(columns=["gene_name", "ko_numbers", "raw_product"])

        if genome["kofam_path"]:
            kofam_ko_set = load_kofamscan_tsv(genome["kofam_path"])
            print(f"      kofam:  {len(kofam_ko_set)} unikalnych K-numerów")
        else:
            print(f"      [OSTRZEŻENIE] Brak pliku KofamScan")
            kofam_ko_set = set()

        result_df = search_genes_in_genome(gene_list, bakta_df, kofam_ko_set)
        n_found = len(result_df)
        print(f"      Znaleziono: {n_found}/{len(gene_list)} genów")

        # Ustal ścieżkę wyjściową
        if genome["genome_dir"]:
            # Tryb podkatalogów: zapisz do katalogu genomu (jak w Twoim skrypcie)
            out_path = Path(genome["genome_dir"]) / args.output_filename
        else:
            # Tryb płaski: zapisz do output_dir lub bieżącego katalogu
            out_dir = Path(args.output_dir) if args.output_dir else Path(".")
            out_dir.mkdir(parents=True, exist_ok=True)
            out_path = out_dir / f"{name}_{args.output_filename}"

        result_df.to_csv(out_path, sep="\t", index=False)
        print(f"      Zapisano:  {out_path}")

        summary.append({
            "genome":      name,
            "genes_found": n_found,
            "genes_total": len(gene_list),
        })

    # Podsumowanie końcowe
    print("\n" + "=" * 60)
    print(f"  {'Genom':<30} {'Wynik':>6}   Wizualizacja")
    print("-" * 60)
    for s in summary:
        n, total = s["genes_found"], s["genes_total"]
        bar = "█" * n + "░" * (total - n)
        print(f"  {s['genome']:<30} {n:>3}/{total:<3}  {bar}")
    print("=" * 60)
    print(f"\nGotowe! Pliki zapisane jako '{args.output_filename}' w katalogach genomów.")


if __name__ == "__main__":
    main()