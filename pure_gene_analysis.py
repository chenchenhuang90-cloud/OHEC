"""Pure gene-analysis calculations; data access and web handling live elsewhere."""

from __future__ import annotations

import math
import re
from itertools import combinations
from typing import Any, Mapping, Sequence


def _row_value(row: Mapping[str, Any], key: str) -> Any:
    """Read a row field without depending on database/model objects."""
    if key in row:
        return row[key]
    wanted = key.casefold()
    return next((value for name, value in row.items() if str(name).casefold() == wanted), None)


def _gene_key(row: Mapping[str, Any]) -> str:
    return str(_row_value(row, "gene") or "").strip().casefold()


def _finite_number(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def search_mutation_gene(
    query: str,
    mutation_data: Mapping[str, Sequence[Mapping[str, Any]]],
    selected_tables: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Find mutation-matching cell lines and shared genes from supplied DB rows.

    ``mutation_data`` is ``{cell_line: rows}``, where each row is a mapping with
    at least a ``gene`` field. If ``selected_tables`` is given, gene
    intersection is limited to those cell lines; otherwise it uses the tables
    whose rows contain the search query.
    """
    query = str(query or "").strip()
    if not query:
        raise ValueError("Gene query must not be empty.")

    query_key = query.casefold()
    found_tables = [
        table_name
        for table_name, rows in mutation_data.items()
        if any(query_key in str(_row_value(row, "gene") or "").casefold() for row in rows)
    ]

    if not selected_tables:
        display_tables = found_tables
    else:
        available = {str(name).casefold(): name for name in mutation_data}
        display_tables = []
        for requested in selected_tables:
            table_name = available.get(str(requested).strip().casefold())
            if table_name is not None and table_name not in display_tables:
                display_tables.append(table_name)
        found_tables = display_tables.copy()

    if not display_tables:
        return {
            "query": query,
            "found_tables": found_tables,
            "display_tables": [],
            "shared_genes": [],
        }

    genes_by_table: list[dict[str, Mapping[str, Any]]] = []
    for table_name in display_tables:
        gene_rows: dict[str, Mapping[str, Any]] = {}
        for row in mutation_data[table_name]:
            key = _gene_key(row)
            if key:
                gene_rows.setdefault(key, row)
        genes_by_table.append(gene_rows)

    shared_keys = set(genes_by_table[0])
    for gene_rows in genes_by_table[1:]:
        shared_keys.intersection_update(gene_rows)

    shared_genes = []
    for gene_key in shared_keys:
        row = dict(genes_by_table[0][gene_key])
        row["gene"] = _row_value(row, "gene")
        row["effects"] = [
            _row_value(gene_rows[gene_key], "gene_effect")
            for gene_rows in genes_by_table
        ]
        shared_genes.append(row)
    shared_genes.sort(
        key=lambda row: (
            min(
                (number for number in (_finite_number(value) for value in row["effects"])
                 if number is not None),
                default=float("inf"),
            ),
            str(row["gene"]).casefold(),
        )
    )

    return {
        "query": query,
        "found_tables": found_tables,
        "display_tables": display_tables,
        "shared_genes": shared_genes,
    }


def search_crispr_gene(
    gene_name: str,
    crispr_data: Mapping[str, Sequence[Mapping[str, Any]]],
    selected_tables: Sequence[str] | None = None,
) -> dict[str, list[Mapping[str, Any]]]:
    """Return exact gene rows from supplied CRISPR cell-line data."""
    gene_name = str(gene_name or "").strip()
    if not gene_name:
        raise ValueError("Gene name must not be empty.")

    if selected_tables is None:
        chosen_tables = list(crispr_data)
    else:
        available = {str(name).casefold(): name for name in crispr_data}
        chosen_tables = []
        for requested in selected_tables:
            table_name = available.get(str(requested).strip().casefold())
            if table_name is not None and table_name not in chosen_tables:
                chosen_tables.append(table_name)

    target = gene_name.casefold()
    return {
        table_name: [
            row for row in crispr_data[table_name]
            if _gene_key(row) == target
        ]
        for table_name in chosen_tables
        if any(_gene_key(row) == target for row in crispr_data[table_name])
    }


def _mann_whitney_u_pvalue(
    sample_a: Sequence[float],
    sample_b: Sequence[float],
) -> float | None:
    """Two-sided normal-approximation Mann–Whitney U with tie correction."""
    n1, n2 = len(sample_a), len(sample_b)
    if n1 < 2 or n2 < 2:
        return None

    tagged = sorted(
        [(float(value), 0) for value in sample_a]
        + [(float(value), 1) for value in sample_b],
        key=lambda item: item[0],
    )
    rank_sum_a = 0.0
    tie_sizes = []
    start = 0
    while start < len(tagged):
        end = start + 1
        while end < len(tagged) and tagged[end][0] == tagged[start][0]:
            end += 1
        average_rank = (start + 1 + end) / 2.0
        rank_sum_a += sum(
            average_rank for _, group in tagged[start:end] if group == 0
        )
        tie_sizes.append(end - start)
        start = end

    u1 = rank_sum_a - n1 * (n1 + 1) / 2.0
    total = n1 + n2
    mean_u = n1 * n2 / 2.0
    tie_term = sum(size**3 - size for size in tie_sizes if size > 1)
    variance = (n1 * n2 / 12.0) * (
        total + 1.0 - tie_term / (total * (total - 1.0))
    )
    if variance <= 0:
        return 1.0

    z_score = max(0.0, (abs(u1 - mean_u) - 0.5) / math.sqrt(variance))
    return max(0.0, min(1.0, math.erfc(z_score / math.sqrt(2.0))))


def compute_synthetic_lethality_pvalue(
    gene_name: str,
    selected_tables: Sequence[str],
    crispr_data: Mapping[str, Sequence[Mapping[str, Any]]],
    exact_combo_size: int | None = None,
) -> dict[str, Any]:
    """Calculate synthetic-lethality p-value using supplied CRISPR DB rows.

    For each candidate combination of selected cell lines, pooled gene_effect
    values are compared with pooled values from all unselected cell lines.
    The test and candidate preference rules match the web application's
    calculation: two-sided Mann–Whitney U normal approximation with continuity
    and tie corrections; prefer the largest significant combination, otherwise
    use the combination with the smallest p-value.
    """
    gene_name = str(gene_name or "").strip()
    if not gene_name:
        raise ValueError("Gene name must not be empty.")

    available = {str(name).casefold(): name for name in crispr_data}
    selected = []
    for requested in selected_tables:
        table_name = available.get(str(requested).strip().casefold())
        if table_name is None:
            raise ValueError(f"No CRISPR rows supplied for cell line {requested!r}.")
        if table_name not in selected:
            selected.append(table_name)
    if not selected:
        raise ValueError("At least one selected CRISPR cell line is required.")

    target = gene_name.casefold()
    values_by_table: dict[str, list[float]] = {}
    for table_name, rows in crispr_data.items():
        values = []
        for row in rows:
            if _gene_key(row) != target:
                continue
            effect = _finite_number(_row_value(row, "gene_effect"))
            if effect is not None:
                values.append(effect)
        if values:
            values_by_table[table_name] = values

    selected_used = [name for name in selected if name in values_by_table]
    other_tables = [
        name for name in crispr_data
        if name not in selected and name in values_by_table
    ]
    other_values = [
        value for name in other_tables for value in values_by_table[name]
    ]
    minimum_combo_size = max(2, math.ceil(len(selected) / 2.0 + 1.0))

    def failure(message: str, tested_count: int = 0) -> dict[str, Any]:
        return {
            "ok": False,
            "error": message,
            "p_value": None,
            "best_combo_tables": [],
            "best_combo_size": 0,
            "combo_mean": None,
            "combo_stddev": None,
            "total_mutant_table_count": len(selected),
            "min_combo_size": minimum_combo_size,
            "tested_combination_count": tested_count,
            "selected_sample_size": 0,
            "other_sample_size": len(other_values),
            "selected_tables_used": selected_used,
            "other_tables_used": other_tables,
        }

    if len(other_values) < 2:
        return failure("Comparison group needs at least 2 valid gene_effect values.")
    if len(selected_used) < minimum_combo_size:
        return failure("Not enough selected cell lines have valid gene_effect values.")

    requested_size = None
    if exact_combo_size is not None:
        try:
            requested_size = int(exact_combo_size)
        except (TypeError, ValueError) as error:
            raise ValueError("exact_combo_size must be an integer.") from error
        requested_size = max(2, min(requested_size, len(selected_used)))
        if requested_size < minimum_combo_size:
            return failure(
                f"Requested combo size {requested_size} is below minimum "
                f"allowed size {minimum_combo_size}."
            )

    combo_sizes = (
        [requested_size]
        if requested_size is not None
        else range(len(selected_used), minimum_combo_size - 1, -1)
    )
    candidates = []
    significant_candidates = []
    tested_count = 0
    for combo_size in combo_sizes:
        for combo_tables in combinations(selected_used, combo_size):
            tested_count += 1
            combo_values = [
                value for table_name in combo_tables
                for value in values_by_table.get(table_name, [])
            ]
            p_value = _mann_whitney_u_pvalue(combo_values, other_values)
            if p_value is None:
                continue
            candidate = (p_value, combo_size, combo_tables, combo_values)
            candidates.append(candidate)
            if p_value < 0.05:
                significant_candidates.append(candidate)

    if requested_size is not None:
        chosen = min(candidates, key=lambda item: item[0]) if candidates else None
    elif significant_candidates:
        chosen = min(significant_candidates, key=lambda item: (-item[1], item[0]))
    else:
        chosen = min(candidates, key=lambda item: item[0]) if candidates else None

    if chosen is None:
        return failure("Insufficient valid samples to calculate a combination p-value.", tested_count)

    p_value, _, best_tables, best_values = chosen
    combo_mean = sum(best_values) / len(best_values)
    variance = sum((value - combo_mean) ** 2 for value in best_values) / len(best_values)
    return {
        "ok": True,
        "error": "",
        "p_value": p_value,
        "best_combo_tables": list(best_tables),
        "best_combo_size": len(best_tables),
        "combo_mean": combo_mean,
        "combo_stddev": math.sqrt(variance),
        "total_mutant_table_count": len(selected),
        "min_combo_size": minimum_combo_size,
        "tested_combination_count": tested_count,
        "selected_sample_size": len(best_values),
        "other_sample_size": len(other_values),
        "selected_tables_used": selected_used,
        "other_tables_used": other_tables,
    }
























if __name__ == "__main__":
    # Minimal demo: example input structure
    demo_mutation_data = {
        "CellLineA": [{"gene": "TP53", "gene_effect": -0.6}, {"gene": "KRAS", "gene_effect": -0.3}],
        "CellLineB": [{"gene": "TP53", "gene_effect": -0.7}, {"gene": "PTEN", "gene_effect": -0.2}],
    }
    res = search_mutation_gene("TP53", demo_mutation_data)
    print("Demo mutation search result:")
    print(res)

    demo_crispr_data = {
        "CellLineA": [{"gene": "TP53", "gene_effect": -0.8}],
        "CellLineB": [{"gene": "TP53", "gene_effect": -0.9}],
        "CellLineC": [{"gene": "TP53", "gene_effect": -0.1}],
    }
    sl_res = compute_synthetic_lethality_pvalue(
        gene_name="TP53",
        selected_tables=["CellLineA", "CellLineB"],
        crispr_data=demo_crispr_data
    )
    print("\nDemo synthetic lethality result:")
    print(sl_res)
