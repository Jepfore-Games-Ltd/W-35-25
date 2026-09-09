#!/usr/bin/env python3
"""Independently verify every representative in the W(35,25) catalogue."""

from __future__ import annotations

import argparse
import csv
import itertools
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from math import comb
from pathlib import Path

import igraph as ig
import numpy as np


ORDER = 35
WEIGHT = 25
EXPECTED_REPRESENTATIVES = 566
EXPECTED_FAMILIES = {
    "published": 41,
    "completion": 431,
    "rebanded": 94,
}
EXPECTED_RANKS = {
    8: 1,
    9: 7,
    10: 28,
    11: 45,
    12: 38,
    13: 18,
    14: 36,
    15: 138,
    16: 180,
    17: 75,
}
REQUIRED_METADATA = (
    "Lambda",
    "Method",
    "Rank_5",
    "pi(0)",
    "pi(1)",
    "Symmetric",
)


class VerificationError(RuntimeError):
    """A catalogue assertion failed."""


@dataclass(frozen=True)
class Representative:
    identifier: str
    path: Path
    matrix: np.ndarray
    rank_5: int
    row_profile: tuple[int, ...]
    column_profile: tuple[int, ...]

    @property
    def invariant_key(self) -> tuple[object, ...]:
        return (self.rank_5, self.row_profile, self.column_profile)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def parse_bool(value: str, context: str) -> bool:
    require(value in {"True", "False"}, f"{context}: invalid Boolean {value!r}")
    return value == "True"


def parse_matrix_file(path: Path) -> tuple[dict[str, str], np.ndarray]:
    metadata: dict[str, str] = {}
    rows: list[list[int]] = []
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(),
        1,
    ):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("#"):
            content = line[1:].strip()
            if " = " in content:
                key, value = content.split(" = ", 1)
                metadata[key] = value
            continue
        try:
            rows.append([int(value) for value in line.split()])
        except ValueError as error:
            raise VerificationError(
                f"{path}:{line_number}: matrix row is not integral"
            ) from error

    for key in REQUIRED_METADATA:
        require(key in metadata, f"{path}: missing '# {key} = ...' metadata")
    require(len(rows) == ORDER, f"{path}: expected {ORDER} rows, found {len(rows)}")
    require(
        all(len(row) == ORDER for row in rows),
        f"{path}: every matrix row must have {ORDER} entries",
    )
    matrix = np.asarray(rows, dtype=np.int64)
    require(
        bool(np.all(np.isin(matrix, (-1, 0, 1)))),
        f"{path}: matrix contains an entry outside -1, 0, 1",
    )
    return metadata, matrix


def rank_mod_5(matrix: np.ndarray) -> int:
    """Exact Gaussian elimination over F_5."""
    reduced = np.asarray(matrix, dtype=np.int64).copy() % 5
    row = 0
    for column in range(reduced.shape[1]):
        pivot = next(
            (
                candidate
                for candidate in range(row, reduced.shape[0])
                if reduced[candidate, column]
            ),
            None,
        )
        if pivot is None:
            continue
        reduced[[row, pivot]] = reduced[[pivot, row]]
        reduced[row] *= pow(int(reduced[row, column]), -1, 5)
        reduced[row] %= 5
        for other in range(reduced.shape[0]):
            if other != row and reduced[other, column]:
                reduced[other] -= reduced[other, column] * reduced[row]
                reduced[other] %= 5
        row += 1
        if row == reduced.shape[0]:
            break
    return row


def four_subsets(order: int) -> np.ndarray:
    return np.asarray(list(itertools.combinations(range(order), 4)), dtype=np.int16)


def four_profile(
    matrix: np.ndarray,
    subsets: np.ndarray,
    *,
    chunk_size: int = 16384,
) -> tuple[int, ...]:
    """Complete absolute four-row product-sum distribution."""
    order = matrix.shape[0]
    histogram = np.zeros(order + 1, dtype=np.int64)
    compact = np.asarray(matrix, dtype=np.int16)
    for start in range(0, len(subsets), chunk_size):
        block = subsets[start : start + chunk_size]
        products = compact[block[:, 0]] * compact[block[:, 1]]
        products *= compact[block[:, 2]]
        products *= compact[block[:, 3]]
        values = np.abs(products.sum(axis=1, dtype=np.int64))
        require(
            bool(np.all(values <= order)),
            "four-profile value exceeded the matrix order",
        )
        histogram += np.bincount(values, minlength=order + 1)
    require(
        int(histogram.sum()) == comb(order, 4),
        "four-profile histogram has the wrong total",
    )
    return tuple(int(value) for value in histogram)


def slow_four_profile(matrix: np.ndarray) -> tuple[int, ...]:
    order = matrix.shape[0]
    histogram = [0] * (order + 1)
    for indices in itertools.combinations(range(order), 4):
        value = abs(
            sum(
                int(matrix[indices[0], column])
                * int(matrix[indices[1], column])
                * int(matrix[indices[2], column])
                * int(matrix[indices[3], column])
                for column in range(order)
            )
        )
        histogram[value] += 1
    return tuple(histogram)


def lambda_description(matrix: np.ndarray) -> str:
    corner = matrix[-5:, -5:]
    if not np.any(corner):
        return "0"
    if np.all(corner == 1):
        return "J"
    if np.all(corner == -1):
        return "-J"
    if np.all(np.abs(corner) == 1) and rank_mod_5(corner) == 1:
        return "signed rank-one (J-equivalent)"
    return "other"


def signed_incidence_graph(matrix: np.ndarray) -> ig.Graph:
    """Graph whose coloured isomorphisms are signed H-equivalences."""
    require(
        matrix.ndim == 2 and matrix.shape[0] == matrix.shape[1],
        "signed-incidence graph requires a square matrix",
    )
    order = matrix.shape[0]
    row_parent = 0
    row_minus = order
    row_plus = 2 * order
    column_minus = 3 * order
    column_plus = 4 * order
    column_parent = 5 * order

    colors = (
        [0] * order
        + [1] * order
        + [1] * order
        + [2] * order
        + [2] * order
        + [3] * order
    )
    edges: list[tuple[int, int]] = []
    for index in range(order):
        edges.append((row_parent + index, row_minus + index))
        edges.append((row_parent + index, row_plus + index))
        edges.append((column_minus + index, column_parent + index))
        edges.append((column_plus + index, column_parent + index))

    for row in range(order):
        for column in range(order):
            entry = int(matrix[row, column])
            if entry == 0:
                continue
            edges.append(
                (
                    row_plus + row,
                    (column_plus if entry == 1 else column_minus) + column,
                )
            )
            edges.append(
                (
                    row_minus + row,
                    (column_minus if entry == 1 else column_plus) + column,
                )
            )

    graph = ig.Graph(n=6 * order, edges=edges, directed=False)
    graph.vs["kind"] = colors
    require(graph.vcount() == 6 * order, "wrong graph vertex count")
    expected_edges = 4 * order + 2 * int(np.count_nonzero(matrix))
    require(graph.ecount() == expected_edges, "wrong graph edge count")
    return graph


def canonical_graph_form(matrix: np.ndarray) -> tuple[object, ...]:
    """Exact coloured Bliss canonical form of the signed-incidence graph."""
    graph = signed_incidence_graph(matrix)
    permutation = graph.canonical_permutation(color=graph.vs["kind"])
    canonical = graph.permute_vertices(permutation)
    colors = tuple(int(value) for value in canonical.vs["kind"])
    edges = tuple(
        sorted(
            tuple(sorted((int(edge.source), int(edge.target))))
            for edge in canonical.es
        )
    )
    return colors, edges


def graph_h_equivalent(first: np.ndarray, second: np.ndarray) -> bool:
    return canonical_graph_form(first) == canonical_graph_form(second)


def brute_signed_canonical(matrix: np.ndarray) -> tuple[int, ...]:
    """Canonical form by exhaustive signed monomial action for tiny matrices."""
    order = matrix.shape[0]
    best: tuple[int, ...] | None = None
    for row_order in itertools.permutations(range(order)):
        for column_order in itertools.permutations(range(order)):
            permuted = matrix[np.ix_(row_order, column_order)]
            for row_signs in itertools.product((-1, 1), repeat=order):
                signed_rows = permuted * np.asarray(row_signs)[:, None]
                for column_signs in itertools.product((-1, 1), repeat=order):
                    candidate = signed_rows * np.asarray(column_signs)[None, :]
                    key = tuple(int(value) for value in candidate.ravel())
                    if best is None or key < best:
                        best = key
    require(best is not None, "failed to canonicalize a tiny matrix")
    return best


def run_small_self_tests() -> None:
    matrices = [
        np.asarray(values, dtype=np.int64).reshape(2, 2)
        for values in itertools.product((-1, 0, 1), repeat=4)
    ]
    brute_keys = [brute_signed_canonical(matrix) for matrix in matrices]
    graph_keys = [canonical_graph_form(matrix) for matrix in matrices]
    for first in range(len(matrices)):
        for second in range(first, len(matrices)):
            expected = brute_keys[first] == brute_keys[second]
            observed = graph_keys[first] == graph_keys[second]
            require(
                observed == expected,
                "signed-incidence graph failed exhaustive 2x2 self-test",
            )

    sample = np.asarray(
        [
            [1, 0, -1, 1, 0, -1],
            [0, 1, 1, -1, 0, 1],
            [-1, 1, 0, 0, 1, -1],
            [1, -1, 0, 1, -1, 0],
            [0, 0, 1, -1, 1, 1],
            [-1, 1, -1, 0, 1, 0],
        ],
        dtype=np.int64,
    )
    subsets = four_subsets(len(sample))
    require(
        four_profile(sample, subsets) == slow_four_profile(sample),
        "vectorized four-profile failed its slow-reference self-test",
    )


def deterministic_scramble(matrix: np.ndarray) -> np.ndarray:
    rng = np.random.default_rng(3525)
    row_order = rng.permutation(len(matrix))
    column_order = rng.permutation(len(matrix))
    row_signs = rng.choice((-1, 1), size=len(matrix))
    column_signs = rng.choice((-1, 1), size=len(matrix))
    return (
        matrix[np.ix_(row_order, column_order)]
        * row_signs[:, None]
        * column_signs[None, :]
    )


def load_manifest(repository: Path) -> list[dict[str, str]]:
    manifest_path = repository / "representatives" / "manifest.csv"
    require(manifest_path.is_file(), f"missing manifest: {manifest_path}")
    with manifest_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    require(
        len(rows) == EXPECTED_REPRESENTATIVES,
        f"manifest has {len(rows)} rows, expected {EXPECTED_REPRESENTATIVES}",
    )
    expected_ids = [f"H{number:03d}" for number in range(1, 567)]
    require(
        [row["id"] for row in rows] == expected_ids,
        "manifest IDs are not the complete ordered range H001-H566",
    )
    require(
        Counter(row["family"] for row in rows) == EXPECTED_FAMILIES,
        "manifest family counts do not match 41 + 431 + 94",
    )

    listed_paths = {Path(row["path"]).as_posix() for row in rows}
    require(len(listed_paths) == len(rows), "manifest contains duplicate paths")
    actual_paths = {
        path.relative_to(repository / "representatives").as_posix()
        for path in (repository / "representatives").rglob("H*.txt")
    }
    require(
        listed_paths == actual_paths,
        "manifest paths and representative files do not agree",
    )
    return rows


def verify_representative(
    repository: Path,
    manifest_row: dict[str, str],
    subsets: np.ndarray,
) -> Representative:
    relative_path = Path(manifest_row["path"])
    require(
        not relative_path.is_absolute() and ".." not in relative_path.parts,
        f"{manifest_row['id']}: unsafe manifest path",
    )
    path = repository / "representatives" / relative_path
    metadata, matrix = parse_matrix_file(path)
    identifier = manifest_row["id"]

    for header, manifest_field in (
        ("Lambda", "lambda"),
        ("Method", "method"),
        ("Rank_5", "rank5"),
        ("pi(0)", "pi0"),
        ("pi(1)", "pi1"),
        ("Symmetric", "symmetric"),
    ):
        require(
            metadata[header] == manifest_row[manifest_field],
            f"{identifier}: header {header} disagrees with manifest",
        )

    target = WEIGHT * np.eye(ORDER, dtype=np.int64)
    require(
        np.array_equal(matrix @ matrix.T, target),
        f"{identifier}: W W^T != 25 I",
    )
    require(
        np.array_equal(matrix.T @ matrix, target),
        f"{identifier}: W^T W != 25 I",
    )

    computed_rank = rank_mod_5(matrix)
    row_profile = four_profile(matrix, subsets)
    column_profile = four_profile(matrix.T, subsets)
    computed_symmetric = bool(np.array_equal(matrix, matrix.T))

    require(
        computed_rank == int(metadata["Rank_5"]),
        f"{identifier}: rank_5 metadata is wrong",
    )
    require(
        row_profile[0] == int(metadata["pi(0)"]),
        f"{identifier}: pi(0) metadata is wrong",
    )
    require(
        row_profile[1] == int(metadata["pi(1)"]),
        f"{identifier}: pi(1) metadata is wrong",
    )
    require(
        computed_symmetric == parse_bool(metadata["Symmetric"], identifier),
        f"{identifier}: symmetry metadata is wrong",
    )
    require(
        metadata["Lambda"] == lambda_description(matrix),
        f"{identifier}: Lambda metadata is wrong",
    )

    return Representative(
        identifier=identifier,
        path=path,
        matrix=matrix,
        rank_5=computed_rank,
        row_profile=row_profile,
        column_profile=column_profile,
    )


def verify_scramble_self_test(
    representative: Representative,
    subsets: np.ndarray,
) -> None:
    scrambled = deterministic_scramble(representative.matrix)
    require(
        graph_h_equivalent(representative.matrix, scrambled),
        "H001 signed-scramble graph self-test failed",
    )
    require(
        rank_mod_5(scrambled) == representative.rank_5,
        "rank changed under signed scramble",
    )
    require(
        four_profile(scrambled, subsets) == representative.row_profile,
        "row four-profile changed under signed scramble",
    )
    require(
        four_profile(scrambled.T, subsets) == representative.column_profile,
        "column four-profile changed under signed scramble",
    )


def verify_pairwise_inequivalence(
    representatives: list[Representative],
) -> tuple[int, int]:
    buckets: dict[tuple[object, ...], list[Representative]] = defaultdict(list)
    for representative in representatives:
        buckets[representative.invariant_key].append(representative)

    total_pairs = len(representatives) * (len(representatives) - 1) // 2
    graph_pairs = sum(
        len(bucket) * (len(bucket) - 1) // 2
        for bucket in buckets.values()
    )
    invariant_pairs = total_pairs - graph_pairs

    for bucket in buckets.values():
        if len(bucket) < 2:
            continue
        names = ", ".join(item.identifier for item in bucket)
        print(f"  exact signed-graph bucket: {names}", flush=True)
        graph_keys = {
            item.identifier: canonical_graph_form(item.matrix)
            for item in bucket
        }
        for first, second in itertools.combinations(bucket, 2):
            equivalent = (
                graph_keys[first.identifier] == graph_keys[second.identifier]
            )
            print(
                f"    {first.identifier} vs {second.identifier}: "
                f"{'EQUIVALENT' if equivalent else 'distinct'}",
                flush=True,
            )
            require(
                not equivalent,
                f"{first.identifier} and {second.identifier} are H-equivalent",
            )
    return invariant_pairs, graph_pairs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Verify all 566 W(35,25) H-class representatives.",
    )
    parser.add_argument(
        "--repository",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="repository root (default: parent of scripts/)",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=25,
        help="print progress after this many matrices (0 disables)",
    )
    arguments = parser.parse_args()
    repository = arguments.repository.resolve()

    started = time.perf_counter()
    print("Running verifier self-tests...", flush=True)
    run_small_self_tests()
    print("  exhaustive 2x2 graph test: passed", flush=True)
    print("  vectorized four-profile test: passed", flush=True)

    rows = load_manifest(repository)
    subsets = four_subsets(ORDER)
    representatives: list[Representative] = []
    print(f"Checking {len(rows)} representatives...", flush=True)
    for number, row in enumerate(rows, 1):
        representatives.append(verify_representative(repository, row, subsets))
        if arguments.progress_every and (
            number % arguments.progress_every == 0 or number == len(rows)
        ):
            print(f"  checked {number}/{len(rows)}", flush=True)

    verify_scramble_self_test(representatives[0], subsets)
    print("  H001 signed-scramble self-test: passed", flush=True)

    ranks = Counter(item.rank_5 for item in representatives)
    require(dict(sorted(ranks.items())) == EXPECTED_RANKS, "wrong rank histogram")
    print("Checking pairwise H-inequivalence...", flush=True)
    invariant_pairs, graph_pairs = verify_pairwise_inequivalence(representatives)
    total_pairs = len(representatives) * (len(representatives) - 1) // 2
    require(
        invariant_pairs + graph_pairs == total_pairs,
        "pair accounting failed",
    )

    elapsed = time.perf_counter() - started
    print()
    print("VERIFIED")
    print(f"  representatives: {len(representatives)}")
    print(f"  weighing identities and metadata: all passed")
    print(f"  pairs separated by exact invariants: {invariant_pairs}")
    print(f"  pairs checked by exact signed graph isomorphism: {graph_pairs}")
    print(f"  total pairwise H-inequivalence certificates: {total_pairs}")
    print(f"  rank histogram: {dict(sorted(ranks.items()))}")
    print(f"  elapsed seconds: {elapsed:.1f}")


if __name__ == "__main__":
    try:
        main()
    except VerificationError as error:
        raise SystemExit(f"VERIFICATION FAILED: {error}") from error
