# Author: Christopher Munro 2026.
# This code is released under the Creative Commons Attribution (CC BY) 4.0 license. 
# This lets others distribute, remix, adapt, and build upon this work, even commercially,
# as long as they credit the original creator.
# https://creativecommons.org/share-your-work/use-remix/cc-licenses/

import argparse
import numpy as np
from ortools.sat.python import cp_model


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def blocks(lambda_block):
    if lambda_block == 0:
        C = np.array(
            [
                [1, 0, 1, -1, -1, 1],
                [1, 1, 0, 1, -1, -1],
                [1, -1, 1, 0, 1, -1],
                [1, -1, -1, 1, 0, 1],
                [1, 1, -1, -1, 1, 0],
            ],
            dtype=np.int64,
        )
        Σ_block_sums = np.zeros((6, 6), dtype=np.int64)
        Σ_block_sums[1:, 1:] = 1
        Λ = np.zeros((5, 5), dtype=np.int64)
    else:
        C = np.array(
            [
                [0, 1, -1, -1, 1, 0],
                [1, 0, 1, -1, -1, 0],
                [-1, 1, 0, 1, -1, 0],
                [-1, -1, 1, 0, 1, 0],
                [1, -1, -1, 1, 0, 0],
            ],
            dtype=np.int64,
        )
        Σ_block_sums = np.zeros((6, 6), dtype=np.int64)
        Σ_block_sums[:5, :5] = 1
        Σ_block_sums[5, 5] = 5
        Λ = np.ones((5, 5), dtype=np.int64)

    # 5CC^T + ΛΛ^T = 25 I_5
    check(
        np.array_equal(
            5 * C @ C.T + Λ @ Λ.T,
            25 * np.eye(5, dtype=np.int64),
        ),
        "C and Λ do not satisfy the fixed-row equation",
    )
    return C, Σ_block_sums, Λ


def solve(lambda_block, workers, seed):
    C, Σ_block_sums, Λ = blocks(lambda_block)

    # Step 1: repeat each column of C five times to form Ψ.
    Ψ = np.repeat(C.T, 5, axis=0)
    C_column_products = C.T @ C
    model = cp_model.CpModel()

    # Step 2: s_ij is the block seed for Cyc(s_ij) in Σ.
    block_seeds = {}
    for block_row in range(6):
        for block_column in range(block_row, 6):
            block_seeds[block_row, block_column] = tuple(
                model.new_int_var(
                    -1,
                    1,
                    f"s_{block_row}_{block_column}_{seed_position}",
                )
                for seed_position in range(5)
            )

    def seed_entry(block_row, block_column, seed_position):
        seed_position %= 5
        if block_row <= block_column:
            return block_seeds[block_row, block_column][seed_position]
        return block_seeds[block_column, block_row][(-seed_position) % 5]

    # Symmetry gives s_ji by reversing s_ij 
    for block in range(6):
        model.add(seed_entry(block, block, 1) == seed_entry(block, block, 4))
        model.add(seed_entry(block, block, 2) == seed_entry(block, block, 3))

    # Every row of Cyc(s_ij) has sum sum(s_ij).
    for block_row in range(6):
        for block_column in range(block_row, 6):
            model.add(
                sum(
                    seed_entry(block_row, block_column, seed_position)
                    for seed_position in range(5)
                )
                == int(Σ_block_sums[block_row, block_column])
            )

    # This is the circulant form of ΣΣ^T + ΨΨ^T = 25 I_30.
    # row_displacement is the second row's position minus the first row's
    # position within their respective 5x5 block rows.
    for first_row_block in range(6):
        for second_row_block in range(first_row_block, 6):
            for row_displacement in range(5):
                products = []
                for column_block in range(6):
                    for seed_position in range(5):
                        product = model.new_int_var(
                            -1,
                            1,
                            "p_"
                            f"{first_row_block}_{second_row_block}_"
                            f"{row_displacement}_{column_block}_{seed_position}",
                        )
                        model.add_multiplication_equality(
                            product,
                            [
                                seed_entry(
                                    first_row_block,
                                    column_block,
                                    seed_position,
                                ),
                                seed_entry(
                                    second_row_block,
                                    column_block,
                                    seed_position - row_displacement,
                                ),
                            ],
                        )
                        products.append(product)

                required_inner_product = (
                    25
                    if first_row_block == second_row_block
                    and row_displacement == 0
                    else 0
                )
                model.add(
                    sum(products)
                    + int(C_column_products[first_row_block, second_row_block])
                    == required_inner_product
                )

    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = seed

    status = solver.solve(model)
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        raise SystemExit(
            f"no W(35,25) found: solver status {solver.status_name(status)}"
        )

    # Step 2, continued: build the circulant blocks from the solved seeds.
    Σ = np.zeros((30, 30), dtype=np.int64)
    for block_row in range(6):
        for block_column in range(6):
            for row_in_block in range(5):
                for column_in_block in range(5):
                    seed_position = (column_in_block - row_in_block) % 5
                    Σ[
                        5 * block_row + row_in_block,
                        5 * block_column + column_in_block,
                    ] = solver.value(
                        seed_entry(block_row, block_column, seed_position)
                    )

    # Step 3: verify the two block equations for W W^T = 25 I_35.
    target_30 = 25 * np.eye(30, dtype=np.int64)
    check(
        np.array_equal(Σ @ Σ.T + Ψ @ Ψ.T, target_30),
        "ΣΣ^T + ΨΨ^T is not 25 I_30",
    )
    check(
        np.array_equal(Σ @ Ψ + Ψ @ Λ, np.zeros((30, 5), dtype=np.int64)),
        "ΣΨ + ΨΛ is not zero",
    )

    # Assemble the block matrix: W = ((Σ, Ψ), (Ψ^T, Λ)).
    W = np.block([[Σ, Ψ], [Ψ.T, Λ]])
    target_35 = 25 * np.eye(35, dtype=np.int64)
    check(np.array_equal(W, W.T), "W is not symmetric")
    check(np.all(np.isin(W, (-1, 0, 1))), "W has a non-ternary entry")
    check(np.array_equal(W @ W.T, target_35), "W W^T is not 25 I_35")
    check(np.array_equal(W.T @ W, target_35), "W^T W is not 25 I_35")
    output_path = f"w35-{lambda_block}.txt"
    np.savetxt(output_path, W, fmt="%d")
    print(f"done: {output_path}", flush=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda-block", type=int, choices=(0, 1), required=True)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=1)
    arguments = parser.parse_args()

    solve(
        arguments.lambda_block,
        arguments.workers,
        arguments.seed,
    )
