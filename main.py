#!/usr/bin/env python3
import argparse
from pathlib import Path

NEG_INF = -10**18

def read_sequence(path):
    chars = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            chars.append("".join(line.split()))
    seq = "".join(chars).upper()
    if not seq:
        raise ValueError(f"File is empty: {path}")
    return seq

def match_or_mismatch(a, b, match_score, mismatch_score):
    if a == b:
        return match_score
    else:
        return mismatch_score
    
def find_best_candidate(candidates):
    best = candidates[0]

    for candidate in candidates:
        if candidate[0] > best[0]:
            best = candidate

    return best

def create_matrix(rows, cols, value):
    matrix = []

    for i in range(rows):
        row = []

        for j in range(cols):
            row.append(value)

        matrix.append(row)

    return matrix

def affine_global_align(seq1, seq2, match_score, mismatch_score, gap_open, gap_extend):
    """
    Global alignment with affine gap penalty.

    This still follows the slide's idea:
      score matrix + backtrack pointer

    But affine gap scoring needs three score matrices:
      match_matrix[i][j]       -> diagonal move: seq1[i-1] with seq2[j-1]
      gap_in_seq2_matrix[i][j] -> down move:     seq1[i-1] with '-'
      gap_in_seq1_matrix[i][j] -> right move:    '-' with seq2[j-1]

    The matching backtrack matrices store which previous state was chosen.
    """
    len1 = len(seq1)
    len2 = len(seq2)

    match_matrix = create_matrix(len1 + 1, len2 + 1, NEG_INF)
    gap_in_seq2_matrix = create_matrix(len1 + 1, len2 + 1, NEG_INF)
    gap_in_seq1_matrix = create_matrix(len1 + 1, len2 + 1, NEG_INF)

    backtrack_match = create_matrix(len1 + 1, len2 + 1, None)
    backtrack_gap_in_seq2 = create_matrix(len1 + 1, len2 + 1, None)
    backtrack_gap_in_seq1 = create_matrix(len1 + 1, len2 + 1, None)

    match_matrix[0][0] = 0

    for i in range(1, len1 + 1):
        if i == 1:
            gap_in_seq2_matrix[i][0] = gap_open
            backtrack_gap_in_seq2[i][0] = "match"
        else:
            gap_in_seq2_matrix[i][0] = gap_in_seq2_matrix[i - 1][0] + gap_extend
            backtrack_gap_in_seq2[i][0] = "gap_in_seq2"

    for j in range(1, len2 + 1):
        if j == 1:
            gap_in_seq1_matrix[0][j] = gap_open
            backtrack_gap_in_seq1[0][j] = "match"
        else:
            gap_in_seq1_matrix[0][j] = gap_in_seq1_matrix[0][j - 1] + gap_extend
            backtrack_gap_in_seq1[0][j] = "gap_in_seq1"

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            diagonal_score = match_or_mismatch(
                seq1[i - 1],
                seq2[j - 1],
                match_score,
                mismatch_score,
            )

            diagonal_candidates = [
                (match_matrix[i - 1][j - 1] + diagonal_score, "match"),
                (gap_in_seq2_matrix[i - 1][j - 1] + diagonal_score, "gap_in_seq2"),
                (gap_in_seq1_matrix[i - 1][j - 1] + diagonal_score, "gap_in_seq1"),
            ]
            best_diagonal = find_best_candidate(diagonal_candidates)
            match_matrix[i][j] = best_diagonal[0]
            backtrack_match[i][j] = best_diagonal[1]

            down_candidates = [
                (match_matrix[i - 1][j] + gap_open, "match"),
                (gap_in_seq2_matrix[i - 1][j] + gap_extend, "gap_in_seq2"),
                (gap_in_seq1_matrix[i - 1][j] + gap_open, "gap_in_seq1"),
            ]
            best_down = find_best_candidate(down_candidates)
            gap_in_seq2_matrix[i][j] = best_down[0]
            backtrack_gap_in_seq2[i][j] = best_down[1]

            right_candidates = [
                (match_matrix[i][j - 1] + gap_open, "match"),
                (gap_in_seq2_matrix[i][j - 1] + gap_open, "gap_in_seq2"),
                (gap_in_seq1_matrix[i][j - 1] + gap_extend, "gap_in_seq1"),
            ]
            best_right = find_best_candidate(right_candidates)
            gap_in_seq1_matrix[i][j] = best_right[0]
            backtrack_gap_in_seq1[i][j] = best_right[1]

    final_candidates = [
        (match_matrix[len1][len2], "match"),
        (gap_in_seq2_matrix[len1][len2], "gap_in_seq2"),
        (gap_in_seq1_matrix[len1][len2], "gap_in_seq1"),
    ]
    best_score, state = find_best_candidate(final_candidates)

    aligned1 = []
    aligned2 = []
    i = len1
    j = len2

    while i > 0 or j > 0:
        if state == "match":
            previous_state = backtrack_match[i][j]
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
            state = previous_state
        elif state == "gap_in_seq2":
            previous_state = backtrack_gap_in_seq2[i][j]
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
            state = previous_state
        elif state == "gap_in_seq1":
            previous_state = backtrack_gap_in_seq1[i][j]
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1
            state = previous_state
        else:
            raise RuntimeError("Backtrack state is missing.")

    aligned1.reverse()
    aligned2.reverse()
    return best_score, "".join(aligned1), "".join(aligned2)

def alignment_middle_line(a1, a2):
    chars = []
    for x, y in zip(a1, a2):
        if x == y and x != "-":
            chars.append("|")
        elif x != y and x != "-" and y != "-":
            chars.append("x")
        else:
            chars.append(" ")
    return "".join(chars)

def format_alignment(a1, a2, width):
    middle = alignment_middle_line(a1, a2)
    blocks = []
    for start in range(0, len(a1), width):
        end = start + width
        blocks.append(a1[start:end])
        blocks.append(middle[start:end])
        blocks.append(a2[start:end])
        blocks.append("")
    return "\n".join(blocks).rstrip()

def pairwise_jobs(items):
    yield items[0], items[1]

def build_report(name1, seq1, name2, seq2, score, a1, a2, width):
    return "\n".join([
        f"-----------------------------------------",
        f"Length for {name1}= {len(seq1)}",
        f"Length for {name2}= {len(seq2)}",
        f"-----------------------------------------",
        f"Alignment score= {score}",
        "",
        format_alignment(a1, a2, width),
    ])

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "files",
        nargs="*",
    )

    return parser.parse_args()

def main():

    args = parse_args()

    match_score = 7
    mismatch_score = -5
    gap_open = -3
    gap_extend = -1
    width = 100

    file_list = args.files
    if len(file_list) != 2:
        raise SystemExit("2 files are needed.")

    items = []
    for filename in file_list:
        path = Path(filename)
        if not path.exists():
            raise SystemExit(f"No such file: {filename}")
        seq = read_sequence(path)
        items.append((path.name, seq))

    blocks = []

    name1 = items[0][0]
    seq1 = items[0][1]
    name2 = items[1][0]
    seq2 = items[1][1]

    score, a1, a2 = affine_global_align(seq1, seq2, match_score, mismatch_score, gap_open, gap_extend)
    blocks.append(build_report(name1, seq1, name2, seq2, score, a1, a2, width))

    header = "\n".join([
        f"-----------------------------------------",
        f"Compared files= {Path(file_list[0]).stem} vs {Path(file_list[1]).stem}",
        f"-----------------------------------------",
        f"match score= {match_score}",
        f"mismatch score= {mismatch_score}",
        f"gap opening penalty= {gap_open}",
        f"gap extension penalty= {gap_extend}",
        "",
    ])

    full_output = header + "\n\n".join(blocks)

    print(full_output)

    names = []
    for filename in file_list:
        names.append(Path(filename).stem)

    output_path = f"alignment_scores_{names[0]}_{names[1]}.txt"
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(full_output)

    print(f"\nResults are written to: {output_path}")

if __name__ == "__main__":
    main()
