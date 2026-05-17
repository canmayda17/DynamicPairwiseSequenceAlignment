import argparse
from pathlib import Path

# The variable to make not reachable all nodes in the beginning.
NEG_INF = -10**18

def read_sequences(path):
    # Load sequences from one file.
    sequences = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            seq = "".join(line.split()).upper()
            if seq:
                sequences.append(seq)
    if len(sequences) != 2:
        raise ValueError(f"File must contain  2 sequences: {path}")
    return sequences[0], sequences[1]

def match_or_mismatch(a, b, match_score, mismatch_score):
    # Return match_score if bases match, else mismatch_score.
    if a == b:
        return match_score
    else:
        return mismatch_score
    
def find_best_candidate(candidates):
    # Pick the candidate with the highest score.
    best = candidates[0]

    for candidate in candidates:
        if candidate[0] > best[0]:
            best = candidate

    return best

def create_matrix(rows, cols, value):
    # Build a 2D matrix filled with the same initial value.
    matrix = []

    for i in range(rows):
        row = []

        for j in range(cols):
            row.append(value)

        matrix.append(row)

    return matrix

def global_align(seq1, seq2, match_score, mismatch_score, gap_open, gap_extend):

    len1 = len(seq1)
    len2 = len(seq2)

    match_matrix = create_matrix(len1 + 1, len2 + 1, NEG_INF)       # match_matrix[i][j] -> diagonal move: seq1[i-1] with seq2[j-1]
    gap_in_seq2_matrix = create_matrix(len1 + 1, len2 + 1, NEG_INF) # gap_in_seq2_matrix[i][j] -> down move: seq1[i-1] with '-'
    gap_in_seq1_matrix = create_matrix(len1 + 1, len2 + 1, NEG_INF) # gap_in_seq1_matrix[i][j] -> right move: '-' with seq2[j-1]

    # The backtrack matrices store which previous state chosen.
    backtrack_match = create_matrix(len1 + 1, len2 + 1, None)
    backtrack_gap_in_seq2 = create_matrix(len1 + 1, len2 + 1, None)
    backtrack_gap_in_seq1 = create_matrix(len1 + 1, len2 + 1, None)

    # Aligning two empty prefixes costs zero.
    match_matrix[0][0] = 0

    # First column: consume seq1 only; gaps appear in seq2.
    for i in range(1, len1 + 1):
        if i == 1:
            # First gap in a sequence has a gap_open penalty.
            gap_in_seq2_matrix[i][0] = gap_open
            backtrack_gap_in_seq2[i][0] = "match"
        else:
            # Extending an existing gap has a gap_extend penalty.
            gap_in_seq2_matrix[i][0] = gap_in_seq2_matrix[i - 1][0] + gap_extend
            backtrack_gap_in_seq2[i][0] = "gap_in_seq2"

    # First row: consume seq2 only; gaps appear in seq1.
    for j in range(1, len2 + 1):
        if j == 1:
            gap_in_seq1_matrix[0][j] = gap_open
            backtrack_gap_in_seq1[0][j] = "match"
        else:
            gap_in_seq1_matrix[0][j] = gap_in_seq1_matrix[0][j - 1] + gap_extend
            backtrack_gap_in_seq1[0][j] = "gap_in_seq1"

    # Fill the table cell by cell.
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            # Score for pairing seq1[i-1] with seq2[j-1].
            diagonal_score = match_or_mismatch(
                seq1[i - 1],
                seq2[j - 1],
                match_score,
                mismatch_score,
            )

            # Step for match/mismatch (diagonal).
            diagonal_candidates = [
                (match_matrix[i - 1][j - 1] + diagonal_score, "match"),
                (gap_in_seq2_matrix[i - 1][j - 1] + diagonal_score, "gap_in_seq2"),
                (gap_in_seq1_matrix[i - 1][j - 1] + diagonal_score, "gap_in_seq1"),
            ]
            best_diagonal = find_best_candidate(diagonal_candidates)
            match_matrix[i][j] = best_diagonal[0]
            backtrack_match[i][j] = best_diagonal[1]

            # Step for that seq1[i-1] aligned to a gap in seq2 (move down).
            down_candidates = [
                (match_matrix[i - 1][j] + gap_open, "match"),
                (gap_in_seq2_matrix[i - 1][j] + gap_extend, "gap_in_seq2"),
                (gap_in_seq1_matrix[i - 1][j] + gap_open, "gap_in_seq1"),
            ]
            best_down = find_best_candidate(down_candidates)
            gap_in_seq2_matrix[i][j] = best_down[0]
            backtrack_gap_in_seq2[i][j] = best_down[1]

            # Step for that seq2[j-1] aligned to a gap in seq1 (move right).
            right_candidates = [
                (match_matrix[i][j - 1] + gap_open, "match"),
                (gap_in_seq2_matrix[i][j - 1] + gap_open, "gap_in_seq2"),
                (gap_in_seq1_matrix[i][j - 1] + gap_extend, "gap_in_seq1"),
            ]
            best_right = find_best_candidate(right_candidates)
            gap_in_seq1_matrix[i][j] = best_right[0]
            backtrack_gap_in_seq1[i][j] = best_right[1]

    # Best score at the bottom-right corner is chosen.
    final_candidates = [
        (match_matrix[len1][len2], "match"),
        (gap_in_seq2_matrix[len1][len2], "gap_in_seq2"),
        (gap_in_seq1_matrix[len1][len2], "gap_in_seq1"),
    ]
    best_score, state = find_best_candidate(final_candidates)

    # Traceback from (len1, len2) to (0, 0) using stored backtrack states.
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
            print("There is an error in the backtracking")

    # Traceback built strings backward, with reverse it returned to normal order (AC-T -> T-CA(to show on output)).
    aligned1.reverse()
    aligned2.reverse()
    return best_score, "".join(aligned1), "".join(aligned2)

def alignment_middle_line(a1, a2):
    # "|" for matches, "x" for mismatches, " " for indels, also update count variables.
    chars = []
    match_count = 0
    mismatch_count = 0
    gap_open_count = 0
    gap_extend_count = 0
    prev_x = None
    prev_y = None
    for x, y in zip(a1, a2):
        if x == y and x != "-":
            chars.append("|")
            match_count += 1
        elif x != y and x != "-" and y != "-":
            chars.append("x")
            mismatch_count += 1
        else:
            if x == "-" and prev_x != "-":
                gap_open_count += 1
            elif y == "-" and prev_y != "-":
                gap_open_count += 1
            elif x == "-" and prev_x == "-":
                gap_extend_count += 1
            elif y == "-" and prev_y == "-":
                gap_extend_count += 1
            chars.append(" ")

        prev_x = x
        prev_y = y

    return "".join(chars), match_count, mismatch_count, gap_open_count, gap_extend_count

def format_alignment(a1, a2, width):
    # Split long alignments into readable blocks of fixed width(100).
    middle, match_count, mismatch_count, gap_open_count, gap_extend_count = alignment_middle_line(a1, a2)
    blocks = []
    for start in range(0, len(a1), width):
        end = start + width
        blocks.append(a1[start:end])
        blocks.append(middle[start:end])
        blocks.append(a2[start:end])
        blocks.append("")
    alignment_text = "\n".join(blocks).rstrip()
    return alignment_text, match_count, mismatch_count, gap_open_count, gap_extend_count

def build_report(name1, seq1, name2, seq2, score, a1, a2, width):
    alignment_text, match_count, mismatch_count, gap_open_count, gap_extend_count = format_alignment(a1, a2, width)
    return "\n".join([
        f"-----------------------------------------",
        f"Length for {name1}= {len(seq1)}",
        f"Length for {name2}= {len(seq2)}",
        f"-----------------------------------------",
        f"Match count= {match_count}",
        f"Mismatch count= {mismatch_count}",
        f"Gap opening count= {gap_open_count}",
        f"Gap extension count= {gap_extend_count}",
        f"-----------------------------------------",
        f"Alignment score= {score}",
        "",
        alignment_text,
    ])

def parse_args():
    parser = argparse.ArgumentParser()

    # To ensure take input as one .seq file with exactly two sequences (two lines).
    parser.add_argument("file")

    return parser.parse_args()

def main():
    """Read two sequences, align them, print results, and write a report file."""

    args = parse_args()

    # Assignment default scoring scheme.
    match_score = 7
    mismatch_score = -5
    gap_open = -3
    gap_extend = -1
    width = 100

    path = Path(args.file)
    if not path.exists():
        raise SystemExit(f"No such file: {args.file}")
    
    # 2 sequences are read from the input file, and their names are set.
    seq1, seq2 = read_sequences(path)
    name1 = "Sequence 1"
    name2 = "Sequence 2"

    blocks = []

    # Run global affine-gap alignment on the two sequences.
    score, a1, a2 = global_align(seq1, seq2, match_score, mismatch_score, gap_open, gap_extend)
    blocks.append(build_report(name1, seq1, name2, seq2, score, a1, a2, width))

    # Header shown before the alignment block which includes file name and scoring information.
    header = "\n".join([
        f"-----------------------------------------",
        f"Input file= {path.name}",
        f"-----------------------------------------",
        f"match score= {match_score}",
        f"mismatch score= {mismatch_score}",
        f"gap opening penalty= {gap_open}",
        f"gap extension penalty= {gap_extend}",
        "",
    ])

    # Adding lengths for both sequences and alignment score with showing alignment.
    full_output = header + "\n\n".join(blocks)

    print(full_output)

    output_path = f"{path.stem}_results.txt"

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(full_output)

    print(f"\nResults are written to: {output_path}") 
    
if __name__ == "__main__":
    main()
