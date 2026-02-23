#include <vector>
#include <cstddef>
#include <algorithm>
#include "rect.h"

static inline bool in_bounds(int r, int c, int R, int C) {
	return (r >= 0 && r < R && c >= 0 && c < C);
}

// Returns maximal (greedy) contiguous, non-overlapping all-true rectangles
// with width >= minW and height >= minH.
void greedy_rects(
	const std::vector<std::vector<bool>>& mat,
	int minW,
	int minH,
	std::vector<Rect> &out
) {
	out.clear();
	const int R = static_cast<int>(mat.size());
	if (R == 0) return;
	const int C = static_cast<int>(mat[0].size());
	if (C == 0) return;

	// Defensive: ensure rectangular input
	for (int r = 1; r < R; ++r) {
		if (static_cast<int>(mat[r].size()) != C) {
			// Non-rectangular; simplest policy: treat as no rectangles.
			return;
		}
	}

	// consumed[r][c] == true means we've already assigned this true cell to a rectangle
	std::vector<std::vector<unsigned char>> consumed(R, std::vector<unsigned char>(C, 0));

	auto is_available_true = [&](int r, int c) -> bool {
		return mat[r][c] && !consumed[r][c];
	};

	auto row_supports = [&](int r, int c0, int w) -> bool {
		// check cells (r, c0..c0+w-1) are true and not consumed
		for (int c = c0; c < c0 + w; ++c) {
			if (!is_available_true(r, c)) return false;
		}
		return true;
	};

	for (int r = 0; r < R; ++r) {
		for (int c = 0; c < C; ++c) {
			if (!is_available_true(r, c)) continue;

			// 1) Maximize width on the top row starting at (r,c)
			int w = 0;
			while (c + w < C && is_available_true(r, c + w)) {
				++w;
			}

			// 2) Maximize height with that fixed width
			int h = 0;
			while (r + h < R && row_supports(r + h, c, w)) {
				++h;
			}

			// Consume the rectangle cells regardless; otherwise you'd revisit them.
			for (int rr = r; rr < r + h; ++rr) {
				for (int cc = c; cc < c + w; ++cc) {
					consumed[rr][cc] = 1;
				}
			}

			if (w >= minW && h >= minH) {
				out.push_back(Rect{r, c, w, h});
			}
		}
	}
}