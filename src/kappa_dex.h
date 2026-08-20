#pragma once

class kappa_mermx;

// Minimal kappa_dex surface for Stage 1 (setup_kappa_qkmer_index).
class kappa_dex
	{
public:
	bool m_AddNeighborhood = false;
	const kappa_mermx *m_ptrScoreMx = 0;
	};
