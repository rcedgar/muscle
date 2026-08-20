#include "myutils.h"
#include "flat_chain_reader.h"
#include "cif.h"
#include "struct_desc.h"
#include <map>

void ChainizeLabel(string &Label, const string &ChainStr);

/***
https://mmcif.wwpdb.org/docs/user-guide/guide.html

PDBx/mmCIF uses data blocks to organize related information and data.
A data block is a logical partition of a data file designated by a data_ record.
A data block may be named by appending a text string after the data_ record and a
data block is terminated by either another data_ record or by the end of the file.

Data items are presented in two styles: key-value and tabular. 

An example of a key-value style where the PDBx/mmCIF item is followed directly by
a corresponding value:

_cell.entry_id	4HHB
_cell.length_a	63.150
_cell.length_b	83.590

In tabular style, a loop_ record is followed by rows of data item names and then
white-space delimited data values (e.g. ATOM).

The hash symbol (#) is used to separate categories to improve readability, 
but is not strictly necessary. It is also used to indicate comments.

https://mmcif.wwpdb.org/docs/tutorials/content/atomic-description.html
_atom_site.group_PDB
    The group of atoms to which the atom site belongs. This data
    item is provided for compatibility with the original Protein
    Data Bank format, and only for that purpose.
_atom_site.id
    A unique identifier for each atom position.
_atom_site.type_symbol
    The atom type.
_atom_site.auth_asym_id
    alternative to _atom_site.label_asym_id. This item holds the PDB chain identifier.
_atom_site.auth_atom_id
    alternative to _atom_site.label_atom_id. This item holds the PDB atom name.
_atom_site.auth_comp_id
    alternative to _atom_site.label_comp_id. This item holds the PDB 3-letter-code residue names
_atom_site.auth_seq_id
    alternative to _atom_site.label_seq_id. This item holds the 1-based PDB residue number.
_atom_site.pdbx_PDB_ins_code
    This data item corresponds to the PDB insertion code. 
_atom_site.pdbx_PDB_model_num
    This data item identifies the model number in an ensemble of coordinate data. 
_atom_site.group_PDB
    This data item is a place holder for the tags used by the PDB to identify coordinate records (e.g. ATOM or HETATM).
_atom_site.label_alt_id
    This item is a uniquely identifies for each alternative site for this atom position.
_atom_site.label_asym_id
    This data item is reference to item _struct_asym.id defined in category STRUCT_ASYM. This item identifies an instance of particular entity in the deposited coordinate set. For a structure determined by crystallographic method this corresponds to a unique identifier within the cyrstallographic asymmetric unit.
_atom_site.label_atom_id
    This data item is a reference to item _chem_comp_atom.atom_id defined in category CHEM_COMP_ATOM which is stored in the Chemical Component Dictionary. This atom identifier uniquely identifies each atom within each chemical component. 
_atom_site.label_comp_id
    This data item is a reference to item _chem_comp.id defined in category CHEM_COMP. This item is the primary identifier for chemical components which may either be mononers in a polymeric entity or complete non-polymer entities.
_atom_site.label_entity_id
    This data item is a reference to _entity.id defined in the ENTITY category. This item is used to identify chemically distinct portions of the molecular structure (e.g. polymer chains, ligands, solvent).
_atom_site.label_seq_id
    This data item is a reference to _entity_poly_seq.num defined in the ENTITY_POLY_SEQ category. This item is used to maintain the correspondence between the chemical sequence of a polymeric entity and the sequence information in the coordinate list and in may other structural categories. This identifier has no meaning for non-polymer entities.
_atom_site.Cartn_x _atom_site.Cartn_y _atom_site.Cartn_z
    Cartesian coordinate components describing the position of this atom site.

#
loop_
_atom_site.group_PDB # 0
_atom_site.id # 1
_atom_site.type_symbol # 2
_atom_site.label_atom_id # 3
_atom_site.label_alt_id # 4
_atom_site.label_comp_id # 5
_atom_site.label_asym_id # 6
_atom_site.label_entity_id # 7
_atom_site.label_seq_id # 8
_atom_site.pdbx_PDB_ins_code # 9
_atom_site.Cartn_x # 10
_atom_site.Cartn_y # 11
_atom_site.Cartn_z # 12
_atom_site.occupancy # 13
_atom_site.B_iso_or_equiv # 14
_atom_site.pdbx_formal_charge # 15
_atom_site.auth_seq_id # 16
_atom_site.auth_comp_id # 17
_atom_site.auth_asym_id # 18
_atom_site.auth_atom_id # 19
_atom_site.pdbx_PDB_model_num # 20

   0   1     2  3   4   5 6 7 8   9      10     11      12    13     141516     171819  20
ATOM   1     N  N   . PRO A 1 1   ? -21.348 -0.482  37.648  1.00 139.11 ? 1    PRO A N   1
ATOM   2     C  CA  . PRO A 1 1   ? -20.794 -1.467  36.686  1.00 152.26 ? 1    PRO A CA  1
ATOM   3     C  C   . PRO A 1 1   ? -19.682 -2.307  37.296  1.00 128.03 ? 1    PRO A C   1
ATOM   4     O  O   . PRO A 1 1   ? -18.771 -2.752  36.594  1.00 102.07 ? 1    PRO A O   1
ATOM   5     C  CB  . PRO A 1 1   ? -20.302 -0.682  35.484  1.00 112.52 ? 1    PRO A CB  1
ATOM   6     C  CG  . PRO A 1 1   ? -21.298 0.482   35.504  1.00 149.55 ? 1    PRO A CG  1
ATOM   7     C  CD  . PRO A 1 1   ? -21.457 0.838   37.004  1.00 109.96 ? 1    PRO A CD  1
***/

#if 0
void ChainReader2::IncFormatErrors()
	{
	m_CRGlobalLock.lock();
	++m_CRGlobalFormatErrors;
	m_CRGlobalLock.unlock();
	}

uint ChainReader2::GetCIFFieldIdx(const map<string, uint> &FieldToIdx, const string &Name)
	{
	map<string, uint>::const_iterator iter = FieldToIdx.find(Name);
	if (iter == FieldToIdx.end())
		{
		IncFormatErrors();
		Log("%s: CIF field %s not found\n",
		  m_CurrentFN.c_str(), Name.c_str());
		return UINT_MAX;
		}
	uint Idx = iter->second;
	return Idx;
	}

void ChainReader2::ChainsFromLines_CIF(const vector<string> &Lines,
  vector<PDBChain *> &Chains, const string &FallbackLabel)
	{
	Chains.clear();
	string TmpBaseLabel = FallbackLabel;
	const string &Line0 = Lines[0];
	if (StartsWith(Line0, "data_"))
		{
		vector<string> Fields;
		Split(Line0, Fields, '_');
		if (SIZE(Fields) == 2 && Fields[0] != "data")
			{
			TmpBaseLabel = Fields[1];
			if (TmpBaseLabel == "")
				TmpBaseLabel = FallbackLabel;
			}
		}
	const string &BaseLabel = TmpBaseLabel;

	string CurrentChainStr;
	PDBChain *Chain = 0;
	const uint N = SIZE(Lines);
	CIF_PARSER_STATE PS = PS_WaitingForLoop;
	vector<string> FieldList;
	vector<string> ATOMLines;
	for (uint i = 0; i < N; ++i)
		{
		const string &Line = Lines[i];
		if (!Line.empty() && Line[0] == '#')
			continue;
		if (PS == PS_Finished)
			break;

		switch (PS)
			{
		case PS_WaitingForLoop:
			{
			if (Line == "loop_")
				{
				PS = PS_AtLoop;
				continue;
				}
			break;
			}

		case PS_AtLoop:
			{
			if (StartsWith(Line, "_atom_site."))
				{
				PS = PS_InFieldList;
				string Name = Line;
				StripWhiteSpace(Name);
				FieldList.push_back(Name);
				}
			else
				PS = PS_WaitingForLoop;
			break;
			}

		case PS_InFieldList:
			{
			if (StartsWith(Line, "_atom_site."))
				{
				string Name = Line;
				StripWhiteSpace(Name);
				FieldList.push_back(Name);
				}
			else if (Line == "loop_")
				PS = PS_AtLoop;
			else if (StartsWith(Line, "ATOM ") ||
					 StartsWith(Line, "HETATM"))
				{
				PS = PS_InATOMs;
				ATOMLines.push_back(Line);
				}
			break;
			}

		case PS_InATOMs:
			{
			if (StartsWith(Line, "ATOM ") ||
				StartsWith(Line, "HETATM"))
				ATOMLines.push_back(Line);
			else
				PS = PS_Finished;
			break;
			}

		default:
			asserta(false);
			}
		}

	const uint FieldCount = SIZE(FieldList);
	map<string, uint> FieldToIdx;
	for (uint Idx = 0; Idx < FieldCount ; ++Idx)
		FieldToIdx[FieldList[Idx]] = Idx;

	uint Chain_FldIdx = GetCIFFieldIdx(FieldToIdx, "_atom_site.auth_asym_id");
	uint CA_FldIdx = GetCIFFieldIdx(FieldToIdx, "_atom_site.label_atom_id");
	uint X_FldIdx = GetCIFFieldIdx(FieldToIdx, "_atom_site.Cartn_x");
	uint Y_FldIdx = GetCIFFieldIdx(FieldToIdx, "_atom_site.Cartn_y");
	uint Z_FldIdx = GetCIFFieldIdx(FieldToIdx, "_atom_site.Cartn_z");
	uint aa_FldIdx = GetCIFFieldIdx(FieldToIdx, "_atom_site.label_comp_id");
	uint ModelNr_FldIdx = GetCIFFieldIdx(FieldToIdx, "_atom_site.pdbx_PDB_model_num");

	if (Chain_FldIdx == UINT_MAX) return;
	if (CA_FldIdx == UINT_MAX) return;
	if (X_FldIdx == UINT_MAX) return;
	if (Y_FldIdx == UINT_MAX) return;
	if (Z_FldIdx == UINT_MAX) return;
	if (aa_FldIdx == UINT_MAX) return;
	uint CurrentModelNr = UINT_MAX;

	const uint ATOMLineCount = SIZE(ATOMLines);
	for (uint i = 0; i < ATOMLineCount; ++i)
		{
		const string &Line = ATOMLines[i];
		vector<string> Fields;
		SplitWhite(Line, Fields);
		uint n = SIZE(Fields);
		if (n != FieldCount)
			{
			IncFormatErrors();
			Log("%s: Expected %u fields got %u in '%s'\n",
			  m_CurrentFN.c_str(), FieldCount, n, Line.c_str());
			if (Chain != 0)
				delete Chain;
			return;
			}

		const string &CA_Fld = Fields[CA_FldIdx];
		if (CA_Fld != "CA")
			continue;

		if (ModelNr_FldIdx != UINT_MAX)
			{
			uint ModelNr = (uint) atoi(Fields[ModelNr_FldIdx].c_str());
			if (CurrentModelNr != UINT_MAX && ModelNr != CurrentModelNr)
				break;
			CurrentModelNr = ModelNr;
			}

		const string &Chain_Fld = Fields[Chain_FldIdx];
		string ChainStr = Chain_Fld;
		if (ChainStr == "")
			ChainStr = "__";
		if (ChainStr != CurrentChainStr)
			{
			if (Chain != 0)
				Chains.push_back(Chain);
			Chain = new PDBChain;
			string Label = BaseLabel;
			ChainizeLabel(Label, ChainStr);
			Chain->m_Label = Label;
			CurrentChainStr = ChainStr;
			}

		const string &aa_Fld = Fields[aa_FldIdx];
		if (SIZE(aa_Fld) != 3)
			continue;
		float X = StrToFloatf(Fields[X_FldIdx]);
		float Y = StrToFloatf(Fields[Y_FldIdx]);
		float Z = StrToFloatf(Fields[Z_FldIdx]);
		char aa = GetOneFromThree(aa_Fld);

		Chain->m_Seq.push_back(aa);
		Chain->m_Xs.push_back(X);
		Chain->m_Ys.push_back(Y);
		Chain->m_Zs.push_back(Z);
		}
	if (Chain != 0)
		Chains.push_back(Chain);
	}
#endif // 0

void flat_chain_reader::ChainsFromLines_CIF(const vector<string> &Lines,
	vector<flat_chain_t *> &Chains, const string &FallbackLabel)
	{
	Chains.clear();

	string Entry;
	string Title;
	map<string, string> MolByEntity;
	map<string, vector<string> > RefsByEntity;
	ExtractCifMeta(Lines, FallbackLabel, Entry, Title, MolByEntity,
	  RefsByEntity);
	const string &BaseLabel = Entry;

	string CurrentChainStr;
	string CurrentEntityId;
	map<string, string> ChainToEntity;
	const uint N = SIZE(Lines);
	CIF_PARSER_STATE PS = PS_WaitingForLoop;
	vector<string> FieldList;
	vector<string> ATOMLines;
	vector<float> Xs, Ys, Zs;
	vector<char> aas;
	Xs.reserve(RESERVE_CHAIN_LENGTH);
	Ys.reserve(RESERVE_CHAIN_LENGTH);
	Zs.reserve(RESERVE_CHAIN_LENGTH);
	aas.reserve(RESERVE_CHAIN_LENGTH);
	for (uint i = 0; i < N; ++i)
		{
		const string &Line = Lines[i];
		if (!Line.empty() && Line[0] == '#')
			continue;
		if (PS == PS_Finished)
			break;

		switch (PS)
			{
		case PS_WaitingForLoop:
			{
			if (Line == "loop_")
				{
				PS = PS_AtLoop;
				continue;
				}
			break;
			}

		case PS_AtLoop:
			{
			if (StartsWith(Line, "_atom_site."))
				{
				PS = PS_InFieldList;
				string Name = Line;
				StripWhiteSpace(Name);
				FieldList.push_back(Name);
				}
			else
				PS = PS_WaitingForLoop;
			break;
			}

		case PS_InFieldList:
			{
			if (StartsWith(Line, "_atom_site."))
				{
				string Name = Line;
				StripWhiteSpace(Name);
				FieldList.push_back(Name);
				}
			else if (Line == "loop_")
				PS = PS_AtLoop;
			else if (StartsWith(Line, "ATOM ") ||
					 StartsWith(Line, "HETATM"))
				{
				PS = PS_InATOMs;
				ATOMLines.push_back(Line);
				}
			break;
			}

		case PS_InATOMs:
			{
			if (StartsWith(Line, "ATOM ") ||
				StartsWith(Line, "HETATM"))
				ATOMLines.push_back(Line);
			else
				PS = PS_Finished;
			break;
			}

		default:
			asserta(false);
			}
		}

	const uint FieldCount = SIZE(FieldList);
	map<string, uint> FieldToIdx;
	for (uint Idx = 0; Idx < FieldCount ; ++Idx)
		FieldToIdx[FieldList[Idx]] = Idx;

	uint Chain_FldIdx = GetCIFFieldIdx(FieldToIdx, "_atom_site.auth_asym_id");
	uint CA_FldIdx = GetCIFFieldIdx(FieldToIdx, "_atom_site.label_atom_id");
	uint X_FldIdx = GetCIFFieldIdx(FieldToIdx, "_atom_site.Cartn_x");
	uint Y_FldIdx = GetCIFFieldIdx(FieldToIdx, "_atom_site.Cartn_y");
	uint Z_FldIdx = GetCIFFieldIdx(FieldToIdx, "_atom_site.Cartn_z");
	uint aa_FldIdx = GetCIFFieldIdx(FieldToIdx, "_atom_site.label_comp_id");
	uint ModelNr_FldIdx = GetCIFFieldIdx(FieldToIdx, "_atom_site.pdbx_PDB_model_num");
	uint Entity_FldIdx = UINT_MAX;
	{
	map<string, uint>::const_iterator eit =
		FieldToIdx.find("_atom_site.label_entity_id");
	if (eit != FieldToIdx.end())
		Entity_FldIdx = eit->second;
	}

	if (Chain_FldIdx == UINT_MAX) return;
	if (CA_FldIdx == UINT_MAX) return;
	if (X_FldIdx == UINT_MAX) return;
	if (Y_FldIdx == UINT_MAX) return;
	if (Z_FldIdx == UINT_MAX) return;
	if (aa_FldIdx == UINT_MAX) return;
	uint CurrentModelNr = UINT_MAX;

	auto FinishChain = [&]()
		{
		if (aas.empty())
			return;
		string Label = BaseLabel;
		ChainizeLabel(Label, CurrentChainStr);
		string eid = CurrentEntityId;
		if (eid.empty())
			{
			map<string, string>::const_iterator it =
				ChainToEntity.find(CurrentChainStr);
			if (it != ChainToEntity.end())
				eid = it->second;
			}
		string db_ref = PickRefByEntity(RefsByEntity, eid);
		string molecule = PickMoleculeByEntity(MolByEntity, eid);
		AppendStructDescToLabel(Label, Entry, db_ref, molecule, Title);
		auto chain = flat_chain_t::newflat(Label, aas, Xs, Ys, Zs);
		Chains.push_back(chain);
		aas.clear();
		Xs.clear();
		Ys.clear();
		Zs.clear();
		};

	const uint ATOMLineCount = SIZE(ATOMLines);
	for (uint i = 0; i < ATOMLineCount; ++i)
		{
		const string &Line = ATOMLines[i];
		vector<string> Fields;
		SplitWhite(Line, Fields);
		uint n = SIZE(Fields);
		if (n != FieldCount)
			{
			IncFormatErrors();
			Log("%s: Expected %u fields got %u in '%s'\n",
			  m_CurrentFN.c_str(), FieldCount, n, Line.c_str());
			return;
			}

		const string &CA_Fld = Fields[CA_FldIdx];
		if (CA_Fld != "CA")
			continue;

		if (ModelNr_FldIdx != UINT_MAX)
			{
			uint ModelNr = (uint) atoi(Fields[ModelNr_FldIdx].c_str());
			if (CurrentModelNr != UINT_MAX && ModelNr != CurrentModelNr)
				break;
			CurrentModelNr = ModelNr;
			}

		const string &Chain_Fld = Fields[Chain_FldIdx];
		string ChainStr = Chain_Fld;
		if (ChainStr == "")
			ChainStr = "__";
		string entity_id;
		if (Entity_FldIdx != UINT_MAX)
			{
			entity_id = Fields[Entity_FldIdx];
			StripWhiteSpace(entity_id);
			if (!CifValuePresent(entity_id))
				entity_id.clear();
			}
		if (ChainStr != CurrentChainStr)
			{
			FinishChain();
			CurrentChainStr = ChainStr;
			CurrentEntityId = entity_id;
			if (!entity_id.empty() &&
			  ChainToEntity.find(ChainStr) == ChainToEntity.end())
				ChainToEntity[ChainStr] = entity_id;
			}
		else if (CurrentEntityId.empty() && !entity_id.empty())
			{
			CurrentEntityId = entity_id;
			if (ChainToEntity.find(ChainStr) == ChainToEntity.end())
				ChainToEntity[ChainStr] = entity_id;
			}

		const string &aa_Fld = Fields[aa_FldIdx];
		if (SIZE(aa_Fld) != 3)
			continue;
		float X = StrToFloatf(Fields[X_FldIdx]);
		float Y = StrToFloatf(Fields[Y_FldIdx]);
		float Z = StrToFloatf(Fields[Z_FldIdx]);
		char aa = GetOneFromThree(aa_Fld);

		aas.push_back(aa);
		Xs.push_back(X);
		Ys.push_back(Y);
		Zs.push_back(Z);
		}
	FinishChain();
	}

void flat_chain_reader::IncFormatErrors()
	{
	m_CRGlobalLock.lock();
	++m_CRGlobalFormatErrors;
	m_CRGlobalLock.unlock();
	}

uint flat_chain_reader::GetCIFFieldIdx(const map<string, uint> &FieldToIdx, const string &Name)
	{
	map<string, uint>::const_iterator iter = FieldToIdx.find(Name);
	if (iter == FieldToIdx.end())
		{
		IncFormatErrors();
		Log("%s: CIF field %s not found\n",
		  m_CurrentFN.c_str(), Name.c_str());
		return UINT_MAX;
		}
	uint Idx = iter->second;
	return Idx;
	}
