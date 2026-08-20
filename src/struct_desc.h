#pragma once

#include <map>
#include <string>
#include <vector>

bool CifValuePresent(const std::string &value);
void CleanStructText(std::string &s);
std::string CleanStructTextCopy(const std::string &s);

std::string FormatDbRef(const std::string &db_name, const std::string &accession);
std::string PreferUnp(const std::vector<std::string> &refs);

void AppendStructDescToLabel(std::string &Label, const std::string &Entry,
	const std::string &DbRef, const std::string &Molecule, const std::string &Title);

std::string PickRefByEntity(
	const std::map<std::string, std::vector<std::string> > &by_entity,
	const std::string &entity_id);
std::string PickMoleculeByEntity(
	const std::map<std::string, std::string> &mol_by_entity,
	const std::string &entity_id);
std::string PickMoleculeByChain(
	const std::map<std::string, std::string> &mol_by_chain,
	const std::string &chain);

void ExtractCifMeta(const std::vector<std::string> &Lines,
	const std::string &FallbackLabel,
	std::string &Entry, std::string &Title,
	std::map<std::string, std::string> &MolByEntity,
	std::map<std::string, std::vector<std::string> > &RefsByEntity);

void ExtractPdbMeta(const std::vector<std::string> &Lines, std::string &Entry,
	std::string &Title,
	std::map<std::string, std::string> &MolByChain,
	std::map<std::string, std::vector<std::string> > &RefsByChain);
