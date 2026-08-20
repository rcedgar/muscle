#pragma once

#include "linereader2.h"
#include <list>

class PDBFileScanner
	{
public:
	mutex m_Lock;
	string m_RootFileName;
	string m_CurrentFileName;
	list<string> m_PendingFiles;
	list<string> m_PendingDirs;
	LineReader2 m_LR;
	string m_Line;
	uint m_FileCount = 0;
	bool m_Trace = false;

public:
	void Open(const string &FN);
	bool GetNext(string &FN);

private:
	void ReadNextDir();
	bool PushFileOrDir(const string &FN);
	bool FileNameHasStructureExt(const string &FN) const;
	bool IsStructureExt(const string &Ext) const;
	};
