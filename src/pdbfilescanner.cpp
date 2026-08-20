#include "myutils.h"
#include "pdbfilescanner.h"

void PDBFileScanner::Open(const string &Path)
	{
	m_Trace = opt(trace_filescanner);
	if (m_Trace) Log("PDBFileScanner::Open(%s)\n", Path.c_str());
	m_RootFileName = Path;
	m_CurrentFileName.clear();
	m_PendingFiles.clear();
	m_PendingDirs.clear();
	bool Ok = PushFileOrDir(Path);
	if (!Ok)
		Die("Not found '%s'", Path.c_str());
	}

bool PDBFileScanner::IsStructureExt(const string &Ext) const
	{
	string LowerExt = Ext;
	ToLower(LowerExt);
	if (Ext == "bca" || Ext == "bcb")
		return true;

#define x(s)	if (Ext == #s || Ext == string(#s) + string(".gz")) return true;
	x(pdb)
	x(ent)
	x(cif)
	x(mmcif)
	x(cal)
	x(can)
#undef x
	return false;
	}

bool PDBFileScanner::FileNameHasStructureExt(const string &FN) const
	{
	string Ext;
	GetExtFromPathName(FN, Ext);
	return IsStructureExt(Ext);
	}

bool PDBFileScanner::PushFileOrDir(const string &Path)
	{
	if (m_Trace) Log("PDBFileScanner::PushFileOrDir(%s)\n", Path.c_str());
	if (Path.empty())
		return false;
	if (IsDirectory(Path))
		{
		m_PendingDirs.push_back(Path);
		return true;
		}
	else if (IsRegularFile(Path))
		{
		m_PendingFiles.push_back(Path);
		return true;
		}
	else
		Die("Invalid file/pathname '%s'", Path.c_str());
	return false;
	}

void PDBFileScanner::ReadNextDir()
	{
	if (m_PendingDirs.empty())
		return;
	const string Dir = m_PendingDirs.front();
	if (m_Trace) Log("PDBFileScanner::ReadNextDir(%s)\n", Dir.c_str());
	m_PendingDirs.pop_front();
	vector<string> FileNames;
	vector<bool> IsSubDirs;
	mylistdir(Dir, FileNames, IsSubDirs);
	for (uint i = 0; i < SIZE(FileNames); ++i)
		{
		const string &FN = FileNames[i];
		if (FN == "." || FN == "..")
			continue; 
		const string &Path = Dir + string("/") + FN;
		if (IsDirectory(Path))
			m_PendingDirs.push_back(Path);
		else if (FileNameHasStructureExt(FN))
			m_PendingFiles.push_back(Path);
		}
	}

bool PDBFileScanner::GetNext(string &FN)
	{
	FN.clear();
	m_Lock.lock();
	bool Ok = false;
	uint SanityCounter = 0;
	for (;;)
		{
		if (m_Trace) Log("PDBFileScanner::GetNext() SanityCounter=%u\n", SanityCounter);
		if (SanityCounter > 100) Die("PDBFileScanner::GetNext() SanityCounter=%u\n", SanityCounter);

		while (m_PendingFiles.empty() && !m_PendingDirs.empty())
			ReadNextDir();
		if (m_PendingFiles.empty())
			{
			Ok = false;
			break;
			}

		m_CurrentFileName = m_PendingFiles.front();
		m_PendingFiles.pop_front();

		string Ext;
		GetExtFromPathName(m_CurrentFileName, Ext);
		ToLower(Ext);
		if (m_Trace) Log("  m_CurrentFileName=%s Ext=%s\n", m_CurrentFileName.c_str(), Ext.c_str());
		if (FileNameHasStructureExt(m_CurrentFileName))
			{
			FN = m_CurrentFileName;
			++m_FileCount;
			Ok = true;
			break;
			}
		else if (Ext == "files")
			{
			vector<string> Paths;
			ReadLinesFromFile(m_CurrentFileName, Paths);
			for (uint i = 0; i < SIZE(Paths); ++i)
				{
				string Path = Paths[i];
				if (optset_filesdir)
					{
					string FilesDir = opt(filesdir);
					Dirize(FilesDir);
					Path = FilesDir + Path;
					}
				PushFileOrDir(Path);
				}
			}
		}
	m_Lock.unlock();
	return Ok;
	}

#if 0
void cmd_scan_files()
	{
	PDBFileScanner FS;
	FS.Open(g_Arg1);
	string FN;
	if (!optset_output)
		Die("-output required");
	FILE *fOut = CreateStdioFile(opt(output));

	time_t lastt = time(0);
	while (FS.GetNext(FN))
		{
		time_t now = time(0);
		if (now > lastt)
			{
			Progress("%u files\r", FS.m_FileCount);
			lastt = now;
			}

		fputs(FN.c_str(), fOut);
		fputc('\n', fOut);
		}
	Progress("%u files total\n", FS.m_FileCount);
	CloseStdioFile(fOut);
	}
#endif
