#include "muscle.h"
#include "tree.h"
#include "textfile.h"

#define TRACE 0

// Tokens in Newick files are:
//		( ) : , ;
//		string
//		'string'
//		"string"
//		[ comment ]
//
// We can't safely distinguish between identifiers and floating point
// numbers at the lexical level (because identifiers may be numeric,
// or start with digits), so both edge lengths and identifiers are
// returned as strings.

const char *Tree::NTTStr(NEWICK_TOKEN_TYPE NTT) const
	{
	switch (NTT)
		{
#define c(x)	case NTT_##x: return #x;
	c(Unknown)
	c(Lparen)
	c(Rparen)
	c(Colon)
	c(Comma)
	c(Semicolon)
	c(String)
	c(SingleQuotedString)
	c(DoubleQuotedString)
	c(Comment)
#undef c
		}
	return "??";
	}

NEWICK_TOKEN_TYPE Tree::GetToken(TextFile &File, char szToken[], unsigned uBytes) const
	{
// Skip leading white space
	File.SkipWhite();

	char c;
	File.GetCharX(c);

// In case a single-character token
	szToken[0] = c;
	szToken[1] = 0;

	unsigned uBytesCopied = 0;
	NEWICK_TOKEN_TYPE TT;
	switch (c)
		{
	case '(':
		return NTT_Lparen;

	case ')':
		return NTT_Rparen;

	case ':':
		return NTT_Colon;

	case ';':
		return NTT_Semicolon;

	case ',':
		return NTT_Comma;

	case '\'':
		TT = NTT_SingleQuotedString;
		File.GetCharX(c);
		break;

	case '"':
		TT = NTT_DoubleQuotedString;
		File.GetCharX(c);
		break;

	case '[':
		TT = NTT_Comment;
		break;

	default:
		TT = NTT_String;
		break;
		}

	for (;;)
		{
		if (TT != NTT_Comment)
			{
			if (uBytesCopied < uBytes - 2)
				{
				szToken[uBytesCopied++] = c;
				szToken[uBytesCopied] = 0;
				}
			else
				Die("Tree::GetToken: input buffer too small, token so far='%s'", szToken);
			}
		bool bEof = File.GetChar(c);
		if (bEof)
			return TT;

		switch (TT)
			{
		case NTT_String:
			if (0 != strchr("():;,", c))
				{
				File.PushBack(c);
				return NTT_String;
				}
			if (isspace(c))
				return NTT_String;
			break;

		case NTT_SingleQuotedString:
			if ('\'' == c)
				return NTT_String;
			break;

		case NTT_DoubleQuotedString:
			if ('"' == c)
				return NTT_String;
			break;

		case NTT_Comment:
			if (']' == c)
				return GetToken(File, szToken, uBytes);
			break;

		default:
			Die("Tree::GetToken, invalid TT=%u", TT);
			}
		}
	}

// NOTE: this hack must come after definition of Tree::GetToken.
#if	TRACE
#define GetToken	GetTokenVerbose
#endif

void Tree::FromFile(const string &FileName)
	{
	TextFile TF(FileName.c_str(), false);
	FromFile(TF);
	TF.Close();
	}

void Tree::FromFile(TextFile &File)
	{
// Assume rooted.
// If we discover that it is unrooted, will convert on the fly.
	CreateRooted();

	double dEdgeLength;
	bool bEdgeLength = GetGroupFromFile(File, 0, &dEdgeLength);

// Next token should be either ';' for rooted tree or ',' for unrooted.
	char szToken[16];
	NEWICK_TOKEN_TYPE NTT = GetToken(File, szToken, sizeof(szToken));

// If rooted, all done.
	if (NTT_Semicolon == NTT)
		{
		if (bEdgeLength)
			Log(" *** Warning *** edge length on root group in Newick file %s\n",
			  File.GetFileName());
		Validate();
		return;
		}

	if (NTT_Comma != NTT)
		Die("Tree::FromFile, expected ';' or ',', got '%s'", szToken);

	const unsigned uThirdNode = UnrootFromFile();
	bEdgeLength = GetGroupFromFile(File, uThirdNode, &dEdgeLength);
	if (bEdgeLength)
		SetEdgeLength(0, uThirdNode, dEdgeLength);
	Validate();
	}

// Return true if edge length for this group.
bool Tree::GetGroupFromFile(TextFile &File, unsigned uNodeIndex,
  double *ptrdEdgeLength)
	{
	char szToken[1024];
	NEWICK_TOKEN_TYPE NTT = GetToken(File, szToken, sizeof(szToken));

// Group is either leaf name or (left, right).
	if (NTT_String == NTT)
		{
		SetLeafName(uNodeIndex, szToken);
#if	TRACE
		Log("Group is leaf '%s'\n", szToken);
#endif
		}
	else if (NTT_Lparen == NTT)
		{
		const unsigned uLeft = AppendBranch(uNodeIndex);
		const unsigned uRight = uLeft + 1;

	// Left sub-group...
#if	TRACE
		Log("Got '(', group is compound, expect left sub-group\n");
#endif
		double dEdgeLength;
		bool bLeftLength = GetGroupFromFile(File, uLeft, &dEdgeLength);
#if	TRACE
		if (bLeftLength)
			Log("Edge length for left sub-group: %.3g\n", dEdgeLength);
		else
			Log("No edge length for left sub-group\n");
#endif
		if (bLeftLength)
			SetEdgeLength(uNodeIndex, uLeft, dEdgeLength);

	// ... then comma ...
#if	TRACE
		Log("Expect comma\n");
#endif
		NTT = GetToken(File, szToken, sizeof(szToken));
		if (NTT_Comma != NTT)
			Die("Tree::GetGroupFromFile, expected ',', got '%s'", szToken);

	// ...then right sub-group...
#if	TRACE
		Log("Expect right sub-group\n");
#endif
		bool bRightLength = GetGroupFromFile(File, uRight, &dEdgeLength);
		if (bRightLength)
			SetEdgeLength(uNodeIndex, uRight, dEdgeLength);

#if	TRACE
		if (bRightLength)
			Log("Edge length for right sub-group: %.3g\n", dEdgeLength);
		else
			Log("No edge length for right sub-group\n");
#endif

	// ... then closing parenthesis.
#if	TRACE
		Log("Expect closing parenthesis (or comma if > 2-ary)\n");
#endif
		NTT = GetToken(File, szToken, sizeof(szToken));
		if (NTT_Rparen == NTT)
			;
		else if (NTT_Comma == NTT)
			{
			File.PushBack(',');
			return false;
			}
		else
			Die("Tree::GetGroupFromFile, expected ')' or ',', got '%s'", szToken);
		}
	else
		Die("Tree::GetGroupFromFile, expected '(' or leaf name, got '%s'",
		  szToken);

// Group may optionally be followed by edge length.
	bool bEof = File.SkipWhiteX();
	if (bEof)
		return false;
	char c;
	File.GetCharX(c);
#if	TRACE
	Log("Character following group, could be colon, is '%c'\n", c);
#endif
	if (':' == c)
		{
		NTT = GetToken(File, szToken, sizeof(szToken));
		if (NTT_String != NTT)
			Die("Tree::GetGroupFromFile, expected edge length, got '%s'", szToken);
		*ptrdEdgeLength = atof(szToken);
		return true;
		}
	File.PushBack(c);
	return false;
	}
