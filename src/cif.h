#pragma once

enum CIF_PARSER_STATE
	{
	PS_WaitingForLoop,
	PS_AtLoop,
	PS_InFieldList,
	PS_InATOMs,
	PS_Finished,
	};
