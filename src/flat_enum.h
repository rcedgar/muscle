#pragma once

enum FE
	{
#define x(name)	FE_##name,
#include "flat_type_names.h"
	FE_N
	};

static inline const char *FE2str(FE fe)
	{
	switch (fe)
		{
#define x(name)	case FE_##name : return #name;
#include "flat_type_names.h"
		}
	return "FE_ERROR";
	}