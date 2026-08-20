#pragma once

enum FAN
	{
#define f(x)	FAN_##x,
#include "flat_featlist.h"
	FAN_COUNT
	};

const char *FAN2str(FAN fan);
FAN str2FAN(const char *s);
FAN str2FAN(const string &s);
bool is_quantized(FAN fan);

FAN parse_alpha_name(
	const string &alpha_name,
	uint &alpha_size);
