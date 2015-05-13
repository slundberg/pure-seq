#ifndef CONTIGSHG38_H
#define CONTIGSHG38_H

#include <array>
#include <string>

class ContigsHg38 {
public:
	static const int count;
	static const std::array<std::string, 455> names;
	static const int64_t sizes[];
};

#endif