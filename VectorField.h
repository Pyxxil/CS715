#pragma once

#include <vector>

#include "Vector.h"

class VectorField {
public:
	VectorField(int width, int height);

	void draw() const;

	std::pair<int, int> forces(int x, int y) const;

private:
	std::vector<Vector> mVectors;
};