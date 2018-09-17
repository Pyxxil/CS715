#include "VectorField.h"

#include <iostream>

#include "Vector.h"

static constexpr double e = 2.718281;
static constexpr double pi = 3.141592654;

VectorField::VectorField(int width, int height) {
	const int a = 10, b = 10;
	double angle;

	for (int i = 0; i < 128; i++) {
		angle = 0.5 * i;

		double x = (a + b * angle) * std::cos(angle);
		double y = (a + b * angle) * std::sin(angle);

		mVectors.emplace_back(
			Vector(x, y,
				(a + b * ((i + 1) * 0.5) * std::cos((i + 1) * 0.5)) - x,
				(a + b * ((i + 1) * 0.5) * std::sin((i + 1) * 0.5)) - y
			)
		);
	}
}

void VectorField::draw() const {
	for (const auto &vec : mVectors) {
		vec.draw();
	}
}

static double distance(int x1, int x2, int y1, int y2) {
	const auto dx = x1 - x2;
	const auto dy = y1 - y2;
	return std::sqrt(dx * dx + dy * dy);
}

std::pair<int, int> VectorField::forces(int x, int y) const {
	std::pair<int, int> force{ 0, 0 };
	auto min = std::numeric_limits<int>::max();
	for (auto &vec : mVectors) {
		const auto dist = distance(x, vec.X(), y, vec.Y());
		if (dist < min) {
			force = { vec.DX(), vec.DY() };
			min = dist;
		}
	}

	std::cout << "Force is " << force.first << ", " << force.second << '\n';
	return force;
}