#pragma once

class Vector {
public:
	Vector(int _x, int _y, int _dx, int _dy);

	void draw() const;

	int X() const { return x; }
	int Y() const { return y; }
	int DX() const { return dx; }
	int DY() const { return dy; }

private:
	int x;
	int y;
	int dx;
	int dy;
};
