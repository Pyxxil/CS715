#include "Vector.h"

#include <windows.h>
#include <GL/GL.h>
#include <GL/GLU.h>
#include <GL/freeglut.h>

#include <iostream>

Vector::Vector(int _x, int _y, int _dx, int _dy)
	: x(_x), y(_y), dx(_dx), dy(_dy)
{}

void Vector::draw() const {
	glColor3f(255, 255, 255);

	glBegin(GL_LINES);

	glVertex2f(x, y);
	glVertex2f(x + (x - (x + dx)) * 0.25, y + (y - (y + dy)) * 0.25);

	glEnd();
}