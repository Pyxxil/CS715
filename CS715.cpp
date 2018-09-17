#include <windows.h>
#include <GL/GL.h>
#include <GL/GLU.h>
#include <GL/freeglut.h>

#include <iostream>
#include <ctime>
#include <vector>

static constexpr int width = 1920;
static constexpr int height = 1080;

typedef struct {
	float x;
	float y;
	float z;
} XYZ;

typedef struct {
	double m;   /* Mass                          */
	XYZ p;      /* Position                      */
	XYZ v;      /* Velocity                      */
	XYZ f;      /* Force                         */
	int fixed;  /* Fixed point or free to move   */
} PARTICLE;

typedef struct {
	XYZ dpdt;
	XYZ dvdt;
} PARTICLEDERIVATIVES;

typedef struct {
	double gravitational;
	double viscousdrag;
} PARTICLEPHYS;

typedef struct {
	int from;
	int to;
	double springconstant;
	double dampingconstant;
	double restlength;
} PARTICLESPRING;

/*
   Update the forces on each particle
*/
void CalculateForces(
	PARTICLE *p, int np,
	PARTICLEPHYS phys,
	PARTICLESPRING *s, int ns)
{
	int i, p1, p2;
	XYZ down = { 0.0,0.5,-1.0 };
	XYZ zero = { 0.0,0.0,0.0 };
	XYZ f;
	double len, dx, dy, dz;

	for (i = 0; i < np; i++) {
		p[i].f = zero;
		if (p[i].fixed)
			continue;

		/* Gravitation */
		p[i].f.x += phys.gravitational * p[i].m * down.x;
		p[i].f.y += phys.gravitational * p[i].m * down.y;
		p[i].f.z += phys.gravitational * p[i].m * down.z;

		/* Viscous drag */
		p[i].f.x -= phys.viscousdrag * p[i].v.x;
		p[i].f.y -= phys.viscousdrag * p[i].v.y;
		p[i].f.z -= phys.viscousdrag * p[i].v.z;
	}

	/* Handle the spring interaction */
	for (i = 0; i < ns; i++) {
		p1 = s[i].from;
		p2 = s[i].to;
		dx = p[p1].p.x - p[p2].p.x;
		dy = p[p1].p.y - p[p2].p.y;
		dz = p[p1].p.z - p[p2].p.z;
		len = sqrt(dx*dx + dy * dy + dz * dz);
		f.x = s[i].springconstant  * (len - s[i].restlength);
		f.x += s[i].dampingconstant * (p[p1].v.x - p[p2].v.x) * dx / len;
		f.x *= -dx / len;
		f.y = s[i].springconstant  * (len - s[i].restlength);
		f.y += s[i].dampingconstant * (p[p1].v.y - p[p2].v.y) * dy / len;
		f.y *= -dy / len;
		f.z = s[i].springconstant  * (len - s[i].restlength);
		f.z += s[i].dampingconstant * (p[p1].v.z - p[p2].v.z) * dz / len;
		f.z *= -dz / len;
		if (!p[p1].fixed) {
			p[p1].f.x += f.x;
			p[p1].f.y += f.y;
			p[p1].f.z += f.z;
		}
		if (!p[p2].fixed) {
			p[p2].f.x -= f.x;
			p[p2].f.y -= f.y;
			p[p2].f.z -= f.z;
		}
	}
}

/*
   Calculate the derivatives
   dp/dt = v
   dv/dt = f / m
*/
void CalculateDerivatives(
	PARTICLE *p, int np,
	PARTICLEDERIVATIVES *deriv)
{
	int i;

	for (i = 0; i < np; i++) {
		deriv[i].dpdt.x = p[i].v.x;
		deriv[i].dpdt.y = p[i].v.y;
		deriv[i].dpdt.z = p[i].v.z;
		deriv[i].dvdt.x = p[i].f.x / p[i].m;
		deriv[i].dvdt.y = p[i].f.y / p[i].m;
		deriv[i].dvdt.z = p[i].f.z / p[i].m;
	}
}

/*
   Perform one step of the solver
*/
void UpdateParticles(
	PARTICLE *p, int np,
	PARTICLEPHYS phys,
	PARTICLESPRING *s, int ns,
	double dt, int method)
{
	int i;
	PARTICLE *ptmp;
	PARTICLEDERIVATIVES *deriv;

	deriv = (PARTICLEDERIVATIVES *)malloc(np * sizeof(PARTICLEDERIVATIVES));

	switch (method) {
	case 0:                                   /* Euler */
		CalculateForces(p, np, phys, s, ns);
		CalculateDerivatives(p, np, deriv);
		for (i = 0; i < np; i++) {
			p[i].p.x += deriv[i].dpdt.x * dt;
			p[i].p.y += deriv[i].dpdt.y * dt;
			p[i].p.z += deriv[i].dpdt.z * dt;
			p[i].v.x += deriv[i].dvdt.x * dt;
			p[i].v.y += deriv[i].dvdt.y * dt;
			p[i].v.z += deriv[i].dvdt.z * dt;
		}
		break;
	case 1:                                   /* Midpoint */
		CalculateForces(p, np, phys, s, ns);
		CalculateDerivatives(p, np, deriv);
		ptmp = (PARTICLE *)malloc(np * sizeof(PARTICLE));
		for (i = 0; i < np; i++) {
			ptmp[i] = p[i];
			ptmp[i].p.x += deriv[i].dpdt.x * dt / 2;
			ptmp[i].p.y += deriv[i].dpdt.y * dt / 2;
			ptmp[i].p.z += deriv[i].dpdt.z * dt / 2;
			ptmp[i].p.x += deriv[i].dvdt.x * dt / 2;
			ptmp[i].p.y += deriv[i].dvdt.y * dt / 2;
			ptmp[i].p.z += deriv[i].dvdt.z * dt / 2;
		}
		CalculateForces(ptmp, np, phys, s, ns);
		CalculateDerivatives(ptmp, np, deriv);
		for (i = 0; i < np; i++) {
			p[i].p.x += deriv[i].dpdt.x * dt;
			p[i].p.y += deriv[i].dpdt.y * dt;
			p[i].p.z += deriv[i].dpdt.z * dt;
			p[i].v.x += deriv[i].dvdt.x * dt;
			p[i].v.y += deriv[i].dvdt.y * dt;
			p[i].v.z += deriv[i].dvdt.z * dt;
		}
		free(ptmp);
		break;
	}

	free(deriv);
}

/*
   A 1st order 1D DE solver.
   Assumes the function is not time dependent.
   Parameters
	  h      - step size
	  y0     - last value
	  method - algorithm to use [0,5]
	  fcn    - evaluate the derivative of the field
*/
double Solver1D(double h, double y0, int method, double(*fcn)(double))
{
	double ynew;
	double k1, k2, k3, k4, k5, k6;

	switch (method) {
	case 0:                          /* Euler method */
		k1 = h * (*fcn)(y0);
		ynew = y0 + k1;
		break;
	case 1:                          /* Modified Euler */
		k1 = h * (*fcn)(y0);
		k2 = h * (*fcn)(y0 + k1);
		ynew = y0 + (k1 + k2) / 2;
		break;
	case 2:                          /* Heuns method */
		k1 = h * (*fcn)(y0);
		k2 = h * (*fcn)(y0 + 2 * k1 / 3);
		ynew = y0 + k1 / 4 + 3 * k2 / 4;
		break;
	case 3:                          /* Midpoint */
		k1 = h * (*fcn)(y0);
		k2 = h * (*fcn)(y0 + k1 / 2);
		ynew = y0 + k2;
		break;
	case 4:                          /* 4'th order Runge-kutta */
		k1 = h * (*fcn)(y0);
		k2 = h * (*fcn)(y0 + k1 / 2);
		k3 = h * (*fcn)(y0 + k2 / 2);
		k4 = h * (*fcn)(y0 + k3);
		ynew = y0 + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6;
		break;
	case 5:                          /* England 4'th order, six stage */
		k1 = h * (*fcn)(y0);
		k2 = h * (*fcn)(y0 + k1 / 2);
		k3 = h * (*fcn)(y0 + (k1 + k2) / 4);
		k4 = h * (*fcn)(y0 - k2 + 2 * k3);
		k5 = h * (*fcn)(y0 + (7 * k1 + 10 * k2 + k4) / 27);
		k6 = h * (*fcn)(y0 + (28 * k1 - 125 * k2 + 546 * k3 + 54 * k4 - 378 * k5) / 625);
		ynew = y0 + k1 / 6 + 4 * k3 / 6 + k4 / 6;
		break;
	}

	return ynew;
}

std::pair<GLfloat, GLfloat> particle{ 50, 500 };

std::vector<std::pair<GLfloat, GLfloat>> positions{};

constexpr double DEGREE_TO_RADIANS = 3.141592654 / 180.f;

constexpr double DT = 0.03125;

double distance(int x1, int x2, int y1, int y2) {
	const auto dx = x1 - x2;
	const auto dy = y1 - y2;
	return std::sqrt(dx * dx + dy * dy);
}

auto radians(int degrees) {
	return degrees * DEGREE_TO_RADIANS;
}

int nparticles = 0;
PARTICLE *particles;
int nsprings = 0;
PARTICLESPRING *springs;
PARTICLEPHYS physical;

void display() {
	glClearColor(29 / 255.f, 42 / 255.f, 60 / 255.f, 1);
	glClear(GL_COLOR_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(-width, width, height, -height);

	static const double h = std::cos(20) * 10;

	glColor3f(255, 128, 0);
	glLineWidth(2);
	//for (int radius = 1; radius < width * 2; radius += width / 16) {
	//	for (int i = 0; i < 360; i++) {
	//		glBegin(GL_LINES);
	//		const auto x0 = cos(radians(i)) * radius;
	//		const auto y0 = sin(radians(i)) * radius;
	//		const auto x1 = cos(radians(i + 1)) * radius;
	//		const auto y1 = sin(radians(i + 1)) * radius;
	//		glVertex2f(x0, y0);
	//		glVertex2f(x1, y1);
	//		glEnd();
	//	}
	//}

	glLineWidth(10);

	glColor3f(0, 128, 255);
	for (int i = 0; i < nparticles; i++) {
		/* Display a particle at particles[i].p */
		glBegin(GL_QUADS);
		glVertex2f(particles[i].p.x + 15, particles[i].p.y + 15);
		glVertex2f(particles[i].p.x - 15, particles[i].p.y + 15);
		glVertex2f(particles[i].p.x - 15, particles[i].p.y - 15);
		glVertex2f(particles[i].p.x + 15, particles[i].p.y - 15);
		glEnd();
	}

	glLineWidth(2);
	glColor3f(255, 128, 0);
	glBegin(GL_LINES);
	for (int i = 0; i < nsprings; i++) {
		const int p1 = springs[i].from;
		const int p2 = springs[i].to;
		/* Display a spring between point p1 and p2 */
		glVertex2f(particles[p1].p.x, particles[p1].p.y);
		glVertex2f(particles[p2].p.x, particles[p2].p.y);
	}
	glEnd();

	glColor3f(0, 128, 255);
	glLineWidth(5);
	glBegin(GL_LINE_STRIP);
	//for (const auto &pos : positions) {
	//	glVertex2f(pos.first, pos.second);
	//}
	glEnd();
	glLineWidth(1);

	glFlush();
}

void moveParticle(int) {
	positions.push_back(particle);
	particle = { particle.first + DT * particle.second, particle.second - DT * particle.first };

	if ((particle.first <= -width || particle.first >= width) && (particle.second <= -height || particle.second >= height)) {
		positions = {};
		particle = { 50, 500 };
	}

	UpdateParticles(particles, nparticles, physical, springs, nsprings, DT, 0);

	glutPostRedisplay();
	glutTimerFunc(10, moveParticle, 0);
}

/*
   Set up the particle system.
   Initialise all variables to default values
*/
void SetupParticles(int np, int ns)
{
	int i;

	nparticles = np;
	nsprings = ns;

	if (particles != NULL)
		free(particles);
	particles = (PARTICLE *)malloc(nparticles * sizeof(PARTICLE));
	if (springs != NULL)
		free(springs);
	springs = (PARTICLESPRING *)malloc(nsprings * sizeof(PARTICLESPRING));

	for (i = 0; i < np; i++) {
		particles[i].m = 1;
		particles[i].fixed = FALSE;
		particles[i].v.x = 0;
		particles[i].v.y = 0;
		particles[i].v.z = 0;
	}
	for (i = 0; i < ns; i++) {
		springs[i].springconstant = 10;
		springs[i].dampingconstant = 4;
		springs[i].restlength = 200;
	}

	physical.gravitational = 9.8;
	physical.viscousdrag = 0.0000000000000000001;
}

/*
   Create 5 particles
   The particles are placed randomly on the interval -5 -> 5
*/
void InitialiseSystem(void)
{
	int i;

	SetupParticles(2, 1);
	particles[0].p = { 0, 0, 0 };
	particles[0].fixed = true;
	particles[1].p = { 0, 600, 0 };
	//particles[2].p = { 300, 0, 0 };
	//particles[2].fixed = true;

	/* Random positions */
	//for (i = 0; i < 6; i++) {
	//	particles[i].p.x = rand() % (width * 2) - width;
	//	particles[i].p.y = rand() % (height * 2) - height;
	//	particles[i].p.z = 10 * (rand() % 1000 - 500) / 1000.0;
	//	std::cout << "Particle " << i << " at (" << particles[i].p.x << ", " << particles[i].p.y << ", " << particles[i].p.z << ")\n";
	//}

	/* Edges */
	springs[0].from = 0; springs[0].to = 1;
	//springs[1].from = 2; springs[1].to = 1;
	//springs[1].from = 1; springs[1].to = 2;
	//springs[2].from = 0; springs[2].to = 3;
	//springs[3].from = 0; springs[3].to = 4;
	//springs[4].from = 0; springs[4].to = 5;
	//springs[5].from = 1; springs[5].to = 2;
	//springs[6].from = 1; springs[6].to = 3;
	//springs[7].from = 1; springs[7].to = 4;
	//springs[8].from = 1; springs[8].to = 5;
	//springs[9].from = 2; springs[9].to = 3;
	//springs[10].from = 2; springs[10].to = 4;
	//springs[11].from = 2; springs[11].to = 5;
	//springs[12].from = 3; springs[12].to = 4;
	//springs[13].from = 3; springs[13].to = 5;
	//springs[14].from = 4; springs[14].to = 5;
}

int main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutCreateWindow("");
	srand(time(NULL));
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutFullScreen();
	InitialiseSystem();
	glutInitWindowSize(width, height);
	glutInitWindowPosition(50, 50);
	glutDisplayFunc(display);
	glutTimerFunc(10, moveParticle, 1);
	glutMainLoop();
	return 0;
}