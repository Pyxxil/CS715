#include <windows.h>
#include <GL/GL.h>
#include <GL/GLU.h>
#include <GL/freeglut.h>
//
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
 *   Update the forces on each particle
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
 *  Calculate the derivatives
 *  dp/dt = v
 *  dv/dt = f / m
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

enum METHOD {
    EULER,
    MIDPOINT,
    RUNGE_KUTTA,
};

/*
 * Perform one step of the solver
 */
void UpdateParticles(
    PARTICLE *p, int np,
    PARTICLEPHYS phys,
    PARTICLESPRING *s, int ns,
    double dt, METHOD method)
{
    int i;
    PARTICLE *ptmp, *ptmp1, *ptmp2, *ptmp3;
    PARTICLEDERIVATIVES *deriv;

    deriv = (PARTICLEDERIVATIVES *)malloc(np * sizeof(PARTICLEDERIVATIVES));

    switch (method) {
    case EULER:
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
    case MIDPOINT:
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
    case RUNGE_KUTTA:
        CalculateForces(p, np, phys, s, ns);
        CalculateDerivatives(p, np, deriv);
        ptmp = (PARTICLE *)malloc(np * sizeof(PARTICLE));
        ptmp1 = (PARTICLE *)malloc(np * sizeof(PARTICLE));
        ptmp2 = (PARTICLE *)malloc(np * sizeof(PARTICLE));
        ptmp3 = (PARTICLE *)malloc(np * sizeof(PARTICLE));
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
            ptmp1[i] = p[i];
            ptmp1[i].p.x += deriv[i].dpdt.x * dt + ptmp[i].p.x / 2;
            ptmp1[i].p.y += deriv[i].dpdt.y * dt + ptmp[i].p.y / 2;
            ptmp1[i].p.z += deriv[i].dpdt.z * dt + ptmp[i].p.z / 2;
            ptmp1[i].p.x += deriv[i].dvdt.x * dt / 2;
            ptmp1[i].p.y += deriv[i].dvdt.y * dt / 2;
            ptmp1[i].p.z += deriv[i].dvdt.z * dt / 2;
        }

        CalculateForces(ptmp1, np, phys, s, ns);
        CalculateDerivatives(ptmp1, np, deriv);
        for (i = 0; i < np; i++) {
            ptmp2[i] = p[i];
            ptmp2[i].p.x += deriv[i].dpdt.x * dt + ptmp1[i].p.x / 2;
            ptmp2[i].p.y += deriv[i].dpdt.y * dt + ptmp1[i].p.y / 2;
            ptmp2[i].p.z += deriv[i].dpdt.z * dt + ptmp1[i].p.z / 2;
            ptmp2[i].p.x += deriv[i].dvdt.x * dt / 2;
            ptmp2[i].p.y += deriv[i].dvdt.y * dt / 2;
            ptmp2[i].p.z += deriv[i].dvdt.z * dt / 2;
        }

        CalculateForces(ptmp2, np, phys, s, ns);
        CalculateDerivatives(ptmp2, np, deriv);
        for (i = 0; i < np; i++) {
            ptmp3[i] = p[i];
            ptmp3[i].p.x += deriv[i].dpdt.x * dt + ptmp2[i].p.x / 2;
            ptmp3[i].p.y += deriv[i].dpdt.y * dt + ptmp2[i].p.y / 2;
            ptmp3[i].p.z += deriv[i].dpdt.z * dt + ptmp2[i].p.z / 2;
            ptmp3[i].p.x += deriv[i].dvdt.x * dt;
            ptmp3[i].p.y += deriv[i].dvdt.y * dt;
            ptmp3[i].p.z += deriv[i].dvdt.z * dt;
        }

        for (i = 0; i < np; i++) {
            p[i].p.x += (ptmp[i].p.x + ptmp1[i].p.x * 2 + ptmp2[i].p.x * 2 + ptmp3[i].p.x) / 6;
            p[i].p.y += (ptmp[i].p.y + ptmp1[i].p.y * 2 + ptmp2[i].p.y * 2 + ptmp3[i].p.y) / 6;
            p[i].p.z += (ptmp[i].p.z + ptmp1[i].p.z * 2 + ptmp2[i].p.z * 2 + ptmp3[i].p.z) / 6;
        }

        CalculateForces(p, np, phys, s, ns);
        CalculateDerivatives(p, np, deriv);
        for (i = 0; i < np; i++) {
            p[i] = p[i];
            p[i].p.x += deriv[i].dpdt.x;
            p[i].p.y += deriv[i].dpdt.y;
            p[i].p.z += deriv[i].dpdt.z;
            p[i].p.x += deriv[i].dvdt.x;
            p[i].p.y += deriv[i].dvdt.y;
            p[i].p.z += deriv[i].dvdt.z;
        }

        free(ptmp);
        free(ptmp1);
        free(ptmp2);
        free(ptmp3);
        break;
    }

    free(deriv);
}

constexpr double DT = 1.f / 64;
int nparticles = 0;
PARTICLE *particles;
int nsprings = 0;
PARTICLESPRING *springs;
PARTICLEPHYS physical;

/*!
 * Display the entire spring system.
 */
void display() {
    glClearColor(29 / 255.f, 42 / 255.f, 60 / 255.f, 1);
    glClear(GL_COLOR_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluOrtho2D(0, width * 2, 0, height * 2);

    glColor3f(255, 128, 0);
    glLineWidth(2);

    glLineWidth(10);

    glColor3f(0, 128, 255);
    glBegin(GL_QUADS);
    for (int i = 0; i < nparticles; i++) {
        /* Display a particle at particles[i].p */
        glVertex2f(particles[i].p.x + 15, particles[i].p.y + 15);
        glVertex2f(particles[i].p.x - 15, particles[i].p.y + 15);
        glVertex2f(particles[i].p.x - 15, particles[i].p.y - 15);
        glVertex2f(particles[i].p.x + 15, particles[i].p.y - 15);
    }
    glEnd();

    glLineWidth(2);
    glColor3f(255, 128, 0);
    for (int i = 0; i < nsprings; i++) {
        const int p1 = springs[i].from;
        const int p2 = springs[i].to;

        glBegin(GL_LINES);
        /* Display a spring between point p1 and p2 */
        glVertex2f(particles[p1].p.x, particles[p1].p.y);
        glVertex2f(particles[p2].p.x, particles[p2].p.y);
        glEnd();
    }

    glLineWidth(1);

    glFlush();
}

/*!
 * Move the masses on every redraw.
 */
void moveParticle(int) {
    UpdateParticles(particles, nparticles, physical, springs, nsprings, DT, MIDPOINT);
    glutPostRedisplay();
    glutTimerFunc(10, moveParticle, 0);
}

/*
 * Set up the particle system.
 * Initialise all variables to default values
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
        springs[i].dampingconstant = 1; // Critically damped
        springs[i].restlength = 50;
    }

    physical.gravitational = -9.8;
    physical.viscousdrag = 0.0000000000000000001;
}

/*
 *
 */
void InitialiseSystem(void)
{
    int i;

    SetupParticles(25, 72);

    // Initialise position of particles
    for (i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            particles[i * 5 + j].p.x = 800 * j + 50;
            particles[i * 5 + j].p.y = 500 * i + 50;
            particles[i * 5 + j].p.z = 200 * j + 50;
        }
    }

    // Randomly fix 5 particles
    for (i = 0; i < 5; i++) {
        particles[rand() % nparticles].fixed = true;
    }

    // Implement all of the springs with their from and to masses, including the diagonal springs.

    for (i = 0; i < 5; i++) {
        for (int j = 0; j < 4; j++) {
            springs[i * 4 + j].from = i * 5 + j;
            springs[i * 4 + j].to = i * 5 + j + 1;
        }
    }

    for (i = 0; i < 4; i++) {
        for (int j = 0; j < 5; j++) {
            springs[i * 5 + 20 + j].from = i * 5 + j;
            springs[i * 5 + 20 + j].to = i * 5 + 5 + j;
        }
    }

    for (i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            springs[(i + 10) * 4 + j].from = i * 5 + j;
            springs[(i + 10) * 4 + j].to = i * 5 + j + 6;

            springs[(i + 14) * 4 + j].from = i * 5 + j + 1;
            springs[(i + 14) * 4 + j].to = i * 5 + j + 5;
        }
    }
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutCreateWindow("");
    srand(time(NULL));
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutFullScreen();
    InitialiseSystem();
    glutInitWindowSize(100, 100);
    glutInitWindowPosition(50, 50);
    glutDisplayFunc(display);
    glutTimerFunc(10, moveParticle, 1);
    glutMainLoop();
    return 0;
}