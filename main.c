// This is a very simple implement of 2D SPH(smoothed particle hydrodynamics) in c
// This privide real-time simulation using OpenGL
// by Btli

#include <GL/glut.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI  3.141592653589793
#define MAX_PT  2500

//fluid parameters
const double H = 16.0;    // kernel radius
const double Acc_gravity = 9.8 * 12000;
const double Rest_density = 1000.0; //rest density, rho_0
const double Gas_const = 2000.0; //const k of equation of state, P = k(rho-rho_0)
const double Mass = 65.0;   //assume all particles have the same mass
const double Viscosity = 250.0; //viscosity constant
const double Dt = 0.00008; //integration time step

// boundary parameters
const double Boundary_thickness = 16.0;
const double Bound_damping = -0.5;  // rebound at half velocity 

// interaction
const int Max_particles = 2500;
const int Large_cube_particles = 494;
const int Small_cube_particles = 100;

//rendering projection parameters
const int Window_width = 800;
const int Window_height = 600;
const double View_width = 1.5 * 800.0;
const double View_height = 1.5 * 800.0;

// particle parameters
struct Particle {
    double x, y;    //position
    double u, v;        //velocity
    double fx, fy;      //forces
    double rho;         //density
    double pressure;    //pressure
} particle[MAX_PT];

int num = 0;    // active number of particles


//optimization parameters
double H_sq = 0;    //square of kernel radius
double Ker_ploy6 = 0;   //coefficient constant of kernel function for density
double Ker_spiky_grad = 0;   //coefficient constant of the gradient of kernel function for pressure
double ker_visc_lap = 0;    //coefficient constant of the laplacian kernel function for viscosity

void initialize_opt()
{
    H_sq = H * H;
    Ker_ploy6 = 315.0 / (64.0 * PI * pow(H, 9.0));
    Ker_spiky_grad = -45.0 / (PI * pow(H, 6.0));
    ker_visc_lap = 45.0 / (PI * pow(H, 6.0));
}

//initialize the particle systems
void init_SPH()
{
    srand(time(0)); // random number
    for (double y = Boundary_thickness; y < View_height - 2.0 * Boundary_thickness; y += H) {
        for (double x = View_width / 4; x < View_width / 2; x += H) {
            //add a little jitter (wait)
            if (num < Large_cube_particles) {
                double jitter = (rand() % 100 - 50) / 50.0; //random number between (-1, 1)
                particle[num].x = x + jitter;
                particle[num].y = y;
                particle[num].u = 0.0;
                particle[num].v = 0.0;
                particle[num].fx = 0.0;
                particle[num].fy = 0.0;
                particle[num].rho = 0;
                particle[num].pressure = 0;
                num++;
            }
            else {
                return;
            }
        }
    }
}

//calculate the density and pressure of each particle
void compute_density_and_pressure()
{
    double r2;
    for (int i = 0; i < num; i++) {
        particle[i].rho = 0.0;
        for (int j = 0; j < num; j++) {
            // compute the square distance between two particles
            r2 = pow(particle[i].x - particle[j].x, 2) + pow(particle[i].y - particle[j].y, 2);

            if (r2 < H_sq)  // inside the support zone
            {
                particle[i].rho += Mass * Ker_ploy6 * pow(H_sq - r2, 3.0);    /*a kind of weighted summation*/
            }
        }

        //using the equation of state of ideal gas
        particle[i].pressure = Gas_const * (particle[i].rho - Rest_density);
    }
}

// calculate the total force of each particle, F = F_pressure + F_viscosity + F_gravity
void compute_forces()
{
    double rx, ry, r, force;
    for (int i = 0; i < num; i++) {
        particle[i].fx = 0, particle[i].fy = 0;
        for (int j = 0; j < num; j++) {
            if (i == j) // force can only added by other particles
            {
                continue;
            }
            //calculate the distance between two particles
            rx = particle[j].x - particle[i].x;
            ry = particle[j].y - particle[i].y;
            r = pow(rx * rx + ry * ry, 0.5);

            if (r < H)  //inside the support zone
            {
                //compute the pressure force
                force = -1.0 * Mass * (particle[i].pressure + particle[j].pressure) / (2.0 * particle[j].rho) \
                    * Ker_spiky_grad * pow(H - r, 2.0);
                particle[i].fx += force * rx / r;
                particle[i].fy += force * ry / r;

                //compute the viscosity force
                force = Viscosity * Mass / particle[j].rho * ker_visc_lap * (H - r);
                particle[i].fx += force * (particle[j].u - particle[i].u);
                particle[i].fy += force * (particle[j].v - particle[i].v);
            }
        }
        // add the gravity, downwards
        particle[i].fy -= Acc_gravity * particle[i].rho;
    }
}

//update the velocity and position based on the force and density, a = F / rho
void integrate()
{
    double u, v;
    for (int i = 0; i < num; i++) {
        u = particle[i].u;
        v = particle[i].v;

        //update the velocity, using Euler's method
        particle[i].u += Dt * particle[i].fx / particle[i].rho;
        particle[i].v += Dt * particle[i].fy / particle[i].rho;

        //update the position
        /*Euler's method*/
        //particle[i].x += Dt * particle[i].u;
        //particle[i].y += Dt * particle[i].v;
        /*Trapezoidal Rule, for more accuracy*/
        particle[i].x += Dt / 2.0 * (u + particle[i].u);
        particle[i].y += Dt / 2.0 * (v + particle[i].v);

        //handle the boundary conditions
        if (particle[i].x - Boundary_thickness < 0.0) {
            // reset the velocity in opposite direction with half magnitude
            particle[i].u *= Bound_damping;
            particle[i].x = Boundary_thickness;
        }
        if (particle[i].x + Boundary_thickness > View_width) {
            particle[i].u *= Bound_damping;
            particle[i].x = View_width - Boundary_thickness;
        }
        if (particle[i].y - Boundary_thickness < 0.0) {
            particle[i].v *= Bound_damping;
            particle[i].y = Boundary_thickness;
        }
        if (particle[i].y + Boundary_thickness > View_height) {
            particle[i].v *= Bound_damping;
            particle[i].y = View_height - Boundary_thickness;
        }
    }
}

void update()
{
    compute_density_and_pressure();
    compute_forces();
    integrate();

    // redisplay until next loop
    glutPostRedisplay();
}

void init_GL()
{
    glClearColor(0.9, 0.9, 0.9, 1);     // background color
    glEnable(GL_POINT_SMOOTH);          // anti-aliasing 
    glPointSize(H / 2.0);               // set the size of points
    glMatrixMode(GL_PROJECTION);        //use a view projection matrix
}

//called by the windowing system each time we redraw the current frame
void render()
{
    // gives us a blank image buffer state, upon which we can draw the contents of our new frame
    glClear(GL_COLOR_BUFFER_BIT);

    // set the coordinate system for drawing to be an orthographic projection in the screen space
    glLoadIdentity();
    glOrtho(0, View_width, 0, View_height, 0, 1);

    //sets the color in RGBA for all subsequent draw calls (blue)
    glColor4f(0.2, 0.6, 1.0, 1);

    //loop through global array of all fluid particles, placing a vertex at the given x and y coordinate.
    glBegin(GL_POINTS);
    for (int i = 0; i < num; i++) {
        glVertex2f(particle[i].x, particle[i].y);
    }
    glEnd();

    //Performs a buffer swap on the layer in use for the current window
    glutSwapBuffers();
}

void add_particles(int x, int y)
{
    //x = Window_width - x;
    y = Window_height - y;    
    srand(time(0)); // random number
    double jitter = (rand() % 100 - 50) / 50.0; //random number between (-1, 1)
    double y0 = View_height / 1.5 + (rand() % 100 - 50);
    double x0 = View_width / 2.0 + (rand() % 100 - 50);
    x0 = x * 1.0 / Window_width * View_width;
    y0 = y * 1.0 / Window_height * View_height;
    double xi, yi;
    for (int i = 0; i < Small_cube_particles; i++) {
        xi = x0 + (i % 10) * H;
        yi = y0 + (i / 10) * H;

        double jitter = (rand() % 100 - 50) / 100.0; //random number between (-0.5, 0.5)
        if (num < Max_particles)
        {
            particle[num].x = xi + jitter;
            particle[num].y = yi;
            particle[num].u = 0.0;
            particle[num].v = 0.0;
            particle[num].fx = 0.0;
            particle[num].fy = 0.0;
            particle[num].rho = 0;
            particle[num].pressure = 0;
            num++;
        }
    }
}

void keyboard(unsigned char c, int x, int y)
{
    switch (c) {
    case 'a':
        if (num >= Max_particles - Small_cube_particles)
            printf("maximum number of particles reached\n");
        else {
            add_particles(x, y);    //add 100 particles at mouse position
        }
        break;
    case 'r':
    case 'R':
        num = 0;    //clear particles
        init_SPH(); //reset simulation
        break;
    }
}



int main(int argc, char** argv)
{
    initialize_opt();   //initialize the optimization parameters

    //setting the window for render
    glutInitWindowSize(Window_width, Window_height);
    glutInit(&argc, argv);
    glutCreateWindow("My SPH");
    glutDisplayFunc(render);
    glutIdleFunc(update);
    glutKeyboardFunc(keyboard);

    init_GL();
    init_SPH();

    glutMainLoop();
    return 0;
}
