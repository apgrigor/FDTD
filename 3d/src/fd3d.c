#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>

#define Nx 100
#define Ny 100
#define Nz 100
#define steps 700
#define total_frames 350
#define mur1 (1 / 3.)
#define mur2 (4 / 3.)

double dz[Nx][Ny][Nz];
double ez[Nx][Ny][Nz];
double iz[Nx][Ny][Nz];
double sz[Nx][Ny][Nz];

double hx[Nx][Ny][Nz];
double hy[Nx][Ny][Nz];
double bx[Nx][Ny][Nz];
double by[Nx][Ny][Nz];
double ix[Nx][Ny][Nz];
double iy[Nx][Ny][Nz];
double sx[Nx][Ny][Nz];
double sy[Nx][Ny][Nz];

double prev_ez_x_lower[2][2][Ny][Nz];
double prev_ez_x_upper[2][2][Ny][Nz];
double prev_ez_y_lower[2][2][Nx][Nz];
double prev_ez_y_upper[2][2][Nx][Nz];
int prev = 0;
int before_prev = 1;
int old_pprev;

const double dx = 2.5e-5;
const double dt = dx / (6e8 + 0.);
const double wp = 2 * 3.14 * 1e12;
const double vc = wp / 1000;
const double w0 = wp / 1.4142135623730951;
const double d1 = wp * wp * dt / vc;
const double d2 = (1 - 0.5 * dt * vc) / (1 + 0.5 * dt * vc);
const unsigned int ylower = 40;
const unsigned int yupper = 70;


void mur_x(
    double prev_ez_x[2][2][Ny][Nz], double ez[Nx][Ny][Nz], 
    int bound, int pbound, int prev, int before_prev
) {
    for(int j = 0; j < Ny; j++) {
        for(int k = 0; k < Nz; k++) {
            ez[bound][j][k] = -prev_ez_x[before_prev][1][j][k] - mur1 * (
                ez[pbound][j][k] + prev_ez_x[before_prev][0][j][k]
            ) + mur2 * (
                prev_ez_x[prev][0][j][k] + prev_ez_x[prev][1][j][k]
            );

            prev_ez_x[before_prev][0][j][k] = ez[bound][j][k];
            prev_ez_x[before_prev][1][j][k] = ez[pbound][j][k];
        }
    }
}

void mur_y(
    double prev_ez_y[2][2][Nx][Nz], double ez[Nx][Ny][Nz], 
    int bound, int pbound, int prev, int before_prev
) {
    for(int i = 0; i < Nx; i++) {
        for(int k = 0; k < Nz; k++) {
            ez[i][bound][k] = -prev_ez_y[before_prev][1][i][k] - mur1 * (
                ez[i][pbound][k] + prev_ez_y[before_prev][0][i][k]
            ) + mur2 * (
                prev_ez_y[prev][0][i][k] + prev_ez_y[prev][1][i][k]
            );

            prev_ez_y[before_prev][0][i][k] = ez[i][bound][k];
            prev_ez_y[before_prev][1][i][k] = ez[i][pbound][k];
        }
    }
}


static PyObject* fd3d_run_sim(PyObject *self, PyObject *args) {
    const unsigned int steps_per_frame = steps / total_frames;
    int current_frame = 0;

    npy_intp z50_dims[] = {total_frames, Nx, Ny};
    PyObject *z50_frames = PyArray_SimpleNew(3, z50_dims, NPY_DOUBLE);
    double *z50_data = PyArray_DATA((PyArrayObject *) z50_frames);

    // initialize to zero
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                dz[i][j][k] = 0;
                ez[i][j][k] = 0;
                iz[i][j][k] = 0;
                sz[i][j][k] = 0;

                hx[i][j][k] = 0;
                hy[i][j][k] = 0;
                bx[i][j][k] = 0;
                by[i][j][k] = 0;
                ix[i][j][k] = 0;
                iy[i][j][k] = 0;
                sx[i][j][k] = 0;
                sy[i][j][k] = 0;
            }
        }
    }

    // main loop
    for(int step = 1; step < steps + 1; step++) {
        for(int i = 1; i < Nx; i++) {
            for(int j = 1; j < Ny; j++) {
                for(int k = 1; k < Nz; k++) {
                    dz[i][j][k] += 0.5 * (
                        hy[i][j][k] - hy[i - 1][j][k] - 
                        hx[i][j][k] + hx[i][j - 1][k]
                    );
                }
            }
        }
        // fprintf(stdout, "hello %d\n", step);
        dz[50][25][50] = sin(w0 * dt * step);

        for(int i = 0; i < Nx; i++) {
            for(int j = 0; j < ylower; j++) {
                for(int k = 0; k < Nz; k++) {
                    ez[i][j][k] = dz[i][j][k];
                }
            }
        }
        for(int i = 0; i < Nx; i++) {
            for(int j = ylower; j < yupper+1; j++) {
                for(int k = 0; k < Nz; k++) {
                    ez[i][j][k] = dz[i][j][k] - iz[i][j][k] - d2 * sz[i][j][k];
                    iz[i][j][k] += d1 * ez[i][j][k];
                    sz[i][j][k] = d2 * sz[i][j][k] - d1 * ez[i][j][k];
                }
            }
        }
        for(int i = 0; i < Nx; i++) {
            for(int j = yupper+1; j < Ny; j++) {
                for(int k = 0; k < Nz; k++) {
                    ez[i][j][k] = dz[i][j][k];
                }
            }
        }

        mur_x(prev_ez_x_lower, ez, 0, 1, prev, before_prev);
        mur_x(prev_ez_x_upper, ez, Nx - 1, Nx - 2, prev, before_prev);

        mur_y(prev_ez_y_lower, ez, 0, 1, prev, before_prev);
        mur_y(prev_ez_y_upper, ez, Ny - 1, Ny - 2, prev, before_prev);

        old_pprev = before_prev;
        before_prev = prev;
        prev = old_pprev;

        for(int i = 0; i < Nx; i++) {
            for(int j = 0; j < Ny-1; j++) {
                for(int k = 0; k < Nz-1; k++) {
                    bx[i][j][k] += 0.5 * (ez[i][j][k] - ez[i][j+1][k]);
                }
            }
        }
        for(int i = 0; i < Nx-1; i++) {
            for(int j = 0; j < Ny; j++) {
                for(int k = 0; k < Nz-1; k++) {
                    by[i][j][k] += 0.5 * (ez[i+1][j][k] - ez[i][j][k]);
                }
            }
        }


        for(int i = 0; i < Nx-1; i++) {
            for(int j = 0; j < ylower; j++) {
                for(int k = 0; k < Nz-1; k++) {
                    hx[i][j][k] = bx[i][j][k];
                    hy[i][j][k] = by[i][j][k];
                }
            }
        }
        for(int i = 0; i < Nx; i++) {
            for(int j = ylower; j < yupper+1; j++) {
                for(int k = 0; k < Nz; k++) {
                    hx[i][j][k] = bx[i][j][k] - ix[i][j][k] - d2 * sx[i][j][k];
                    ix[i][j][k] += d1 * hx[i][j][k];
                    sx[i][j][k] = d2 * sx[i][j][k] - d1 * hx[i][j][k];

                    hy[i][j][k] = by[i][j][k] - iy[i][j][k] - d2 * sy[i][j][k];
                    iy[i][j][k] += d1 * hy[i][j][k];
                    sy[i][j][k] = d2 * sy[i][j][k] - d1 * hy[i][j][k];
                }
            }
        }
        for(int i = 0; i < Nx; i++) {
            for(int j = yupper+1; j < Ny; j++) {
                for(int k = 0; k < Nz; k++) {
                    hx[i][j][k] = bx[i][j][k];
                    hy[i][j][k] = by[i][j][k];
                }
            }
        }
        
        if(step % steps_per_frame == 0) {
            for(int i = 0; i < Nx; i++) {
                for(int j = 0; j < Ny; j++) {
                    z50_data[
                        Ny * Nx * current_frame + Nx * i + j
                    ] = ez[i][j][50];
                }
            }
            if(current_frame % 10 == 0) {
                fprintf(stdout, "creating frame %d\n", current_frame);
            }
            current_frame++;
        }
    }

    return z50_frames;
}



// CPython boilerplate

static PyMethodDef FD3DMethods[] = {
    {"run_sim", fd3d_run_sim, METH_VARARGS, "runs FD3D simulation"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef fd3dmodule = {
    PyModuleDef_HEAD_INIT,
    "fd3d",
    "fd3d simulation module",
    -1,
    FD3DMethods
};

PyMODINIT_FUNC PyInit_fd3d(void) {
    PyObject *m;

    m = PyModule_Create(&fd3dmodule);
    if(m == NULL) {
        printf("PyModule_Create failed\n");
        return NULL;
    }

    import_array();

    return m;
}