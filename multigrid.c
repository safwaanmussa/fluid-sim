#include <math.h>
#include <stdlib.h>

#define MAX_LEVELS 6

typedef struct {
    int size;
    float *data;
} Grid;

typedef struct {
    int numLevels;
    Grid *grids;
} MultigridSolver;

MultigridSolver* createMultigridSolver(int size) {
    MultigridSolver *solver = malloc(sizeof(MultigridSolver));
    solver->numLevels = 0;
    while (size >= 3 && solver->numLevels < MAX_LEVELS) {
        solver->numLevels++;
        size = (size - 1) / 2 + 1;
    }
    
    solver->grids = malloc(solver->numLevels * sizeof(Grid));
    
    size = (1 << (solver->numLevels - 1)) + 1;
    for (int i = 0; i < solver->numLevels; i++) {
        solver->grids[i].size = size;
        solver->grids[i].data = calloc(size * size, sizeof(float));
        size = (size - 1) / 2 + 1;
    }
    
    return solver;
}

void freeMultigridSolver(MultigridSolver *solver) {
    for (int i = 0; i < solver->numLevels; i++) {
        free(solver->grids[i].data);
    }
    free(solver->grids);
    free(solver);
}

void restrictGrid(Grid *fine, Grid *coarse) {
    int n = coarse->size;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            coarse->data[i + j * n] = 
                0.25f * fine->data[2*i + 2*j * fine->size] +
                0.125f * (fine->data[2*i-1 + 2*j * fine->size] +
                          fine->data[2*i+1 + 2*j * fine->size] +
                          fine->data[2*i + (2*j-1) * fine->size] +
                          fine->data[2*i + (2*j+1) * fine->size]) +
                0.0625f * (fine->data[2*i-1 + (2*j-1) * fine->size] +
                           fine->data[2*i+1 + (2*j-1) * fine->size] +
                           fine->data[2*i-1 + (2*j+1) * fine->size] +
                           fine->data[2*i+1 + (2*j+1) * fine->size]);
        }
    }
}

void prolongate(Grid *coarse, Grid *fine) {
    int n = coarse->size;
    for (int j = 0; j < n - 1; j++) {
        for (int i = 0; i < n - 1; i++) {
            fine->data[2*i + 2*j * fine->size] += coarse->data[i + j * n];
            fine->data[2*i+1 + 2*j * fine->size] += 0.5f * (coarse->data[i + j * n] + coarse->data[i+1 + j * n]);
            fine->data[2*i + (2*j+1) * fine->size] += 0.5f * (coarse->data[i + j * n] + coarse->data[i + (j+1) * n]);
            fine->data[2*i+1 + (2*j+1) * fine->size] += 0.25f * (coarse->data[i + j * n] + coarse->data[i+1 + j * n] +
                                                                 coarse->data[i + (j+1) * n] + coarse->data[i+1 + (j+1) * n]);
        }
    }
}

void relax(Grid *grid, Grid *rhs, float alpha, int iterations) {
    int n = grid->size;
    float recip = 1.0f / (1.0f + 4.0f * alpha);
    for (int iter = 0; iter < iterations; iter++) {
        for (int j = 1; j < n - 1; j++) {
            for (int i = 1; i < n - 1; i++) {
                grid->data[i + j * n] = (rhs->data[i + j * n] + 
                    alpha * (grid->data[i-1 + j * n] + grid->data[i+1 + j * n] +
                             grid->data[i + (j-1) * n] + grid->data[i + (j+1) * n])) * recip;
            }
        }
    }
}

void vCycle(MultigridSolver *solver, Grid *rhs, int level) {
    Grid *current = &solver->grids[level];
    
    if (level == solver->numLevels - 1) {
        relax(current, rhs, 1.0f, 50);  // More iterations at coarsest level
    } else {
        relax(current, rhs, 1.0f, 3);
        
        Grid *coarse = &solver->grids[level + 1];
        Grid coarseRHS = {coarse->size, malloc(coarse->size * coarse->size * sizeof(float))};
        
        // Compute residual and restrict to coarser grid
        for (int j = 1; j < current->size - 1; j++) {
            for (int i = 1; i < current->size - 1; i++) {
                float res = rhs->data[i + j * current->size] - 
                    (current->data[i + j * current->size] -
                     0.25f * (current->data[i-1 + j * current->size] + current->data[i+1 + j * current->size] +
                              current->data[i + (j-1) * current->size] + current->data[i + (j+1) * current->size]));
                coarseRHS.data[(i/2) + (j/2) * coarse->size] = res;
            }
        }
        
        for (int i = 0; i < coarse->size * coarse->size; i++) {
            coarse->data[i] = 0;  // Initial guess
        }
        
        vCycle(solver, &coarseRHS, level + 1);
        
        prolongate(coarse, current);
        free(coarseRHS.data);
        
        relax(current, rhs, 1.0f, 3);
    }
}

void solvePoissonMultigrid(MultigridSolver *solver, float *x, float *b) {
    int n = solver->grids[0].size;
    Grid *finest = &solver->grids[0];
    Grid rhs = {n, b};
    
    for (int i = 0; i < n * n; i++) {
        finest->data[i] = x[i];  // Initial guess
    }
    
    for (int i = 0; i < 10; i++) {  // V-cycle iterations
        vCycle(solver, &rhs, 0);
    }
    
    for (int i = 0; i < n * n; i++) {
        x[i] = finest->data[i];
    }
}