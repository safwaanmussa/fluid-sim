// multigrid.h
#ifndef MULTIGRID_H
#define MULTIGRID_H

// Grid structure for multigrid solver
typedef struct {
    int size;
    float *data;
} Grid;

// Multigrid solver structure
typedef struct {
    int numLevels;
    Grid *grids;
} MultigridSolver;

// Fluid grid structure
typedef struct {
    float *velocityX, *velocityY, *prevVelocityX, *prevVelocityY;
    float *density, *prevDensity;
    float *pressure, *divergence;
    MultigridSolver *solver;
} FluidGrid;

// Multigrid solver functions
MultigridSolver* createMultigridSolver(int size);
void freeMultigridSolver(MultigridSolver *solver);
void solvePoissonMultigrid(MultigridSolver *solver, float *x, float *b);

// Fluid grid functions
FluidGrid* createFluidGrid(void);
void freeFluidGrid(FluidGrid *grid);
void setBoundaryConditions(int boundaryType, float *x);
void diffuse(int boundaryType, float *x, float *x0, float diffRate, float dt);
void advect(int boundaryType, float *d, float *d0, float *velocityX, float *velocityY, float dt);
void projectMultigrid(MultigridSolver *solver, float *velocityX, float *velocityY, float *pressure, float *divergence);
void updateFluidGrid(FluidGrid *grid);
void addDensity(FluidGrid *grid, int x, int y, float amount);
void addVelocity(FluidGrid *grid, int x, int y, float amountX, float amountY);
void renderFluidGrid(FluidGrid *grid);

#endif // MULTIGRID_H