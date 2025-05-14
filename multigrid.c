// multigrid.c
#include "multigrid.h"
#include <math.h>
#include <stdlib.h>
#include <string.h> // For memset

#define MAX_LEVELS 6
#define IDX(i, j, n) ((i) + (j) * (n))

// Forward declarations of internal functions
static void restrictGrid(Grid *fine, Grid *coarse);
static void prolongate(Grid *coarse, Grid *fine);
static void relax(Grid *grid, Grid *rhs, float alpha, int iterations);
static void vCycle(MultigridSolver *solver, Grid *rhs, int level);

MultigridSolver* createMultigridSolver(int size) {
    MultigridSolver *solver = malloc(sizeof(MultigridSolver));
    if (!solver) return NULL;
    
    // Calculate number of levels based on grid size
    solver->numLevels = 0;
    int currentSize = size;
    while (currentSize >= 3 && solver->numLevels < MAX_LEVELS) {
        solver->numLevels++;
        currentSize = (currentSize - 1) / 2 + 1;
    }
    
    solver->grids = malloc(solver->numLevels * sizeof(Grid));
    if (!solver->grids) {
        free(solver);
        return NULL;
    }
    
    // Initialize grids for each level
    currentSize = size;
    for (int i = 0; i < solver->numLevels; i++) {
        solver->grids[i].size = currentSize;
        solver->grids[i].data = calloc(currentSize * currentSize, sizeof(float));
        
        if (!solver->grids[i].data) {
            // Clean up on allocation failure
            for (int j = 0; j < i; j++) {
                free(solver->grids[j].data);
            }
            free(solver->grids);
            free(solver);
            return NULL;
        }
        
        currentSize = (currentSize - 1) / 2 + 1;
    }
    
    return solver;
}

void freeMultigridSolver(MultigridSolver *solver) {
    if (!solver) return;
    
    if (solver->grids) {
        for (int i = 0; i < solver->numLevels; i++) {
            free(solver->grids[i].data);
        }
        free(solver->grids);
    }
    
    free(solver);
}

// Optimized restrict operation for transferring from fine to coarse grid
static void restrictGrid(Grid *fine, Grid *coarse) {
    int fineSize = fine->size;
    int coarseSize = coarse->size;
    
    // Clear coarse grid before restriction
    memset(coarse->data, 0, coarseSize * coarseSize * sizeof(float));
    
    // Perform weighted restriction (Full-weighting operator)
    for (int j = 1; j < coarseSize - 1; j++) {
        int fineJ = j * 2;
        
        for (int i = 1; i < coarseSize - 1; i++) {
            int fineI = i * 2;
            int coarseIdx = IDX(i, j, coarseSize);
            
            // Center point
            coarse->data[coarseIdx] = 0.25f * fine->data[IDX(fineI, fineJ, fineSize)];
            
            // Direct neighbors (orthogonal)
            if (fineI > 0)
                coarse->data[coarseIdx] += 0.125f * fine->data[IDX(fineI-1, fineJ, fineSize)];
            if (fineI < fineSize-1)
                coarse->data[coarseIdx] += 0.125f * fine->data[IDX(fineI+1, fineJ, fineSize)];
            if (fineJ > 0)
                coarse->data[coarseIdx] += 0.125f * fine->data[IDX(fineI, fineJ-1, fineSize)];
            if (fineJ < fineSize-1)
                coarse->data[coarseIdx] += 0.125f * fine->data[IDX(fineI, fineJ+1, fineSize)];
            
            // Diagonal neighbors
            if (fineI > 0 && fineJ > 0)
                coarse->data[coarseIdx] += 0.0625f * fine->data[IDX(fineI-1, fineJ-1, fineSize)];
            if (fineI < fineSize-1 && fineJ > 0)
                coarse->data[coarseIdx] += 0.0625f * fine->data[IDX(fineI+1, fineJ-1, fineSize)];
            if (fineI > 0 && fineJ < fineSize-1)
                coarse->data[coarseIdx] += 0.0625f * fine->data[IDX(fineI-1, fineJ+1, fineSize)];
            if (fineI < fineSize-1 && fineJ < fineSize-1)
                coarse->data[coarseIdx] += 0.0625f * fine->data[IDX(fineI+1, fineJ+1, fineSize)];
        }
    }
    
    // Handle boundaries
    for (int i = 0; i < coarseSize; i++) {
        coarse->data[IDX(i, 0, coarseSize)] = coarse->data[IDX(i, 1, coarseSize)];
        coarse->data[IDX(i, coarseSize-1, coarseSize)] = coarse->data[IDX(i, coarseSize-2, coarseSize)];
        coarse->data[IDX(0, i, coarseSize)] = coarse->data[IDX(1, i, coarseSize)];
        coarse->data[IDX(coarseSize-1, i, coarseSize)] = coarse->data[IDX(coarseSize-2, i, coarseSize)];
    }
}

// Optimized prolongation (interpolation) from coarse to fine grid
static void prolongate(Grid *coarse, Grid *fine) {
    int coarseSize = coarse->size;
    int fineSize = fine->size;
    
    // For even-indexed cells in fine grid (direct transfer)
    for (int j = 0; j < coarseSize; j++) {
        int fineJ = j * 2;
        if (fineJ >= fineSize) continue;
        
        for (int i = 0; i < coarseSize; i++) {
            int fineI = i * 2;
            if (fineI >= fineSize) continue;
            
            fine->data[IDX(fineI, fineJ, fineSize)] += coarse->data[IDX(i, j, coarseSize)];
        }
    }
    
    // For odd-indexed cells in x-direction
    for (int j = 0; j < coarseSize; j++) {
        int fineJ = j * 2;
        if (fineJ >= fineSize) continue;
        
        for (int i = 0; i < coarseSize - 1; i++) {
            int fineI = i * 2 + 1;
            if (fineI >= fineSize) continue;
            
            fine->data[IDX(fineI, fineJ, fineSize)] +=
                0.5f * (coarse->data[IDX(i, j, coarseSize)] + 
                        coarse->data[IDX(i+1, j, coarseSize)]);
        }
    }
    
    // For odd-indexed cells in y-direction
    for (int j = 0; j < coarseSize - 1; j++) {
        int fineJ = j * 2 + 1;
        if (fineJ >= fineSize) continue;
        
        for (int i = 0; i < coarseSize; i++) {
            int fineI = i * 2;
            if (fineI >= fineSize) continue;
            
            fine->data[IDX(fineI, fineJ, fineSize)] +=
                0.5f * (coarse->data[IDX(i, j, coarseSize)] + 
                        coarse->data[IDX(i, j+1, coarseSize)]);
        }
    }
    
    // For odd-indexed cells in both directions
    for (int j = 0; j < coarseSize - 1; j++) {
        int fineJ = j * 2 + 1;
        if (fineJ >= fineSize) continue;
        
        for (int i = 0; i < coarseSize - 1; i++) {
            int fineI = i * 2 + 1;
            if (fineI >= fineSize) continue;
            
            fine->data[IDX(fineI, fineJ, fineSize)] +=
                0.25f * (coarse->data[IDX(i, j, coarseSize)] + 
                         coarse->data[IDX(i+1, j, coarseSize)] +
                         coarse->data[IDX(i, j+1, coarseSize)] + 
                         coarse->data[IDX(i+1, j+1, coarseSize)]);
        }
    }
}

// Optimized Gauss-Seidel red-black relaxation for solving linear systems
static void relax(Grid *grid, Grid *rhs, float alpha, int iterations) {
    const int n = grid->size;
    const float recip = 1.0f / (1.0f + 4.0f * alpha);
    int i, j, iter;
    
    // Red-black Gauss-Seidel iteration
    for (iter = 0; iter < iterations; iter++) {
        // Red cells (i+j is even)
        for (j = 1; j < n - 1; j++) {
            int startI = (j % 2) == 0 ? 1 : 2;
            for (i = startI; i < n - 1; i += 2) {
                int idx = IDX(i, j, n);
                grid->data[idx] = (rhs->data[idx] +
                    alpha * (grid->data[IDX(i-1, j, n)] + grid->data[IDX(i+1, j, n)] +
                             grid->data[IDX(i, j-1, n)] + grid->data[IDX(i, j+1, n)])) * recip;
            }
        }
        
        // Black cells (i+j is odd)
        for (j = 1; j < n - 1; j++) {
            int startI = (j % 2) == 0 ? 2 : 1;
            for (i = startI; i < n - 1; i += 2) {
                int idx = IDX(i, j, n);
                grid->data[idx] = (rhs->data[idx] +
                    alpha * (grid->data[IDX(i-1, j, n)] + grid->data[IDX(i+1, j, n)] +
                             grid->data[IDX(i, j-1, n)] + grid->data[IDX(i, j+1, n)])) * recip;
            }
        }
        
        // Apply boundary conditions by copying adjacent interior cells
        for (i = 1; i < n - 1; i++) {
            grid->data[IDX(i, 0, n)] = grid->data[IDX(i, 1, n)];
            grid->data[IDX(i, n-1, n)] = grid->data[IDX(i, n-2, n)];
        }
        
        for (j = 0; j < n; j++) {
            grid->data[IDX(0, j, n)] = grid->data[IDX(1, j, n)];
            grid->data[IDX(n-1, j, n)] = grid->data[IDX(n-2, j, n)];
        }
    }
}

// V-cycle multigrid algorithm for solving Poisson equation
static void vCycle(MultigridSolver *solver, Grid *rhs, int level) {
    Grid *current = &solver->grids[level];
    
    // Direct solve at coarsest level
    if (level == solver->numLevels - 1) {
        relax(current, rhs, 1.0f, 30);  // Fewer iterations needed at coarse level
        return;
    }
    
    // Pre-smoothing
    relax(current, rhs, 1.0f, 3);
    
    // Calculate residual
    Grid *coarse = &solver->grids[level + 1];
    int n = current->size;
    
    // Allocate temporary storage for residual
    float *residual = (float*)malloc(n * n * sizeof(float));
    if (!residual) return;
    
    // Calculate residual: r = f - A*u
    for (int j = 1; j < n - 1; j++) {
        for (int i = 1; i < n - 1; i++) {
            int idx = IDX(i, j, n);
            residual[idx] = rhs->data[idx] - (
                current->data[idx] - 0.25f * (
                    current->data[IDX(i-1, j, n)] + 
                    current->data[IDX(i+1, j, n)] +
                    current->data[IDX(i, j-1, n)] + 
                    current->data[IDX(i, j+1, n)]
                )
            );
        }
    }
    
    // Set up coarse grid right-hand side
    Grid residualGrid = {n, residual};
    
    // Restrict residual to coarser grid
    restrictGrid(&residualGrid, coarse);
    
    // Clear coarse grid solution
    memset(coarse->data, 0, coarse->size * coarse->size * sizeof(float));
    
    // Recursive call to solve coarser problem
    vCycle(solver, coarse, level + 1);
    
    // Prolongate correction from coarse to fine grid
    prolongate(coarse, current);
    
    // Post-smoothing
    relax(current, rhs, 1.0f, 3);
    
    // Free temporary storage
    free(residual);
}

// Main function to solve Poisson equation using multigrid method
void solvePoissonMultigrid(MultigridSolver *solver, float *x, float *b) {
    int n = solver->grids[0].size;
    Grid *finest = &solver->grids[0];
    Grid rhs = {n, b};
    
    // Copy initial guess to finest grid
    memcpy(finest->data, x, n * n * sizeof(float));
    
    // Perform multiple V-cycles for convergence
    for (int i = 0; i < 5; i++) {  // Reduced from 10 cycles to 5
        vCycle(solver, &rhs, 0);
    }
    
    // Copy solution back to output array
    memcpy(x, finest->data, n * n * sizeof(float));
}

void freeFluidGrid(FluidGrid *grid) {
    if (!grid) return;
    
    free(grid->velocityX);
    free(grid->velocityY);
    free(grid->prevVelocityX);
    free(grid->prevVelocityY);
    free(grid->density);
    free(grid->prevDensity);
    free(grid->pressure);
    free(grid->divergence);
    
    if (grid->solver)
        freeMultigridSolver(grid->solver);
        
    free(grid);
}