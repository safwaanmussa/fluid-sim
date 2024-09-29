#include <raylib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "multigrid.c"

#define GRID_INDEX(x, y) ((x) + (y) * GRID_SIZE)
#define SWAP_POINTERS(x0, x) {float *tmp = x0; x0 = x; x = tmp;}

const int GRID_SIZE = 128;
const int CELL_SIZE = 4;
const float TIME_STEP = 0.1f;
const float DIFFUSION_RATE = 0.0f;
const float VISCOSITY = 0.0001f;

typedef struct {
    float *velocityX, *velocityY, *prevVelocityX, *prevVelocityY;
    float *density, *prevDensity;
    float *pressure, *divergence;
    MultigridSolver *solver;
} FluidGrid;

FluidGrid *createFluidGrid(void) {
    FluidGrid *grid = malloc(sizeof(*grid));
    
    int totalCells = GRID_SIZE * GRID_SIZE;
    
    grid->velocityX = calloc(totalCells, sizeof(float));
    grid->velocityY = calloc(totalCells, sizeof(float));
    grid->prevVelocityX = calloc(totalCells, sizeof(float));
    grid->prevVelocityY = calloc(totalCells, sizeof(float));
    grid->density = calloc(totalCells, sizeof(float));
    grid->prevDensity = calloc(totalCells, sizeof(float));
    grid->pressure = calloc(totalCells, sizeof(float));
    grid->divergence = calloc(totalCells, sizeof(float));
    
    grid->solver = createMultigridSolver(GRID_SIZE);
    
    return grid;
}

void freeFluidGrid(FluidGrid *grid) {
    free(grid->velocityX);
    free(grid->velocityY);
    free(grid->prevVelocityX);
    free(grid->prevVelocityY);
    free(grid->density);
    free(grid->prevDensity);
    free(grid->pressure);
    free(grid->divergence);
    freeMultigridSolver(grid->solver);
    free(grid);
}

void setBoundaryConditions(int boundaryType, float *x) {
    for (int i = 1; i < GRID_SIZE - 1; i++) {
        x[GRID_INDEX(i, 0)] = boundaryType == 2 ? -x[GRID_INDEX(i, 1)] : x[GRID_INDEX(i, 1)];
        x[GRID_INDEX(i, GRID_SIZE-1)] = boundaryType == 2 ? -x[GRID_INDEX(i, GRID_SIZE-2)] : x[GRID_INDEX(i, GRID_SIZE-2)];
    }
    for (int j = 1; j < GRID_SIZE - 1; j++) {
        x[GRID_INDEX(0, j)] = boundaryType == 1 ? -x[GRID_INDEX(1, j)] : x[GRID_INDEX(1, j)];
        x[GRID_INDEX(GRID_SIZE-1, j)] = boundaryType == 1 ? -x[GRID_INDEX(GRID_SIZE-2, j)] : x[GRID_INDEX(GRID_SIZE-2, j)];
    }
    
    x[GRID_INDEX(0, 0)] = 0.5f * (x[GRID_INDEX(1, 0)] + x[GRID_INDEX(0, 1)]);
    x[GRID_INDEX(0, GRID_SIZE-1)] = 0.5f * (x[GRID_INDEX(1, GRID_SIZE-1)] + x[GRID_INDEX(0, GRID_SIZE-2)]);
    x[GRID_INDEX(GRID_SIZE-1, 0)] = 0.5f * (x[GRID_INDEX(GRID_SIZE-2, 0)] + x[GRID_INDEX(GRID_SIZE-1, 1)]);
    x[GRID_INDEX(GRID_SIZE-1, GRID_SIZE-1)] = 0.5f * (x[GRID_INDEX(GRID_SIZE-2, GRID_SIZE-1)] + x[GRID_INDEX(GRID_SIZE-1, GRID_SIZE-2)]);
}

void diffuse(int boundaryType, float *x, float *x0, float diffRate, float dt) {
    float a = dt * diffRate * (GRID_SIZE - 2) * (GRID_SIZE - 2);
    // Instead of solveLinearSystem, use the multigrid solver
    for (int i = 0; i < 20; i++) {  // Perform multiple iterations
        setBoundaryConditions(boundaryType, x);
        for (int j = 1; j < GRID_SIZE - 1; j++) {
            for (int i = 1; i < GRID_SIZE - 1; i++) {
                x[GRID_INDEX(i, j)] = (x0[GRID_INDEX(i, j)] + a * (
                    x[GRID_INDEX(i-1, j)] + x[GRID_INDEX(i+1, j)] +
                    x[GRID_INDEX(i, j-1)] + x[GRID_INDEX(i, j+1)]
                )) / (1 + 4 * a);
            }
        }
    }
}

void advect(int boundaryType, float *d, float *d0, float *velocityX, float *velocityY, float dt) {
    float dt0 = dt * GRID_SIZE;
    for (int j = 1; j < GRID_SIZE - 1; j++) {
        for (int i = 1; i < GRID_SIZE - 1; i++) {
            float x = i - dt0 * velocityX[GRID_INDEX(i, j)];
            float y = j - dt0 * velocityY[GRID_INDEX(i, j)];
            
            if (x < 0.5f) x = 0.5f; 
            if (x > GRID_SIZE + 0.5f) x = GRID_SIZE + 0.5f;
            int i0 = (int)x; 
            int i1 = i0 + 1;
            
            if (y < 0.5f) y = 0.5f; 
            if (y > GRID_SIZE + 0.5f) y = GRID_SIZE + 0.5f;
            int j0 = (int)y; 
            int j1 = j0 + 1;
            
            float s1 = x - i0; 
            float s0 = 1 - s1;
            float t1 = y - j0; 
            float t0 = 1 - t1;
            
            d[GRID_INDEX(i, j)] = s0 * (t0 * d0[GRID_INDEX(i0, j0)] + t1 * d0[GRID_INDEX(i0, j1)]) +
                                  s1 * (t0 * d0[GRID_INDEX(i1, j0)] + t1 * d0[GRID_INDEX(i1, j1)]);
        }
    }
    setBoundaryConditions(boundaryType, d);
}

void projectMultigrid(MultigridSolver *solver, float *velocityX, float *velocityY, float *pressure, float *divergence) {
    int n = GRID_SIZE;
    for (int j = 1; j < n - 1; j++) {
        for (int i = 1; i < n - 1; i++) {
            divergence[GRID_INDEX(i, j)] = -0.5f * (
                velocityX[GRID_INDEX(i+1, j)] - velocityX[GRID_INDEX(i-1, j)] +
                velocityY[GRID_INDEX(i, j+1)] - velocityY[GRID_INDEX(i, j-1)]
            ) / n;
            pressure[GRID_INDEX(i, j)] = 0;
        }
    }
    setBoundaryConditions(0, divergence);
    setBoundaryConditions(0, pressure);
    
    solvePoissonMultigrid(solver, pressure, divergence);
    
    for (int j = 1; j < n - 1; j++) {
        for (int i = 1; i < n - 1; i++) {
            velocityX[GRID_INDEX(i, j)] -= 0.5f * n * (pressure[GRID_INDEX(i+1, j)] - pressure[GRID_INDEX(i-1, j)]);
            velocityY[GRID_INDEX(i, j)] -= 0.5f * n * (pressure[GRID_INDEX(i, j+1)] - pressure[GRID_INDEX(i, j-1)]);
        }
    }
    setBoundaryConditions(1, velocityX);
    setBoundaryConditions(2, velocityY);
}

void updateFluidGrid(FluidGrid *grid) {
    // Velocity step
    SWAP_POINTERS(grid->velocityX, grid->prevVelocityX);
    SWAP_POINTERS(grid->velocityY, grid->prevVelocityY);
    
    diffuse(1, grid->velocityX, grid->prevVelocityX, VISCOSITY, TIME_STEP);
    diffuse(2, grid->velocityY, grid->prevVelocityY, VISCOSITY, TIME_STEP);
    
    projectMultigrid(grid->solver, grid->velocityX, grid->velocityY, grid->pressure, grid->divergence);
    
    SWAP_POINTERS(grid->velocityX, grid->prevVelocityX);
    SWAP_POINTERS(grid->velocityY, grid->prevVelocityY);
    
    advect(1, grid->velocityX, grid->prevVelocityX, grid->prevVelocityX, grid->prevVelocityY, TIME_STEP);
    advect(2, grid->velocityY, grid->prevVelocityY, grid->prevVelocityX, grid->prevVelocityY, TIME_STEP);
    
    projectMultigrid(grid->solver, grid->velocityX, grid->velocityY, grid->pressure, grid->divergence);
    
    // Density step
    SWAP_POINTERS(grid->density, grid->prevDensity);
    diffuse(0, grid->density, grid->prevDensity, DIFFUSION_RATE, TIME_STEP);
    SWAP_POINTERS(grid->density, grid->prevDensity);
    advect(0, grid->density, grid->prevDensity, grid->velocityX, grid->velocityY, TIME_STEP);
}

void addDensity(FluidGrid *grid, int x, int y, float amount) {
    grid->density[GRID_INDEX(x, y)] += amount;
}

void addVelocity(FluidGrid *grid, int x, int y, float amountX, float amountY) {
    int index = GRID_INDEX(x, y);
    grid->velocityX[index] += amountX;
    grid->velocityY[index] += amountY;
}

int main(void) {
    const int screenWidth = GRID_SIZE * CELL_SIZE;
    const int screenHeight = GRID_SIZE * CELL_SIZE;
    
    InitWindow(screenWidth, screenHeight, "Fluid Simulation");
    SetTargetFPS(60);
    
    FluidGrid *fluidGrid = createFluidGrid();
    
    Vector2 prevMousePos = {0, 0};
    
    while (!WindowShouldClose()) {
        Vector2 mousePos = GetMousePosition();
        
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            int cellX = (int)(mousePos.x / CELL_SIZE);
            int cellY = (int)(mousePos.y / CELL_SIZE);
            addDensity(fluidGrid, cellX, cellY, 100);
            
            float velX = (mousePos.x - prevMousePos.x) * 2;
            float velY = (mousePos.y - prevMousePos.y) * 2;
            addVelocity(fluidGrid, cellX, cellY, velX, velY);
        }
        
        updateFluidGrid(fluidGrid);
        
        BeginDrawing();
        ClearBackground(BLACK);
        
        for (int y = 0; y < GRID_SIZE; y++) {
            for (int x = 0; x < GRID_SIZE; x++) {
                float d = fluidGrid->density[GRID_INDEX(x, y)];
                Color color = {(unsigned char)(d * 255), 0, 0, 255};
                DrawRectangle(x * CELL_SIZE, y * CELL_SIZE, CELL_SIZE, CELL_SIZE, color);
            }
        }
        
        DrawFPS(10, 10);
        EndDrawing();
        
        prevMousePos = mousePos;
    }
    
    freeFluidGrid(fluidGrid);
    CloseWindow();
    
    return 0;
}
