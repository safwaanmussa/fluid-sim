// main.c

#include <raylib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "multigrid.h"

#define GRID_SIZE 128
#define CELL_SIZE 4
#define TIME_STEP 0.1f
#define DIFFUSION_RATE 0.0f
#define VISCOSITY 0.0001f
#define GRID_INDEX(x, y) ((x) + (y) * GRID_SIZE)
#define SWAP_POINTERS(x0, x) {float *tmp = x0; x0 = x; x = tmp;}

FluidGrid *createFluidGrid(void) {
    // Fix: Use correct memory allocation
    FluidGrid *grid = (FluidGrid*)malloc(sizeof(FluidGrid));
    if (!grid) return NULL;
   
    int totalCells = GRID_SIZE * GRID_SIZE;
   
    grid->velocityX = (float*)calloc(totalCells, sizeof(float));
    grid->velocityY = (float*)calloc(totalCells, sizeof(float));
    grid->prevVelocityX = (float*)calloc(totalCells, sizeof(float));
    grid->prevVelocityY = (float*)calloc(totalCells, sizeof(float));
    grid->density = (float*)calloc(totalCells, sizeof(float));
    grid->prevDensity = (float*)calloc(totalCells, sizeof(float));
    grid->pressure = (float*)calloc(totalCells, sizeof(float));
    grid->divergence = (float*)calloc(totalCells, sizeof(float));
   
    // Check for allocation failures
    if (!grid->velocityX || !grid->velocityY || !grid->prevVelocityX || 
        !grid->prevVelocityY || !grid->density || !grid->prevDensity || 
        !grid->pressure || !grid->divergence) {
        freeFluidGrid(grid);
        return NULL;
    }
   
    grid->solver = createMultigridSolver(GRID_SIZE);
    if (!grid->solver) {
        freeFluidGrid(grid);
        return NULL;
    }
   
    return grid;
}

void setBoundaryConditions(int boundaryType, float *x) {
    // Fix: Pre-calculate indices to avoid repeated index calculations
    int i, j;
    int idx1, idx2;
    
    for (i = 1; i < GRID_SIZE - 1; i++) {
        idx1 = GRID_INDEX(i, 0);
        idx2 = GRID_INDEX(i, 1);
        x[idx1] = boundaryType == 2 ? -x[idx2] : x[idx2];
        
        idx1 = GRID_INDEX(i, GRID_SIZE-1);
        idx2 = GRID_INDEX(i, GRID_SIZE-2);
        x[idx1] = boundaryType == 2 ? -x[idx2] : x[idx2];
    }

    for (j = 1; j < GRID_SIZE - 1; j++) {
        idx1 = GRID_INDEX(0, j);
        idx2 = GRID_INDEX(1, j);
        x[idx1] = boundaryType == 1 ? -x[idx2] : x[idx2];
        
        idx1 = GRID_INDEX(GRID_SIZE-1, j);
        idx2 = GRID_INDEX(GRID_SIZE-2, j);
        x[idx1] = boundaryType == 1 ? -x[idx2] : x[idx2];
    }
   
    // Fix: Corner cases
    int idx_00 = GRID_INDEX(0, 0);
    int idx_01 = GRID_INDEX(0, 1);
    int idx_10 = GRID_INDEX(1, 0);
    int idx_0n = GRID_INDEX(0, GRID_SIZE-1);
    int idx_0n1 = GRID_INDEX(0, GRID_SIZE-2);
    int idx_1n = GRID_INDEX(1, GRID_SIZE-1);
    int idx_n0 = GRID_INDEX(GRID_SIZE-1, 0);
    int idx_n1 = GRID_INDEX(GRID_SIZE-1, 1);
    int idx_n10 = GRID_INDEX(GRID_SIZE-2, 0);
    int idx_nn = GRID_INDEX(GRID_SIZE-1, GRID_SIZE-1);
    int idx_nn1 = GRID_INDEX(GRID_SIZE-1, GRID_SIZE-2);
    int idx_n1n = GRID_INDEX(GRID_SIZE-2, GRID_SIZE-1);
    
    x[idx_00] = 0.5f * (x[idx_10] + x[idx_01]);
    x[idx_0n] = 0.5f * (x[idx_1n] + x[idx_0n1]);
    x[idx_n0] = 0.5f * (x[idx_n10] + x[idx_n1]);
    x[idx_nn] = 0.5f * (x[idx_n1n] + x[idx_nn1]);
}

// Optimized diffuse function using fixed iterations
void diffuse(int boundaryType, float *x, float *x0, float diffRate, float dt) {
    // Fix: Pre-calculate constants outside the loop
    float a = dt * diffRate * (GRID_SIZE - 2) * (GRID_SIZE - 2);
    float reciprocal = 1.0f / (1.0f + 4.0f * a);
    int i, j, k, idx, idx_n, idx_s, idx_e, idx_w;
    
    // Fix: Use a reasonable number of iterations (reduced from 20)
    for (k = 0; k < 10; k++) {
        setBoundaryConditions(boundaryType, x);
        
        for (j = 1; j < GRID_SIZE - 1; j++) {
            // Pre-calculate row indices
            int row = j * GRID_SIZE;
            int row_n = (j-1) * GRID_SIZE;
            int row_s = (j+1) * GRID_SIZE;
            
            for (i = 1; i < GRID_SIZE - 1; i++) {
                idx = i + row;
                idx_n = i + row_n;
                idx_s = i + row_s;
                idx_w = (i-1) + row;
                idx_e = (i+1) + row;
                
                x[idx] = (x0[idx] + a * (
                    x[idx_w] + x[idx_e] +
                    x[idx_n] + x[idx_s]
                )) * reciprocal;
            }
        }
    }
}

void advect(int boundaryType, float *d, float *d0, float *velocityX, float *velocityY, float dt) {
    float dt0 = dt * (float)GRID_SIZE;
    int i, j, i0, i1, j0, j1;
    float x, y, s1, s0, t1, t0;
    
    for (j = 1; j < GRID_SIZE - 1; j++) {
        for (i = 1; i < GRID_SIZE - 1; i++) {
            int idx = GRID_INDEX(i, j);
            
            // Fix: backtracing
            x = i - dt0 * velocityX[idx];
            y = j - dt0 * velocityY[idx];
           
            // Fix: Improved boundary clamping
            if (x < 0.5f) x = 0.5f;
            if (x > GRID_SIZE - 1.5f) x = GRID_SIZE - 1.5f;
            i0 = (int)x;
            i1 = i0 + 1;
           
            if (y < 0.5f) y = 0.5f;
            if (y > GRID_SIZE - 1.5f) y = GRID_SIZE - 1.5f;
            j0 = (int)y;
            j1 = j0 + 1;
           
            s1 = x - i0;
            s0 = 1.0f - s1;
            t1 = y - j0;
            t0 = 1.0f - t1;
            
            // Fix: Pre-calculate indices
            int idx00 = GRID_INDEX(i0, j0);
            int idx01 = GRID_INDEX(i0, j1);
            int idx10 = GRID_INDEX(i1, j0);
            int idx11 = GRID_INDEX(i1, j1);
           
            d[idx] = s0 * (t0 * d0[idx00] + t1 * d0[idx01]) +
                     s1 * (t0 * d0[idx10] + t1 * d0[idx11]);
        }
    }
    
    setBoundaryConditions(boundaryType, d);
}

void projectMultigrid(MultigridSolver *solver, float *velocityX, float *velocityY, float *pressure, float *divergence) {
    int i, j, idx, idx_n, idx_s, idx_e, idx_w;
    float halfRecip = -0.5f / GRID_SIZE;
    
    // Calculate divergence
    for (j = 1; j < GRID_SIZE - 1; j++) {
        // Pre-calculate row indices
        int row = j * GRID_SIZE;
        int row_n = (j-1) * GRID_SIZE;
        int row_s = (j+1) * GRID_SIZE;
        
        for (i = 1; i < GRID_SIZE - 1; i++) {
            idx = i + row;
            idx_n = i + row_n;
            idx_s = i + row_s;
            idx_w = (i-1) + row;
            idx_e = (i+1) + row;
            
            divergence[idx] = halfRecip * (
                velocityX[idx_e] - velocityX[idx_w] +
                velocityY[idx_s] - velocityY[idx_n]
            );
            pressure[idx] = 0.0f;
        }
    }
    
    setBoundaryConditions(0, divergence);
    setBoundaryConditions(0, pressure);
   
    // Fix: Use optimized multigrid solver
    solvePoissonMultigrid(solver, pressure, divergence);
   
    // Apply pressure forces to velocity
    float halfGridSize = 0.5f * GRID_SIZE;
    for (j = 1; j < GRID_SIZE - 1; j++) {
        // Pre-calculate row indices
        int row = j * GRID_SIZE;
        int row_n = (j-1) * GRID_SIZE;
        int row_s = (j+1) * GRID_SIZE;
        
        for (i = 1; i < GRID_SIZE - 1; i++) {
            idx = i + row;
            idx_n = i + row_n;
            idx_s = i + row_s;
            idx_w = (i-1) + row;
            idx_e = (i+1) + row;
            
            velocityX[idx] -= halfGridSize * (pressure[idx_e] - pressure[idx_w]);
            velocityY[idx] -= halfGridSize * (pressure[idx_s] - pressure[idx_n]);
        }
    }
    
    setBoundaryConditions(1, velocityX);
    setBoundaryConditions(2, velocityY);
}

void updateFluidGrid(FluidGrid *grid) {
    // Velocity step
    SWAP_POINTERS(grid->prevVelocityX, grid->velocityX);
    SWAP_POINTERS(grid->prevVelocityY, grid->velocityY);
   
    diffuse(1, grid->velocityX, grid->prevVelocityX, VISCOSITY, TIME_STEP);
    diffuse(2, grid->velocityY, grid->prevVelocityY, VISCOSITY, TIME_STEP);
   
    projectMultigrid(grid->solver, grid->velocityX, grid->velocityY, 
                     grid->pressure, grid->divergence);
   
    SWAP_POINTERS(grid->prevVelocityX, grid->velocityX);
    SWAP_POINTERS(grid->prevVelocityY, grid->velocityY);
   
    advect(1, grid->velocityX, grid->prevVelocityX, 
           grid->prevVelocityX, grid->prevVelocityY, TIME_STEP);
    advect(2, grid->velocityY, grid->prevVelocityY, 
           grid->prevVelocityX, grid->prevVelocityY, TIME_STEP);
   
    projectMultigrid(grid->solver, grid->velocityX, grid->velocityY, 
                     grid->pressure, grid->divergence);
   
    // Density step
    SWAP_POINTERS(grid->prevDensity, grid->density);
    diffuse(0, grid->density, grid->prevDensity, DIFFUSION_RATE, TIME_STEP);
    SWAP_POINTERS(grid->prevDensity, grid->density);
    advect(0, grid->density, grid->prevDensity, 
           grid->velocityX, grid->velocityY, TIME_STEP);
}

void addDensity(FluidGrid *grid, int x, int y, float amount) {
    // Fix: Add boundary checking
    if (x < 0 || x >= GRID_SIZE || y < 0 || y >= GRID_SIZE) return;
    
    grid->density[GRID_INDEX(x, y)] += amount;
}

void addVelocity(FluidGrid *grid, int x, int y, float amountX, float amountY) {
    // Fix: Add boundary checking
    if (x < 0 || x >= GRID_SIZE || y < 0 || y >= GRID_SIZE) return;
    
    int index = GRID_INDEX(x, y);
    grid->velocityX[index] += amountX;
    grid->velocityY[index] += amountY;
}

// Fix: Add rendering optimization with multi-color visualization
void renderFluidGrid(FluidGrid *grid) {
    static Image fluidImage = {0};
    static Texture2D fluidTexture = {0};
    static Color *pixels = NULL;
    
    // Initialize texture if needed
    if (!pixels) {
        fluidImage = GenImageColor(GRID_SIZE, GRID_SIZE, BLACK);
        pixels = fluidImage.data;
        fluidTexture = LoadTextureFromImage(fluidImage);
    }
    
    // Update pixels based on density with color gradient for better visualization
    for (int y = 0; y < GRID_SIZE; y++) {
        for (int x = 0; x < GRID_SIZE; x++) {
            int idx = GRID_INDEX(x, y);
            float d = grid->density[idx];
            
            // Clamp density for visualization
            if (d > 1.0f) d = 1.0f;
            
            // Create color gradient: black -> red -> yellow -> white
            unsigned char r = (unsigned char)(d * 255.0f);
            unsigned char g = d > 0.5f ? (unsigned char)((d - 0.5f) * 2.0f * 255.0f) : 0;
            unsigned char b = d > 0.75f ? (unsigned char)((d - 0.75f) * 4.0f * 255.0f) : 0;
            
            pixels[idx] = (Color){r, g, b, 255};
        }
    }
    
    // Update texture and draw - use one draw call instead of per-cell drawing
    UpdateTexture(fluidTexture, pixels);
    DrawTextureEx(fluidTexture, (Vector2){0, 0}, 0.0f, CELL_SIZE, WHITE);
}

int main(void) {
    const int screenWidth = GRID_SIZE * CELL_SIZE;
    const int screenHeight = GRID_SIZE * CELL_SIZE;
   
    InitWindow(screenWidth, screenHeight, "Fluid Simulation");
    SetTargetFPS(60);
   
    FluidGrid *fluidGrid = createFluidGrid();
    if (!fluidGrid) {
        CloseWindow();
        return 1;
    }
   
    Vector2 prevMousePos = {0, 0};
    bool firstFrame = true;
   
    while (!WindowShouldClose()) {
        Vector2 mousePos = GetMousePosition();
       
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            int cellX = (int)(mousePos.x / CELL_SIZE);
            int cellY = (int)(mousePos.y / CELL_SIZE);
            
            // Add density at mouse position
            addDensity(fluidGrid, cellX, cellY, 100.0f);
           
            // Add velocity based on mouse movement
            if (!firstFrame) {
                float velX = (mousePos.x - prevMousePos.x) * 0.5f;  // Reduced scale factor
                float velY = (mousePos.y - prevMousePos.y) * 0.5f;
                addVelocity(fluidGrid, cellX, cellY, velX, velY);
            }
        }
        
        firstFrame = false;
       
        // Update simulation
        updateFluidGrid(fluidGrid);
       
        // Draw
        BeginDrawing();
        ClearBackground(BLACK);
       
        // Optimized rendering
        renderFluidGrid(fluidGrid);
       
        DrawFPS(10, 10);
        EndDrawing();
       
        prevMousePos = mousePos;
    }
   
    freeFluidGrid(fluidGrid);
    CloseWindow();
   
    return 0;
}