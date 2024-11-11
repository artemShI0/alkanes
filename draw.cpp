#include "raylib.h"

int main()
{
    InitWindow(800, 600, "Hello, World!");
    
    while (!WindowShouldClose())
    {
        BeginDrawing();
        
        ClearBackground(RAYWHITE);
        DrawText("Hello, Raylib!", 400, 300, 20, LIGHTGRAY);
        
        EndDrawing();
    }

    CloseWindow();

    return 0;
}

