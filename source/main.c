#include "main.h"
#include "raylib.h"
#include "body.h"
#include <stdlib.h>

int SCREEN_WIDTH;
int SCREEN_HEIGHT;
float ELAPSED;
float EPSILON = 0.01;

float numOfnumsOnScreen;
Vector2 off;

int main(void)
{
    SCREEN_WIDTH = 1200;
    SCREEN_HEIGHT = 800;
    ELAPSED = 1;

    int bodiesLen = 1;

    BODY_Player player;
    Vector2 points0[3] = {(Vector2){300, 300}, (Vector2){400, 400}, (Vector2){500, 300}};
    BODY_CreatePlayer(&player, points0, 3, 100, 0.5);

    BODY_Body* bodies = (BODY_Body*)malloc(sizeof(BODY_Body) * bodiesLen);
    Vector2 points1[3] = {(Vector2){100, 100}, (Vector2){150, 250}, (Vector2){200, 100}};
    BODY_CreateBody(&bodies[0], points1, 3, (Vector2){50, 50}, 0);

    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "2d physics engine");

    SetTargetFPS(60);

    // Main game loop
    while (!WindowShouldClose())
    {
        ELAPSED = GetFrameTime();

        BeginDrawing();
        ClearBackground(BLACK);

        BODY_MoveAllAndResolveAllCollisions(&player, bodies, bodiesLen);

        EndDrawing();
    }

    CloseWindow();

    return 0;
}