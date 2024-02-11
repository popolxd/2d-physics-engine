#ifndef BODY_H
#define BODY_H

#include "raylib.h"
#include "polygon.h"

typedef struct {
    POLYGON_Polygon polygon;
    Vector2 vel;
    float angularVel;
    float mass;
    float momentOfInertia;
} BODY_Body;

typedef struct {
    BODY_Body body;
    float acc;
    float angularAcc;
} BODY_Player;

Vector2 BODY_GetNormal(Vector2 vec);
float BODY_DotProdVec2(Vector2 vec1, Vector2 vec2);
float BODY_CrossProdVec2(Vector2 vec1, Vector2 vec2);

void BODY_CreateBody(BODY_Body* body, Vector2* points, int pointsLen, Vector2 initVel, float initAngularVel);
void BODY_CreatePlayer(BODY_Player* player, Vector2* points, int pointsLen, float acc, float angularAcc);

void BODY_UpdateBody(BODY_Body* body);
void BODY_UpdatePlayerBody(BODY_Player* player, BODY_Body* bodies);

void BODY_ResolveCollisionWithWall(BODY_Body* body);

float BODY_CheckCollision(BODY_Body* body1, BODY_Body* body2, Vector2* body1pointsBefore, Vector2* body2PointsBefore, Vector2* pointOfIntersection, Vector2* normal);
void BODY_MoveAllAndResolveAllCollisions(BODY_Player* player, BODY_Body* bodies, int bodiesLen);
void BODY_ResolveCollision(BODY_Body* body1, BODY_Body* body2, Vector2 pointOfIntersection, Vector2 normal, float time, float elasticity);
#endif