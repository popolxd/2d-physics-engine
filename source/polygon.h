#ifndef POLYGON_H
#define POLYGON_H

#include "raylib.h"

typedef struct {
    Vector2* points;
    Vector2 centerOfMass;
    int pointsLen;
    // Vector2* edges;
    // Vector2* normals; toto mozno neskor
} POLYGON_Polygon; // counterclockwise

void POLYGON_CreatePolygon(POLYGON_Polygon* polygon, Vector2* points, int pointsLen);
void POLYGON_DrawPolygon(POLYGON_Polygon* polygon);
void POLYGON_MovePolygon(POLYGON_Polygon* polygon, Vector2 vel);
void POLYGON_RotatePolygon(POLYGON_Polygon* polygon, float angularVel);
float POLYGON_CheckPointLineIntersection(Vector2 p1, Vector2 p2, Vector2 l1, Vector2 l2, Vector2 l3, Vector2 l4);
int POLYGON_IsT2Between0and1(float t1, Vector2 p1, Vector2 p2, Vector2 l1, Vector2 l2, Vector2 l3, Vector2 l4);

#endif