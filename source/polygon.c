#include "polygon.h"
#include "raylib.h"
#include "main.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>

void POLYGON_CreatePolygon(POLYGON_Polygon* polygon, Vector2* points, int pointsLen)
{
    polygon->points = (Vector2*)malloc(sizeof(Vector2) * pointsLen);
    polygon->pointsLen = pointsLen;
    
    // points
    for (int i = 0; i < pointsLen; i++) {
        polygon->points[i] = points[i];
    }

    // center of mass
    polygon->centerOfMass = (Vector2){0, 0};
    for (int i = 0; i < pointsLen; i++) {
        polygon->centerOfMass.x += polygon->points[i].x;
        polygon->centerOfMass.y += polygon->points[i].y;
    }
    polygon->centerOfMass.x /= pointsLen;
    polygon->centerOfMass.y /= pointsLen;
}

void POLYGON_DrawPolygon(POLYGON_Polygon* polygon)
{
    // printf("%d\n", polygon->pointsLen);
    for (int i = 0; i < polygon->pointsLen - 1; i++) {
        DrawLineV(polygon->points[i], polygon->points[i + 1], WHITE);
    }
    DrawLineV(polygon->points[polygon->pointsLen - 1], polygon->points[0], WHITE);
    DrawCircleV(polygon->centerOfMass, 3, WHITE);
}

void POLYGON_MovePolygon(POLYGON_Polygon* polygon, Vector2 vel)
{
    for (int i = 0; i < polygon->pointsLen; i++) {
        polygon->points[i].x += vel.x*ELAPSED;
        polygon->points[i].y += vel.y*ELAPSED;
    }
    polygon->centerOfMass.x += vel.x*ELAPSED;
    polygon->centerOfMass.y += vel.y*ELAPSED;
}

void POLYGON_RotatePolygon(POLYGON_Polygon* polygon, float angularVel)
{
    for (int i = 0; i < polygon->pointsLen; i++) {
        polygon->points[i] = (Vector2){
            (polygon->points[i].x - polygon->centerOfMass.x) * cos(angularVel*ELAPSED) - (polygon->points[i].y - polygon->centerOfMass.y) * sin(angularVel*ELAPSED) + polygon->centerOfMass.x,
            (polygon->points[i].x - polygon->centerOfMass.x) * sin(angularVel*ELAPSED) + (polygon->points[i].y - polygon->centerOfMass.y) * cos(angularVel*ELAPSED) + polygon->centerOfMass.y
        };
    }
}

// int POLYGON_IsT2Between0and1(float t1, Vector2 p1, Vector2 p2, Vector2 l1, Vector2 l2, Vector2 l3, Vector2 l4)
// {
//     float t2;
//     float denominator = (l2.y + t1*(l4.y - l2.y)) - (l1.y + t1*(l3.y - l1.y));

//     if (denominator != 0) { // mozem ceknut rovno s y-om
//         t2 = ((p1.y + t1*(p2.y - p1.y)) - (l1.y + t1*(l3.y - l1.y))) / denominator;

//         if (t2 >= 0 && t2 <= 1) return 1;
//         else return 0;

//     } else { // skusim ceknut s x-om
//         denominator = (l2.x + t1*(l4.x - l2.x)) - (l1.x + t1*(l3.x - l1.x));

//         if (denominator != 0) { // mozem ceknut s x-om
//             t2 = ((p1.x + t1*(p2.x - p1.x)) - (l1.x + t1*(l3.x - l1.x))) / denominator;

//             if (t2 >= 0 && t2 <= 1) return 1;
//             else return 0;
//         }
//     }
    
//     return 0;
// }

// float POLYGON_CheckPointLineIntersection(Vector2 p1, Vector2 p2, Vector2 l1, Vector2 l2, Vector2 l3, Vector2 l4) // toto sa da urcite optimiznut
// {
//     float a = (l3.y + l2.y - l1.y - l4.y)*(p2.x + l2.x - p1.x - l4.x) - (p2.y + l2.y - p1.y - l4.y)*(l3.x + l2.x - l1.x - l4.x);
//     float b = (l1.y - l2.y)*(p2.x + l2.x - p1.x - l4.x) + (p1.x - l2.x)*(l3.y + l2.y - l1.y - l4.y) - (p1.y - l2.y)*(l3.x + l2.x - l1.x - l4.x) - (l1.x - l2.x)*(p2.y + l2.y - p1.y - l4.y);
//     float c = (l1.y - l2.y)*(p1.x - l2.x) - (p1.y - l2.y)*(l1.x - l2.x);

//     // if (fabs(a) > 0.00001) {  // almost linear (b*t1 + c = 0)
//     //     float t1 = - c / b;
//     //     if (t1 >= 0 && t1 <= 1) {
//     //         if (POLYGON_IsT2Between0and1(t1, p1, p2, l1, l2, l3, l4)) return t1;
//     //     }

//     // } else {

//     float discriminant = pow(b, 2) - 4*a*c;
//     if (discriminant < 0) return 2;

//     float t11 = (-b - sqrt(discriminant)) / (2*a);

//     if (t11 >= 0 && t11 <= 1) {
//         if (POLYGON_IsT2Between0and1(t11, p1, p2, l1, l2, l3, l4)) return t11;
//     }

//     float t12 = (-b + sqrt(discriminant)) / (2*a); // vzdy vacsie ako t11
//     if (t12 >= 0 && t12 <= 1) {
//         if (POLYGON_IsT2Between0and1(t12, p1, p2, l1, l2, l3, l4)) return t12;
//     }
//     // }

//     return 2;
// }

double POLYGON_CheckPointLineIntersection(Vector2 pointPosNow, Vector2 linePosNow[2], Vector2 pointPosBefore, Vector2 linePosBefore[2])
{
    double a = (linePosNow[0].y - linePosBefore[0].y + pointPosBefore.y - pointPosNow.y)*(linePosNow[1].x - linePosBefore[1].x + linePosBefore[0].x - linePosNow[0].x)
    + (pointPosNow.x - pointPosBefore.x + linePosBefore[0].x - linePosNow[0].x)*(linePosNow[1].y - linePosBefore[1].y + linePosBefore[0].y - linePosNow[0].y);

    double b = (linePosNow[0].y - linePosBefore[0].y + pointPosBefore.y - pointPosNow.y)*(linePosBefore[1].x - linePosBefore[0].x)
    + (pointPosNow.x - pointPosBefore.x + linePosBefore[0].x - linePosNow[0].x)*(linePosBefore[1].y - linePosBefore[0].y)
    + (pointPosBefore.x - linePosBefore[0].x)*(linePosNow[1].y - linePosBefore[1].y + linePosBefore[0].y - linePosNow[0].y)
    - (pointPosBefore.y - linePosBefore[0].y)*(linePosNow[1].x - linePosBefore[1].x + linePosBefore[0].x - linePosNow[0].x);

    double c = (linePosBefore[1].y - linePosBefore[0].y)*(pointPosBefore.x - linePosBefore[0].x) - (linePosBefore[1].x - linePosBefore[0].x)*(pointPosBefore.y - linePosBefore[0].y);

    if (a == 0) { // je to linearna a ne kvadraticka
        double t = - c / b;

        if (t >= 0 && t <= 1 && POLYGON_IsT2OnLine(pointPosNow, linePosNow, pointPosBefore, linePosBefore, t)) return t; // tu sa to teoreticky moze posrat
        return INFINITY;
    }

    double discriminant = pow(b, 2) - 4*a*c;
    if (discriminant < 0) return INFINITY;

    double t1 = (-b - sqrt(discriminant)) / (2*a);

    if (t1 >= 0 && t1 <= 1 && POLYGON_IsT2OnLine(pointPosNow, linePosNow, pointPosBefore, linePosBefore, t1)) return t1;

    double t2 = (-b + sqrt(discriminant)) / (2*a);

    if (t2 >= 0 && t2 <= 1 && POLYGON_IsT2OnLine(pointPosNow, linePosNow, pointPosBefore, linePosBefore, t2)) return t2;

    return INFINITY;
}

int POLYGON_IsT2OnLine(Vector2 pointPosNow, Vector2 linePosNow[2], Vector2 pointPosBefore, Vector2 linePosBefore[2], double t1)
{
    double t2;

    if ((linePosBefore[1].x - linePosBefore[0].x) + (linePosNow[1].x - linePosBefore[1].x + linePosBefore[0].x - linePosNow[0].x)*t1 != 0) {
        t2 = ((pointPosBefore.x - linePosBefore[0].x) + (pointPosNow.x - pointPosBefore.x + linePosBefore[0].x - linePosNow[0].x)*t1) / ((linePosBefore[1].x - linePosBefore[0].x) + (linePosNow[1].x - linePosBefore[1].x + linePosBefore[0].x - linePosNow[0].x)*t1);
        if (t2 >= 0 && t2 <= 1) return 1;
        return 0;

    } else if ((linePosBefore[1].y - linePosBefore[0].y) + (linePosNow[1].y - linePosBefore[1].y + linePosBefore[0].y - linePosNow[0].y)*t1 != 0) {
        t2 = ((pointPosBefore.y - linePosBefore[0].y) + (pointPosNow.y - pointPosBefore.y + linePosBefore[0].y - linePosNow[0].y)*t1) / ((linePosBefore[1].y - linePosBefore[0].y) + (linePosNow[1].y - linePosBefore[1].y + linePosBefore[0].y - linePosNow[0].y)*t1);
        if (t2 >= 0 && t2 <= 1) return 1;
        return 0;

    } else {
        return 0;
    }
}