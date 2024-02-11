#include "body.h"
#include "main.h"
#include "polygon.h"
#include "raylib.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

Vector2 BODY_GetNormal(Vector2 vec)
{
    float divisor = sqrt(pow(vec.x, 2) + pow(vec.y, 2));
    
    Vector2 result = (Vector2){-vec.y/divisor, vec.x/divisor};

    return result;
}

float BODY_DotProdVec2(Vector2 vec1, Vector2 vec2)
{
    return vec1.x*vec2.x + vec1.y*vec2.y;
}

float BODY_CrossProdVec2(Vector2 vec1, Vector2 vec2) // technicaly not a cross product
{
    return vec1.x*vec2.y - vec1.y*vec2.x;
}

void BODY_CreateBody(BODY_Body* body, Vector2* points, int pointsLen, Vector2 initVel, float initAngularVel)
{
    POLYGON_CreatePolygon(&body->polygon, points, pointsLen);
    body->vel = initVel;
    body->angularVel= initAngularVel;

    // mass
    body->mass = 0;
    float height, side;
    for (int i = 0; i < pointsLen - 1; i++) { // only for convex polygon
        height = BODY_DotProdVec2(
            (Vector2){points[i].x - body->polygon.centerOfMass.x, points[i].y - body->polygon.centerOfMass.y}, BODY_GetNormal((Vector2){points[i + 1].x - points[i].x, points[i + 1].y - points[i].y})
        );
        side = sqrt(pow(points[i + 1].x - points[i].x, 2) + pow(points[i + 1].y - points[i].y , 2));

        body->mass += height*side/2;
    }
    height = BODY_DotProdVec2(
        (Vector2){points[0].x - body->polygon.centerOfMass.x, points[pointsLen - 1].y - body->polygon.centerOfMass.y}, BODY_GetNormal((Vector2){points[0].x - points[pointsLen - 1].x, points[0].y - points[pointsLen - 1].y})
    );
    side = sqrt(pow(points[0].x - points[pointsLen - 1].x, 2) + pow(points[0].y - points[pointsLen - 1].y , 2));
    body->mass += height*side/2;

    // moment of inertia (approximation with thin disk)
    float inertia = 0;
    for (int i = 0; i < pointsLen; i++) {
        inertia += sqrt(pow(body->polygon.centerOfMass.x - body->polygon.points[i].x, 2) + pow(body->polygon.centerOfMass.y - body->polygon.points[i].y, 2));
    }
    inertia = pow(inertia, 2);
    inertia *= (body->mass/(2*pow(pointsLen, 2)));
    body->momentOfInertia = inertia;
}

void BODY_CreatePlayer(BODY_Player* player, Vector2* points, int pointsLen, float acc, float angularAcc)
{
    BODY_CreateBody(&player->body, points, pointsLen, (Vector2){0, 0}, 0);
    player->acc = acc;
    player->angularAcc = angularAcc;
}

void BODY_UpdateBody(BODY_Body* body)
{
    POLYGON_DrawPolygon(&body->polygon);
    POLYGON_MovePolygon(&body->polygon, body->vel);
    POLYGON_RotatePolygon(&body->polygon, body->angularVel);
}

void BODY_UpdatePlayerBody(BODY_Player* player, BODY_Body* bodies)
{
    POLYGON_DrawPolygon(&player->body.polygon);

    if (IsKeyDown(KEY_W)) {
        player->body.vel.y -= player->acc*ELAPSED;
    } else if (IsKeyDown(KEY_S)) {
        player->body.vel.y += player->acc*ELAPSED;
    }

    if (IsKeyDown(KEY_A)) {
        player->body.vel.x -= player->acc*ELAPSED;
    } else if (IsKeyDown(KEY_D)) {
        player->body.vel.x += player->acc*ELAPSED;
    }

    if (IsKeyDown(KEY_Q)) {
        player->body.angularVel -= player->angularAcc*ELAPSED;
    } else if (IsKeyDown(KEY_E)) {
        player->body.angularVel += player->angularAcc*ELAPSED;
    }

    POLYGON_MovePolygon(&player->body.polygon, player->body.vel);
    POLYGON_RotatePolygon(&player->body.polygon, player->body.angularVel);
}

void BODY_ResolveCollisionWithWall(BODY_Body* body) // primitivne zatial
{
    int collision = 0;
    for (int i = 0; i < body->polygon.pointsLen; i++) {
        
        if (body->polygon.points[i].x <= 0 || body->polygon.points[i].x >= SCREEN_WIDTH) {
            body->vel.x = -body->vel.x;
            body->angularVel = -body->angularVel;
            collision = 1;
        }
        if (body->polygon.points[i].y <= 0 || body->polygon.points[i].y >= SCREEN_HEIGHT) {
            body->vel.y = -body->vel.y;
            body->angularVel = -body->angularVel;
            collision = 1;
        }
        if (collision) break;
    }
}

float BODY_CheckCollision(BODY_Body* body1, BODY_Body* body2, Vector2* body1PointsBefore, Vector2* body2PointsBefore, Vector2* pointOfIntersection, Vector2* normal)
{
    for (int j = 0; j < body1->polygon.pointsLen; j++) { // body1 pointy
        for (int k = 0; k < body2->polygon.pointsLen - 1; k++) {

            float time = POLYGON_CheckPointLineIntersection(body1PointsBefore[j], body1->polygon.points[j], body2PointsBefore[k], body2PointsBefore[k + 1], body2->polygon.points[k], body2->polygon.points[k + 1]);
            if (time >= 0 && time <= 1) {
                *pointOfIntersection = (Vector2){body1PointsBefore[j].x + time*(body1->polygon.points[j].x - body1PointsBefore[j].x), body1PointsBefore[j].y + time*(body1->polygon.points[j].y - body1PointsBefore[j].y)};
                *normal = BODY_GetNormal((Vector2){
                    (body2PointsBefore[k + 1].x + time*(body2->polygon.points[k + 1].x - body2PointsBefore[k + 1].x)) - (body2PointsBefore[k].x + time*(body2->polygon.points[k].x - body2PointsBefore[k].x)),
                    (body2PointsBefore[k + 1].y + time*(body2->polygon.points[k + 1].y - body2PointsBefore[k + 1].y)) - (body2PointsBefore[k].y + time*(body2->polygon.points[k].y - body2PointsBefore[k].y))
                });
                return time;
            }
        }

        float time = POLYGON_CheckPointLineIntersection(body1PointsBefore[j], body1->polygon.points[j], body2PointsBefore[body2->polygon.pointsLen - 1], body2PointsBefore[0], body2->polygon.points[body2->polygon.pointsLen - 1], body2->polygon.points[0]);
        if (time >= 0 && time <= 1) {
            *pointOfIntersection = (Vector2){body1PointsBefore[j].x + time*(body1->polygon.points[j].x - body1PointsBefore[j].x), body1PointsBefore[j].y + time*(body1->polygon.points[j].y - body1PointsBefore[j].y)};
            *normal = BODY_GetNormal((Vector2){
                (body2PointsBefore[0].x + time*(body2->polygon.points[0].x - body2PointsBefore[0].x)) - (body2PointsBefore[body2->polygon.pointsLen - 1].x + time*(body2->polygon.points[body2->polygon.pointsLen - 1].x - body2PointsBefore[body2->polygon.pointsLen - 1].x)),
                (body2PointsBefore[0].y + time*(body2->polygon.points[0].y - body2PointsBefore[0].y)) - (body2PointsBefore[body2->polygon.pointsLen - 1].y + time*(body2->polygon.points[body2->polygon.pointsLen - 1].y - body2PointsBefore[body2->polygon.pointsLen - 1].y))
            });
            return time;
        }
    }

    for (int j = 0; j < body2->polygon.pointsLen; j++) { // body2 pointy
        for (int k = 0; k < body1->polygon.pointsLen - 1; k++) {

            float time = POLYGON_CheckPointLineIntersection(body2PointsBefore[j], body2->polygon.points[j], body1PointsBefore[k], body1PointsBefore[k + 1], body1->polygon.points[k], body1->polygon.points[k + 1]);
            if (time >= 0 && time <= 1) {
                *pointOfIntersection = (Vector2){body2PointsBefore[j].x + time*(body2->polygon.points[j].x - body2PointsBefore[j].x), body2PointsBefore[j].y + time*(body2->polygon.points[j].y - body2PointsBefore[j].y)};
                *normal = BODY_GetNormal((Vector2){
                    (body1PointsBefore[k + 1].x + time*(body1->polygon.points[k + 1].x - body1PointsBefore[k + 1].x)) - (body1PointsBefore[k].x + time*(body1->polygon.points[k].x - body1PointsBefore[k].x)),
                    (body1PointsBefore[k + 1].y + time*(body1->polygon.points[k + 1].y - body1PointsBefore[k + 1].y)) - (body1PointsBefore[k].y + time*(body1->polygon.points[k].y - body1PointsBefore[k].y))
                });
                return time;
            }
        }

        float time = POLYGON_CheckPointLineIntersection(body2PointsBefore[j], body2->polygon.points[j], body1PointsBefore[body1->polygon.pointsLen - 1], body1PointsBefore[0], body1->polygon.points[body1->polygon.pointsLen - 1], body1->polygon.points[0]);
        if (time >= 0 && time <= 1) {
            *pointOfIntersection = (Vector2){body2PointsBefore[j].x + time*(body2->polygon.points[j].x - body2PointsBefore[j].x), body2PointsBefore[j].y + time*(body2->polygon.points[j].y - body2PointsBefore[j].y)};
            *normal = BODY_GetNormal((Vector2){
                (body1PointsBefore[0].x + time*(body1->polygon.points[0].x - body1PointsBefore[0].x)) - (body1PointsBefore[body1->polygon.pointsLen - 1].x + time*(body1->polygon.points[body1->polygon.pointsLen - 1].x - body1PointsBefore[body2->polygon.pointsLen - 1].x)),
                (body1PointsBefore[0].y + time*(body1->polygon.points[0].y - body1PointsBefore[0].y)) - (body1PointsBefore[body1->polygon.pointsLen - 1].y + time*(body1->polygon.points[body1->polygon.pointsLen - 1].y - body1PointsBefore[body2->polygon.pointsLen - 1].y))
            });
            return time;
        }
    }
    return 2;
}

void BODY_MoveAllAndResolveAllCollisions(BODY_Player* player, BODY_Body* bodies, int bodiesLen)
{
    Vector2 playerPointsBefore[player->body.polygon.pointsLen];
    Vector2** bodiesPointsBefore = (Vector2**)malloc(sizeof(Vector2*)*bodiesLen); // odpad

    for (int i = 0; i < player->body.polygon.pointsLen; i++) {
        playerPointsBefore[i] = player->body.polygon.points[i];
    }

    BODY_UpdatePlayerBody(player, bodies);

    for (int i = 0; i < bodiesLen; i++) {
        bodiesPointsBefore[i] = (Vector2*)malloc(sizeof(Vector2)*bodies[i].polygon.pointsLen);
        for (int j = 0; j < bodies[i].polygon.pointsLen; j++) {
            bodiesPointsBefore[i][j] = bodies[i].polygon.points[j];
        }

        BODY_UpdateBody(&bodies[i]);
    }

    Vector2 pointOfIntersection;
    Vector2 normal;
    for (int i = 0; i < bodiesLen; i++) { // zatial iba pre playera checkujem
        float time = BODY_CheckCollision(&player->body, &bodies[i], playerPointsBefore, bodiesPointsBefore[i], &pointOfIntersection, &normal);
        if (time >= 0 && time <= 1) {
            // printf("collision! %d\n", GetRandomValue(0, 100));
            BODY_ResolveCollision(&player->body, &bodies[i], pointOfIntersection, normal, time, 1);
        }
        BODY_ResolveCollisionWithWall(&bodies[i]);
    }
    BODY_ResolveCollisionWithWall(&player->body);

    // free bullshit
    for (int i = 0; i < bodiesLen; i++) {
        free(bodiesPointsBefore[i]);
    }
    free(bodiesPointsBefore);
}

void BODY_ResolveCollision(BODY_Body* body1, BODY_Body* body2, Vector2 pointOfIntersection, Vector2 normal, float time, float elasticity)
{
    // rewind a little
    POLYGON_MovePolygon(&body1->polygon, (Vector2){-body1->vel.x * (1 - time + EPSILON), -body1->vel.y * (1 - time + EPSILON)});
    POLYGON_RotatePolygon(&body1->polygon, -body1->angularVel * (1 - time + EPSILON));

    POLYGON_MovePolygon(&body2->polygon, (Vector2){-body2->vel.x * (1 - time + EPSILON), -body2->vel.y * (1 - time + EPSILON)});
    POLYGON_RotatePolygon(&body2->polygon, -body2->angularVel * (1 - time + EPSILON));

    // get important variables
    Vector2 ra = (Vector2){pointOfIntersection.x - body1->polygon.centerOfMass.x, pointOfIntersection.y - body1->polygon.centerOfMass.y}; // vector from center to point of collision
    Vector2 rb = (Vector2){pointOfIntersection.x - body2->polygon.centerOfMass.x, pointOfIntersection.y - body2->polygon.centerOfMass.y};

    float ra_nb_crossprod = BODY_CrossProdVec2(ra, normal);
    float rb_nb_crossprod = BODY_CrossProdVec2(rb, normal);
    float va_nb_dotprod = BODY_DotProdVec2(body1->vel, normal);
    float vb_nb_dotprod = BODY_DotProdVec2(body2->vel, normal);

    // impulse magnitute                                                                         // je sanca ze tu ma byt ra_nb
    float j = (- 1 - elasticity)*(va_nb_dotprod - vb_nb_dotprod + body1->angularVel*ra_nb_crossprod - body2->angularVel*rb_nb_crossprod) / (1/body1->mass + 1/body2->mass + pow(ra_nb_crossprod, 2)/body1->momentOfInertia + pow(rb_nb_crossprod, 2)/body2->momentOfInertia);
    // printf("%f\n", j);

    // apply impulse
    body1->vel.x += j*normal.x / body1->mass;
    body1->vel.y += j*normal.y / body1->mass;
    body1->angularVel += j*ra_nb_crossprod / body1->momentOfInertia;

    body2->vel.x -= j*normal.x / body2->mass;
    body2->vel.y -= j*normal.y / body2->mass;
    body2->angularVel -= j*ra_nb_crossprod / body2->momentOfInertia;
}