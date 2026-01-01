#include<vector>
#include <cmath>
#include <utility>
#include <random>
#include <chrono>

#include "Triangulation.h"

using namespace std;

mt19937 mt_rand{chrono::high_resolution_clock::now().time_since_epoch().count()};
uniform_int_distribution<int> d_dist2{0,2};
uniform_int_distribution<int> d_distTriangles;
uniform_real_distribution<double> d_dist1(0.0, 1.0);


// Flips the edge N° c of triangle t in the triangulation. The flip is supposed to be possible.
void Triangulation::flip(Triangle *t, int c)
{
    Triangle* u=t->neighbor[c], *v1, *v3;
    int d=u->edge(t), c1, c3;
    v1=t->neighbor[(c+1)%3];
    c1=v1->edge(t);
    v3=u->neighbor[(d+1)%3];
    c3=v3->edge(u);
    v1->neighbor[c1] = u;
    v3->neighbor[c3] = t;
    t->neighbor[c] = v3;
    t->neighbor[(c+1)%3] = u;
    u->neighbor[d] = v1;
    u->neighbor[(d+1)%3] = t;
    t->vertex[(c+2)%3] = u->vertex[d];
    u->vertex[(d+2)%3] = t->vertex[c];

    if (t->vertex[(c+1)%3]->t == u)
    {
        t->vertex[(c+1)%3]->t = t;
        t->vertex[(c+1)%3]->e = c;
    }
    if (u->vertex[(d+1)%3]->t == t)
    {
        u->vertex[(d+1)%3]->t = u;
        u->vertex[(d+1)%3]->e = d;
    }
}

// Swaps the colors of the triangles t1 and t2 which are incident in the vertex s.
void Triangulation::triangleColorSwap(Triangle *t1, Triangle *t2, Point *s)
{
    int s1=0, s2=0;
    while (t1->vertex[s1]!=s)
        s1++;
    while (t2->vertex[s2]!=s)
        s2++;
    swap(t1->color,t2->color);
    if (t1->color>0)
    {
        t1->vertex[(s1+1)%3]->t=t1;
        t1->vertex[(s1+1)%3]->e=s1;
        Triangle *t3=t2->neighbor[(s2+1)%3];
        t2->vertex[s2]->t=t3;
        t2->vertex[s2]->e=t3->edge(t2);
        t3=t2->neighbor[(s2+2)%3];
        t2->vertex[(s2+1)%3]->t=t3;
        t2->vertex[(s2+1)%3]->e=t3->edge(t2);
    }
    else
    {
        Triangle *t3=t1->neighbor[s1];
        t1->vertex[(s1+2)%3]->t=t3;
        t1->vertex[(s1+2)%3]->e=t3->edge(t1);
        t2->vertex[s2]->t=t2;
        t2->vertex[s2]->e=(s2+2)%3;
        t2->vertex[(s2+2)%3]->t=t2;
        t2->vertex[(s2+2)%3]->e=(s2+1)%3;
    }
}

// Lists the color swaps that are possible at vertex s. A color swap is possible only if there exists either exactly
// one black or exactly one white triangle with vertex s. t1 stores this triangle (if it exists).
// tri stores the triangles that can be swapped with t1 (if some exist).
void Triangulation::possibleColorSwaps(Point *s, Triangle *&t1, vector<Triangle *> &tri)
{
    tri.clear();
    t1=s->t;
    int e=s->e;
    Triangle *t2;
    if (t1->neighbor[(e+2)%3]->color<=0)  // t1 is the unique black triangle with vertex s
    {
        t2=t1->neighbor[(e+2)%3];
        e=t2->edge(t1);

        while (t2->color<0)
        {
            Triangle *tnext=t2->neighbor[(e+2)%3];
            e=tnext->edge(t2);
            t2=tnext;
        }
        if (t2==t1)         // all other triangles are infinite
            return;
        do
        {
            if (t2->neighbor[(e+1)%3]->color==1)
                tri.push_back(t2);
            Triangle *tnext=t2->neighbor[(e+2)%3];
            e=tnext->edge(t2);
            t2=tnext;
        }
        while (t2!=t1);
    }
    else    // t1 is not the unique black triangle with vertex s
    {
        t2=t1->neighbor[e];

        if (t2->color<0)
            return;

        e=t2->edge(t1);
        if (t2->neighbor[(e+1)%3]->color <=0)
            return;

        // here t2 is the unique white triangle and t1=s->t is its successor in counterclockwise direction

        e=s->e;
        do
        {
            if (t1->neighbor[(e+1)%3]->color<=0)
                tri.push_back(t1);
            Triangle *tnext=t1->neighbor[(e+2)%3];
            e=tnext->edge(t1);
            t1=tnext;
        }
        while (t1!=t2);
    }
}

// Performs a swap color sequence in the given direction: 1 left to right, -1 right to left.
// Returns true if the best polygonization TBest has been improved.
// cmpOp is the operator < for minimization and > for maximization.
template <typename C>
bool Triangulation::colorSwapSequence(double temp, Triangulation &TBest, int direction, C cmpOp)
{
    Triangle *t;
    vector<Triangle*> tri;
    tri.reserve(6);
    bool improved = false;

    int loopBegin, loopEnd;
    if (direction > 0)
    {
        loopBegin=1;
        loopEnd=LP.size();
    }
    else
    {
        loopBegin=LP.size()-1;
        loopEnd=0;
    }
    for (int i=loopBegin; i!=loopEnd; i+=direction)
    {
        possibleColorSwaps(LP[i],t,tri);
        if (!tri.empty())
        {
            double areat2=t->area2();
            do
            {
                int k = mt_rand()%tri.size();
                long long areaVariation2 = t->color==0 ? areat2-tri[k]->area2() : tri[k]->area2()-areat2;
                if (cmpOp(areaVariation2, 0) || ((temp > 0) && (d_dist1(mt_rand) < 1.0/(1.0 + exp(PARAM_SA_C_CS*areaVariation2/temp)))))
                {
                    triangleColorSwap(t,tri[k],LP[i]);
                    area2+=areaVariation2;
                    if (cmpOp(area2,TBest.area2))
                    {
                        TBest.overwriteBy(*this);
                        improved=true;
                    }
                    break;
                }
                swap(tri[k], tri.back());
                tri.pop_back();
            }
            while (!tri.empty());
        }
    }
    return improved;
}

// Performs a sequence of edge flips in the triangulation
void Triangulation::edgeFlipSequence()
{
    int sequenceSize = static_cast<int>(PARAM_SA_C_EF * (LP.size()-1));
    int k=1;
    while (k<=sequenceSize)
    {
        int i=d_distTriangles(mt_rand),
            c=d_dist2(mt_rand);
        Triangle *t=LT[i],
                 *v=t->neighbor[c];
        if (t->color==-1)
        {
            if (v->color!=-1)  // only the edges that are not incident in the point at infinity are counted
                k++;           // the expected number of passes in the loop remains in O(n)
        }
        else
        {
             if (v->color==t->color)
             {
                Point* s=v->vertex[v->edge(t)];
                if (s->onLeft(t->vertex[c],t->vertex[(c+1)%3])>0 && s->onLeft(t->vertex[c],t->vertex[(c+2)%3])<0)
                    flip(t, c);
             }
             k++;
        }
    }
}

// Optimizes the area of the current triangulation.
// The area is minimized if minimization=true and maximized otherwise.
void Triangulation::areaOptimization(bool minimization)
{
    vector<int> Polygon;
    initialPolygonization(Polygon);
    triangulatePolygon(Polygon);
    initArea(Polygon);

    d_distTriangles=uniform_int_distribution<int>{0,static_cast<int>(LT.size()-1)};

    double coeff_c_cs = minimization ? -(PARAM_SA_C_0/areaStdDeviation2())*area2 : (PARAM_SA_C_0/areaStdDeviation2())/area2;
    // coefficient used to compute the signed value of PARAM_SA_C_CS: for minimization it will be used with a positive area variation,
    // and for maximization with a negative one

    auto cmpOp = minimization ? [](double x, double y) { return x<y; } : [](double x, double y) { return x>y; };
    // Comparison function: < for minimization and > for maximization

    int minTempReducSteps = static_cast<int>(round(PARAM_SA_C_SAR / log(PARAM_SA_C_TR)));
        // = (minimum number of iterations between 2 reinitializations of the SA)/PARAM_SA_NB_TR

    Triangulation TBest{*this}; // current best triangulation
    long long bestArea2=area2;  // twice the current best area

    int direction=1; // direction in which the points are traversed in the color-swap sequence: 1 from left to right, -1 from right to left
    int noImprovementBestArea=0;    // number of consecutive steps with no improvement of the best area

    do
    {
        PARAM_SA_C_CS = minimization ? coeff_c_cs/bestArea2 : coeff_c_cs*bestArea2;

        double temp=1.0;                // temperature of the SA
        int nonImprovingIterations = 0; // number of consecutive iterations in loop i with no improvement of the best area
        int tempReducSteps = 0;         // (number of iterations between 2 reinitializations of the SA)/PARAM_SA_NB_TR

        for (int i=1; i<=PARAM_SA_NB_IT; i++)
        {
            for (int j=1;j<=PARAM_SA_NB_CS;j++)
                if (colorSwapSequence(temp, TBest, direction, cmpOp))
                    nonImprovingIterations = 0;
                else
                    nonImprovingIterations++;

            edgeFlipSequence();

            if (i%PARAM_SA_NB_TR==0)
            {
                tempReducSteps++;
                if (nonImprovingIterations>=PARAM_SA_NB_SC && tempReducSteps>=minTempReducSteps)
                {
                    tempReducSteps = 0;     // the SA is reinitialized
                    temp = exp((-i)*PARAM_SA_C_IT / PARAM_SA_NB_IT);
                    if (area2 != TBest.area2)
                        overwriteBy(TBest);  // the current triangulation is set to the best one found so far
                }
                else
                    temp*=PARAM_SA_C_TR;    // the temperature is reduced
            }
        }

        overwriteBy(TBest);
        if (area2 != bestArea2)
        {
            noImprovementBestArea=0;
            bestArea2=area2;
        }
        else
            noImprovementBestArea++;

        if (noImprovementBestArea<2)
        {
            extractPolygon(Polygon);
            deleteTriangles();
            triangulatePolygon(Polygon);
            direction=-direction;
        }
    }
    while (noImprovementBestArea<2);
}

