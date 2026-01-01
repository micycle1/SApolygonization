#include <iostream>
#include <fstream>
#include <string>
#include<vector>
#include <cmath>
#include <set>

#include "Triangulation.h"

using namespace std;


// Returns 1, -1, or 0 depending on whether the current point is on the left of, on the right of, or on the oriented straight-line (ab)
int Point::onLeft(Point *a, Point *b) const
{
	long long det = (b->x-a->x)*(y-a->y) - (x-a->x)*(b->y-a->y);
	return det > 0 ? 1 : det < 0 ? -1 : 0;
}

// Returns the N° of the edge that the current triangle shares with triangle a
// or -1 if the two triangles share no edge
int Triangle::edge(Triangle *t) const
{
		int i=2;
		while(i>=0 && neighbor[i] != t)
			i--;
		return i;
}

// Copy constructor
Triangulation::Triangulation(const Triangulation &T)
{
    LP.reserve(T.LP.size());
    addPoint(INT_MAX,INT_MAX);
    for (int i=1;i<T.LP.size();i++)
    {
        addPoint(0,0);
    }
    LT.reserve(T.LT.size());
    for (int i=0;i<T.LT.size();i++)
    {
        addTriangle(nullptr,nullptr,nullptr);
    }
    overwriteBy(T);
}

// Overwrites the current triangulation by T.
// The triangulations must have the same numbers of vertices and triangles
void Triangulation::overwriteBy(const Triangulation &T)
{
    area2=T.area2;

    for (int i=1;i<LP.size();i++)
       *(LP[i])=*(T.LP[i]);

    for (int i=0;i<LT.size();i++)
    {
       *(LT[i])=*(T.LT[i]);
        for (int j=0;j<3;j++)
        {
            LT[i]->neighbor[j]=LT[T.LT[i]->neighbor[j]->num];
            LT[i]->vertex[j]=LP[T.LT[i]->vertex[j]->num];
        }
    }
    for (int i=1;i<LP.size();i++)
       LP[i]->t=LT[T.LP[i]->t->num];
}

// Destructor
Triangulation::~Triangulation()
{
    deleteVertices();
    deleteTriangles();
}

// Deletes the vertices of the triangulation
void Triangulation::deleteVertices()
{
    for (int i=0;i<LP.size();i++)
        delete LP[i];
    LP.clear();
}

// Deletes the triangles of the triangulation
void Triangulation::deleteTriangles()
{
    for (int i=0;i<LT.size();i++)
        delete LT[i];
    LT.clear();
}

// Returns twice the standard deviation of the area of the colored triangles of the triangulation
double Triangulation::areaStdDeviation2() const
{
    double area_mean = 0.0;
    int nb_triangles = 0;
    for (int i=0; i < LT.size(); i++) {
        if (LT[i]->color != -1) {
            area_mean += LT[i]->area2();
            nb_triangles++;
        }
    }
    area_mean /= nb_triangles;

    double area_stdev = 0.0;
    for (int i=0; i < LT.size(); i++) {
        if (LT[i]->color != -1) {
            double tmp = area_mean - LT[i]->area2();
            area_stdev += tmp * tmp;
        }
    }
    return sqrt(area_stdev / nb_triangles);
}

// Stores in the triangulation twice the area of the Polygon
void Triangulation::initArea(const vector<int>& Polygon)
{
    area2 = (LP[Polygon.back()]->x + LP[Polygon[0]]->x) * (LP[Polygon[0]]->y - LP[Polygon.back()]->y);
    for (int i = 0; i < Polygon.size()-1; i++)
    {
        area2 += (LP[Polygon[i]]->x + LP[Polygon[i + 1]]->x) * (LP[Polygon[i+1]]->y - LP[Polygon[i]]->y);
    }
}

// Adds a new point (x,y) to the list of points of the triangulation.
// Returns a pointer to this point.
Point* Triangulation::addPoint(long long x, long long y)
{
	Point *p = new Point(x,y,LP.size());
	LP.push_back(p);
	return p;
}

// Adds a new triangle to the list of triangles of the triangulation.
// The points s0, s1, and s2 are supposed to already belong to the list of points
// of the triangulation and become respectively the vertices N°0, 1, and 2
// of the new triangle.
// Returns a pointer to the new triangle.
Triangle* Triangulation::addTriangle(Point *s0, Point *s1, Point *s2, int color)
{
	Triangle *t = new Triangle(s0,s1,s2,color,LT.size());
	LT.push_back(t);
	return t;
}

// Pastes edge N°i of triangle a on edge N°j of triangle b.
// a and b are supposed to be already present in the triangle list
// and their edges i and j are supposed to be "pastable"
void Triangulation::pasteTriangles(Triangle *a, int i, Triangle *b, int j)
{
	a->neighbor[i] = b;
	b->neighbor[j] = a;
}

// Initializes the vector of points from a file with the format given for the competition.
// The vector is sorted in lexicographic order with respect to (x,y).
// The point at index 0 is the point at infinity.
// Returns the number of (finite) points, or -1 if the file could not be opened.
long long Triangulation::initLP(const std::string &fileName)
{
    std::ifstream f(fileName);
    if (!f)
    {
        std::cout << "Cannot open file " << fileName << std::endl;
        return -1;
    }
    std::string ch,chNum;

    // first line contains a comment with the number of points between brackets
    std::getline(f,ch);
    int i=1;
    while(!isdigit(ch[i])||ch[i-1]!='(')
        i++;
    while(isdigit(ch[i]))
    {
        chNum.push_back(ch[i]);
        i++;
    }
    int n=std::stoi(chNum);

    LP.reserve(n+1);

    // second line contains a comment with the area of the convex hull of the set of points
    std::getline(f,ch);

    addPoint(INT_MAX,INT_MAX);  // point at infinity:

    // every line contains the index of a point (from 0 to n-1) and its x and y coordinates
    // the points are already sorted in lexicographic order
    for (int i=1;i<=n;i++)
    {
        long long x,y;
        int j;
        f>>j;
        f>>x;
        f>>y;
        addPoint(x,y);
    }

    return n;
}

// Saves the ordering of the vertices of the triangulated polygon in a file
void Triangulation::savePolygon(const string &fileName)
{
    vector<int> Polygon;
    extractPolygon(Polygon);

    ofstream f(fileName);
    for (int i=0; i<Polygon.size(); i++)
        f << Polygon[i]-1 << " ";
}

// Stores in Polygon the ordering of the vertices of the triangulated polygon
void Triangulation::extractPolygon(vector<int>& Polygon)
{
    Polygon.clear();
    Polygon.reserve(LP.size() - 1);
    Point* p = LP[1];
    do
    {
        Polygon.push_back(p->num);
        p = p->t->vertex[(p->e+2)%3];
    } while (p != LP[1]);
}

// Constructs in Polygon an initial polygonization of the point set LP
void Triangulation::initialPolygonization(vector<int> &Polygon) const
{
    Polygon.resize(LP.size()-1);
    Point *leftmostVertex = LP[1], *rightmostVertex = LP[LP.size()-1];
    Polygon[0] = 1;
    int i = 0, j = Polygon.size();
    for (int k=2; k<LP.size()-1; k++)
        if (LP[k]->onLeft(leftmostVertex,rightmostVertex)>=0)
            Polygon[--j] = k;
        else
            Polygon[++i] = k;
    Polygon[++i] = LP.size()-1;
}





///////////////////////////////////////////////////////////////////////////////////////////////////
////          Sweep algorithm to construct a colored triangulation of a given polygon          ////
///////////////////////////////////////////////////////////////////////////////////////////////////


// Constructs a colored triangulation of the Polygon.
// Polygon is an ordering of the points of LP (already initialized).
void Triangulation::triangulatePolygon(const std::vector<int>& Polygon)
{
    int n = LP.size();
    std::vector<int> Event(n);
    for (int i = 0; i < n - 1; i++)
        Event[Polygon[i]] = i;

    std::multiset<SweeplineSegment, bool(*)(const SweeplineSegment&, const SweeplineSegment&)> SweepLine(lessSweepline);
    std::vector<std::set<SweeplineSegment, bool(*)(const SweeplineSegment&, const SweeplineSegment&)>::iterator> L_It(n - 1);
    std::multiset<SweeplineSegment, bool(*)(const SweeplineSegment&, const SweeplineSegment&)>::iterator it, itDown, itUp;

    it = SweepLine.insert(SweeplineSegment{ new Point{LP[1]->x - 1,INT_MIN}, new Point{LP[n - 1]->x + 1,INT_MIN} });

    int current = Event[1],                                   // the vertex to be handled
        previous = current == 0 ? n - 2 : current - 1,        // its predecessor on the polygon
        next = current == n - 2 ? 0 : current + 1;            // its successor on the polygon

    itUp = SweepLine.insert({ LP[Polygon[current]], LP[Polygon[previous]],0 });
    itDown = SweepLine.insert({ LP[Polygon[current]], LP[Polygon[next]],1 });
    L_It[previous] = itUp;
    L_It[current] = itDown;



    const_cast<Point*&>(SweepLine.begin()->rightmost) = LP[Polygon[current]]; // should be const, but we are aware of what we are doing

    Triangle* tUp = nullptr, * tDown = nullptr;

    for (int i = 2; i < n; i++)
    {
        current = Event[i];                                   // the vertex to be handled
        previous = current == 0 ? n - 2 : current - 1;        // its predecessor on the polygon
        next = current == n - 2 ? 0 : current + 1;            // its successor on the polygon

        if (Polygon[previous] > Polygon[current] && Polygon[next] > Polygon[current])   // the segments [current,previous] and [current,next] have to be inserted in the SweepLine
        {
            itDown = SweepLine.insert({ LP[Polygon[current]], LP[Polygon[previous]] });
            itUp = SweepLine.insert({ LP[Polygon[current]], LP[Polygon[next]] });
            L_It[previous] = itDown;
            L_It[current] = itUp;
            if (std::prev(itUp) != itDown)
                std::swap(itDown, itUp);

            const_cast<int&>(itDown->color) = (std::prev(itDown)->color == 0 ? 1 : 0);
            const_cast<int&>(itUp->color) = std::prev(itDown)->color;

            beginApexEvent(prev(itDown)->succ, prev(itDown)->pred, prev(itDown)->rightmost, LP[Polygon[current]], tUp, tDown);
            tUp->color = tDown->color = itUp->color;

            const_cast<Point*&>(std::prev(itDown)->rightmost) = LP[Polygon[current]];
            const_cast<Triangle*&>(std::prev(itDown)->succ) = nullptr;
            const_cast<Triangle*&>(std::prev(itDown)->pred) = tDown;

            const_cast<Triangle*&>(itUp->succ) = tUp;
        }
        else
            if (Polygon[previous] < Polygon[current] && Polygon[next] < Polygon[current])   // the segments [previous,current] and [next,current] have to be removed from the SweepLine
            {
                itUp = L_It[current];
                itDown = L_It[previous];
                bool previousDown = true;
                if (std::prev(itUp) != itDown)
                {
                    std::swap(itDown, itUp);
                    previousDown = false;
                }
                endApexEvent(itUp->succ, itUp->pred, itUp->rightmost,
                    prev(itDown)->succ, prev(itDown)->pred, prev(itDown)->rightmost,
                    itDown->succ, itDown->pred, itDown->rightmost,
                    LP[Polygon[current]], tUp, tDown, previousDown);
                tUp->color = tDown->color = itUp->color;

                const_cast<Point*&>(std::prev(itDown)->rightmost) = LP[Polygon[current]];
                const_cast<Triangle*&>(std::prev(itDown)->succ) = tUp;
                const_cast<Triangle*&>(std::prev(itDown)->pred) = tDown;

                SweepLine.erase(itUp);
                SweepLine.erase(itDown);
            }
            else   // the one of the segments [previous,current] and [next,current] already in the SweepLine has to be replaced by the other
            {
                bool previousBefore;
                if (Polygon[previous] < Polygon[current])
                {
                    it = L_It[previous];
                    const_cast<Point*&>(it->first) = LP[Polygon[current]];
                    const_cast<Point*&>(it->last) = LP[Polygon[next]];
                    L_It[current] = it;
                    previousBefore = true;
                }
                else
                {
                    it = L_It[current];
                    const_cast<Point*&>(it->first) = LP[Polygon[current]];
                    const_cast<Point*&>(it->last) = LP[Polygon[previous]];
                    L_It[previous] = it;
                    previousBefore = false;
                }

                intermediateVertexEvent(it->succ, it->pred, it->rightmost,
                    prev(it)->succ, prev(it)->pred, prev(it)->rightmost,
                    LP[Polygon[current]], tUp, tDown, previousBefore);
                tUp->color = it->color;
                tDown->color = prev(it)->color;

                const_cast<Point*&>(std::prev(it)->rightmost) = LP[Polygon[current]];
                const_cast<Triangle*&>(std::prev(it)->succ) = nullptr;
                const_cast<Triangle*&>(std::prev(it)->pred) = tDown;

                const_cast<Point*&>(it->rightmost) = LP[Polygon[current]];
                const_cast<Triangle*&>(it->succ) = tUp;
                const_cast<Triangle*&>(it->pred) = nullptr;
            }
    }

    tUp->color = -1;
    while (tUp->neighbor[2] != nullptr)
    {
        tUp = tUp->neighbor[2];
        tUp->color = -1;
    }
    tDown->color = -1;
    while (tDown->neighbor[1] != nullptr)
    {
        tDown = tDown->neighbor[1];
        tDown->color = -1;
    }
    pasteTriangles(tUp, 2, tDown, 1);
}


// Starting from frontCurrent, traverses all sub-front triangles visible from rightmostVertex in CCW (resp., CW) order when order=2 (resp., 1).
// Constructs a new sub-front triangle after the last traversed visible sub-front triangle.
// frontCurrent may be nullptr and frontFirst will be the first triangle of the traversed sub-front.
// The vertex 0 of every traversed sub-front triangle is set to currentVertex.
void Triangulation::allVisibleFront(Triangle* frontCurrent, Point* rightmostVertex, Point* currentVertex, int order, Triangle*& frontFirst, Triangle*& frontNew)
{
    Triangle* lastVisible= nullptr, *firstInvisible = frontCurrent;
    if (firstInvisible == nullptr || currentVertex->onLeft(firstInvisible->vertex[2], firstInvisible->vertex[1]) >= 0)
    {
        if (order == 2)
            frontNew = addTriangle(LP[0], rightmostVertex, currentVertex);
        else
            frontNew = addTriangle(LP[0], currentVertex, rightmostVertex);
        frontFirst = frontNew;
    }
    else
    {
       frontFirst = firstInvisible;
        do
        {
            lastVisible = firstInvisible;
            lastVisible->vertex[0] = currentVertex;
            firstInvisible = lastVisible->neighbor[order];
        } while (firstInvisible != nullptr && currentVertex->onLeft(firstInvisible->vertex[2], firstInvisible->vertex[1]) < 0);
        if (order == 2)
            frontNew = addTriangle(LP[0], lastVisible->vertex[1], currentVertex);
        else
            frontNew = addTriangle(LP[0], currentVertex, lastVisible->vertex[2]);
        pasteTriangles(frontNew, 0, lastVisible, order);
    }
    if (firstInvisible != nullptr)
        pasteTriangles(frontNew, order, firstInvisible, order==2?1:2);
}


// Traverses all sub-front triangles starting from frontCurrent in CCW (resp., CW) order when order = 2 (resp., 1).
// frontCurrent is supposed to be non nullptr and frontLast will be the last traversed sub-front triangle
// The vertex 0 of every traversed sub-front triangle is set to currentVertex.
void Triangulation::allFront(Triangle* frontCurrent, Point* currentVertex, int order, Triangle*& frontLast)
{
    frontCurrent->vertex[0] = currentVertex;
    frontLast = frontCurrent;
    while (frontLast->neighbor[order] != nullptr)
    {
        frontLast = frontLast->neighbor[order];
        frontLast->vertex[0] = currentVertex;
    }
}

void Triangulation::beginApexEvent(Triangle* frontSucc, Triangle* frontPred, Point* rightmostVertex, Point* currentVertex, Triangle*& tUp, Triangle*& tDown)
{
    Triangle* lowestUp, *highestDown;
    allVisibleFront(frontSucc, rightmostVertex, currentVertex, 2, lowestUp, tUp);
    allVisibleFront(frontPred, rightmostVertex, currentVertex, 1, highestDown, tDown);
    pasteTriangles(lowestUp, lowestUp == tUp ? 0 : 1, highestDown, highestDown == tDown ? 0 : 2);
}

void Triangulation::endApexEvent(Triangle* frontUpperSucc, Triangle* frontUpperPred, Point* rightmostVertexUpper,
    Triangle* frontLowerSucc, Triangle* frontLowerPred, Point* rightmostVertexLower,
    Triangle* frontMiddleSucc, Triangle* frontMiddlePred, Point* rightmostVertexMiddle,
    Point* currentVertex, Triangle*& tUp, Triangle*& tDown, bool previousDown)
{
    Triangle* lowestUp;
    allVisibleFront(frontUpperSucc, rightmostVertexUpper, currentVertex, 2, lowestUp, tUp);
    if (frontUpperPred != nullptr)
    {
        if (lowestUp == tUp)
            pasteTriangles(lowestUp, 0, frontUpperPred, 2);
        allFront(frontUpperPred, currentVertex, 1, lowestUp);
    }

    Triangle* highestDown;
    allVisibleFront(frontLowerPred, rightmostVertexLower, currentVertex, 1, highestDown, tDown);
    if (frontLowerSucc != nullptr)
    {
        if (highestDown == tDown)
            pasteTriangles(highestDown, 0, frontLowerSucc, 1);
        allFront(frontLowerSucc, currentVertex, 2, highestDown);
    }

    Triangle* highestMiddle;
    if (frontMiddleSucc == nullptr)
        highestMiddle = frontMiddlePred;
    else
        allFront(frontMiddleSucc, currentVertex, 2, highestMiddle);
    pasteTriangles(highestMiddle, 2, lowestUp, lowestUp == tUp ? 0 : 1);

    Triangle* lowestMiddle;
    if (frontMiddlePred == nullptr)
        lowestMiddle = frontMiddleSucc;
    else
        allFront(frontMiddlePred, currentVertex, 1, lowestMiddle);
    pasteTriangles(lowestMiddle, 1, highestDown, highestDown == tDown ? 0:2);

    pasteTriangles(tUp, 1, tDown, 2);

    if (previousDown)
    {
        currentVertex->t = highestMiddle;
        currentVertex->e = 2;
        lowestMiddle->vertex[2]->t = lowestMiddle;
        lowestMiddle->vertex[2]->e = 1;
    }
    else
    {
        currentVertex->t = highestDown;
        currentVertex->e = highestDown == tDown ? 0 : 2;
        highestMiddle->vertex[1]->t = lowestUp;
        highestMiddle->vertex[1]->e = lowestUp == tUp ? 0 : 1;
    }
}

void Triangulation::intermediateVertexEvent(Triangle* frontUpperSucc, Triangle* frontUpperPred, Point* rightmostVertexUpper,
    Triangle* frontLowerSucc, Triangle* frontLowerPred, Point* rightmostVertexLower,
    Point* currentVertex, Triangle*& tUp, Triangle*& tDown, bool previousBefore)
{
    Triangle* lowestUp;
    allVisibleFront(frontUpperSucc, rightmostVertexUpper, currentVertex, 2, lowestUp, tUp);
    if (frontUpperPred != nullptr)
    {
        if (lowestUp == tUp)
            pasteTriangles(lowestUp, 0, frontUpperPred, 2);
        allFront(frontUpperPred, currentVertex, 1, lowestUp);
    }

    Triangle* highestDown;
    allVisibleFront(frontLowerPred, rightmostVertexLower, currentVertex, 1, highestDown, tDown);
    if (frontLowerSucc != nullptr)
    {
        if (highestDown == tDown)
            pasteTriangles(highestDown, 0, frontLowerSucc, 1);
        allFront(frontLowerSucc, currentVertex, 2, highestDown);
    }
    int lowestUpEdge = lowestUp == tUp ? 0 : 1, highestDownEdge = highestDown == tDown ? 0 : 2;
    pasteTriangles(lowestUp, lowestUpEdge, highestDown, highestDownEdge);
    if (previousBefore)
    {
        lowestUp->vertex[lowestUpEdge+1]->t = lowestUp;
        lowestUp->vertex[lowestUpEdge+1]->e = lowestUpEdge;
    }
    else
    {
        currentVertex->t = highestDown;
        currentVertex->e = highestDownEdge;
    }
}

// Returns 1, -1, or 0 depending on whether p is on the left of, on the right of, or on the oriented straight-line (s.first,s.last)
int left(const SweeplineSegment &s, const Point &p)
{
    long long det=(s.last->x-s.first->x)*(p.y-s.first->y) - (p.x-s.first->x)*(s.last->y-s.first->y);
	return det > 0 ? 1 : det < 0 ? -1 : 0;
}

// Returns true if s1 < s2 on the sweepline
bool lessSweepline(const SweeplineSegment &s1, const SweeplineSegment &s2)
{
    if (s1.first->x > s2.first->x)
         return left(s2, *s1.first)<0;
    else
        if (s2.first->x > s1.first->x)
            return left(s1, *s2.first)>0;
        else
            if (s2.first->y > s1.first->y)
                return true;
            else
                if (s2.first->y < s1.first->y)
                    return false;
                else
                    return left(s2, *s1.last)<0;
}


