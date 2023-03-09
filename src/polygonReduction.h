#pragma once
#include "List.h"
#include "types.h"

class Triangle;
class Vertex;

class Triangle {
public:
	Vertex * vertex[3];
	point3d normal;
    Triangle(Vertex *v0, Vertex *v1, Vertex *v2);
	~Triangle();
	void ComputeNormal();
	void ReplaceVertex(Vertex *vold, Vertex *vnew);
	int HasVertex(Vertex *v);
    polygonIndexes getVertices(List<Vertex*>& vertices);
};

class Vertex {
public:
	point3d position;
    int id;
	List<Vertex *> neighbor;
	List<Triangle *> face;
	real objdist;
	Vertex * collapse;
	Vertex(point3d v, int id_);
    ~Vertex();
	void RemoveIfNonNeighbor(Vertex *n);
};

Triangle::Triangle(Vertex *v0,Vertex *v1,Vertex *v2){
	assert(v0!=v1 && v1!=v2 && v2!=v0);  //#mod1
	vertex[0]=v0;
	vertex[1]=v1;
	vertex[2]=v2;
	ComputeNormal();
	for(int i=0;i<3;i++) {
		vertex[i]->face.Add(this);
		for(int j=0;j<3;j++) if(i!=j) {
			vertex[i]->neighbor.AddUnique(vertex[j]);
		}
	}
}

Triangle::~Triangle(){
	int i;
	for(i=0;i<3;i++) {
		if(vertex[i]) vertex[i]->face.Remove(this);
	}
	for(i=0;i<3;i++) {
		int i2 = (i+1)%3;
		if(!vertex[i] || !vertex[i2]) continue;
		vertex[i]->RemoveIfNonNeighbor(vertex[i2]);
		vertex[i2]->RemoveIfNonNeighbor(vertex[i]);
	}
}

int Triangle::HasVertex(Vertex *v) {
	return (v==vertex[0] ||v==vertex[1] || v==vertex[2]);
}

void Triangle::ComputeNormal(){
	point3d v0=vertex[0]->position;
	point3d v1=vertex[1]->position;
	point3d v2=vertex[2]->position;
    point3d a = v1 - v0;
	normal = (v1-v0).crossProduct(v2-v1);
	if(magnitude(normal)==0)return;
	normal = normalizeNormal(normal);
}

void Triangle::ReplaceVertex(Vertex *vold,Vertex *vnew) {
	assert(vold && vnew);
	assert(vold==vertex[0] || vold==vertex[1] || vold==vertex[2]);
	assert(vnew!=vertex[0] && vnew!=vertex[1] && vnew!=vertex[2]);
	if(vold==vertex[0]){
		vertex[0]=vnew;
	}
	else if(vold==vertex[1]){
		vertex[1]=vnew;
	}
	else {
		assert(vold==vertex[2]);
		vertex[2]=vnew;
	}
	int i;
	vold->face.Remove(this);
	assert(!vnew->face.Contains(this));
	vnew->face.Add(this);
	for(i=0;i<3;i++) {
		vold->RemoveIfNonNeighbor(vertex[i]);
		vertex[i]->RemoveIfNonNeighbor(vold);
	}
	for(i=0;i<3;i++) {
		assert(vertex[i]->face.Contains(this)==1);
		for(int j=0;j<3;j++) if(i!=j) {
			vertex[i]->neighbor.AddUnique(vertex[j]);
		}
	}
	ComputeNormal();
}

polygonIndexes Triangle::getVertices(List<Vertex*>& vertices){
    return {(unsigned long)vertices.Find(vertex[0]), (unsigned long)vertices.Find(vertex[1]), (unsigned long)vertices.Find(vertex[2])};
}

Vertex::Vertex(point3d v, int id_) {
	position =v;
    id = id_;
}

Vertex::~Vertex(){
	assert(face.num==0);
	while(neighbor.num) {
		neighbor[0]->neighbor.Remove(this);
		neighbor.Remove(neighbor[0]);
	}
	// vertices.Remove(this);
}
void Vertex::RemoveIfNonNeighbor(Vertex *n) {
	if(!neighbor.Contains(n)) return;
	for(int i=0;i<face.num;i++) {
		if(face[i]->HasVertex(n)) return;
	}
	neighbor.Remove(n);
}

real magnitude(point3d v) {
    return (real)sqrt((v[0] * v[0]) + (v[1] * v[1])+ (v[2] * v[2]));
}

point3d normalizeNormal(point3d v) {
    float d=magnitude(v);
    if (d==0) {
		printf("Cant normalize ZERO vector\n");
		assert(0);
		d=0.1f;
	}
    v[0]/=d;
    v[1]/=d;
    v[2]/=d;
    return v;
}

float ComputeEdgeCollapseCost(Vertex *u,Vertex *v) {
	int i;
	float edgelength = magnitude(v->position - u->position);
	float curvature=0;

	List<Triangle *> sides;
	for(i=0;i<u->face.num;i++) {
		if(u->face[i]->HasVertex(v)){
			sides.Add(u->face[i]);
		}
	}
	for(i=0;i<u->face.num;i++) {
		float mincurv=1;
		for(int j=0;j<sides.num;j++) {
			point3d normal1 = u->face[i]->normal;
            point3d normal2 = sides[j]->normal;
            float dotprod = normal1[0]*normal2[0] + normal1[1]*normal2[1] + normal1[2]*normal2[2];
			mincurv = std::min(mincurv,(1-dotprod)/2.0f);
		}
		curvature = std::max(curvature,mincurv);
	}
	return edgelength * curvature;
}

void ComputeEdgeCostAtVertex(Vertex *v) {
	if(v->neighbor.num==0) {
		v->collapse=NULL;
		v->objdist=-0.01f;
		return;
	}
	v->objdist = 1000000;
	v->collapse=NULL;
	for(int i=0;i<v->neighbor.num;i++) {
		float dist;
		dist = ComputeEdgeCollapseCost(v,v->neighbor[i]);
		if(dist<v->objdist) {
			v->collapse=v->neighbor[i];
			v->objdist=dist;
		}
	}
}
void ComputeAllEdgeCollapseCosts(List<Vertex*>& vertices) {
	for(int i=0;i<vertices.num;i++) {
		ComputeEdgeCostAtVertex(vertices[i]);
	}
}

void Collapse(Vertex *u,Vertex *v, List<Vertex*>& vertices, List<Triangle*>& triangles){
    for(int i=0; i<u->face.num; i++){
        if(u->face.num <= i) continue;
        if(u->face[i]->HasVertex(v)){
            triangles.Remove(u->face[i]);
            delete u->face[i];
        }
        u->face[i]->ReplaceVertex(u, v);
    }
    vertices.Remove(u);
    delete u;
}

int getParentVertex(int org, std::map<int, int>& optimizedVertices, std::map<int, int>& pointsVerticesMapper){
    while(optimizedVertices[org] != -1){
        org = optimizedVertices[org];
    }
    return pointsVerticesMapper[org];
}

// Check if collapse is already collapsed to tail of curr
bool isOnTail(int curr, int collapse, std::map<int, int>& optimizedVetrices){
    int i = collapse;
    while(i != -1){
        i = optimizedVetrices[i];
        if(i == curr) return true;
    };
    return false;

}