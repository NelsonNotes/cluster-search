#pragma once

using namespace std;

class Point {
private:
	double x, y;
	long long int cluster_label;
public:
	Point() { // Default constructor
		x = 0;
		y = 0;
		cluster_label = -1;
	}
	Point(double _x, double _y) { // New point
		x = _x;
		y = _y;
		cluster_label = -1;
	}
	Point(double _x, double _y, long long int _cluster_label) { // Adding point to the cluster
		x = _x;
		y = _y;
		cluster_label = _cluster_label;
	}
	double GetX() {
		return x;
	}
	double GetY() {
		return y;
	}
	void SetX(double _x) {
		x = _x;
	}
	void SetY(double _y) {
		y = _y;
	}
	void SetLabel(long long int _label) {
		cluster_label = _label;
	}
	long long int GetLabel() {
		return cluster_label;
	}
	void MoveX(double dx) {
		x += dx;
	}
	void MoveY(double dy) {
		y += dy;
	}
	Point& operator=(const Point& point) {
		x = point.x;
		y = point.y;
		cluster_label = point.cluster_label;
		return *this;
	}
};