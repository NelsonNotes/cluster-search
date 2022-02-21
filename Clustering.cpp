#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <queue>
#include <cmath>
#include <random>
#include <iomanip>
#include <chrono>
#include "Point.h"

using namespace std;

#define ERR_OPEN_FILE (-1)
#define ERR_INPUT_INT (-2)
#define ERR_INPUT (-3)

double FoolCheckDouble(bool negative = false) {
	string input;
	bool flag = true, flag_point = false;
	int j = 0;
	while (flag) {
		try {
			cin >> input;
			if (input[0] == '-' && negative) j = 1;
			if (j >= input.length()) throw ERR_INPUT;
			for (int i = j; i < input.length(); i++) {
				if (((int(input[i]) > 57 || int(input[i]) < 48) && int(input[i]) != 46) || (int(input[i]) == 46 && flag_point)) {
					throw ERR_INPUT;
				}
				if (int(input[i]) == 46 && !flag_point) {
					flag_point = true;
				}
			}
			flag = false;
		}
		catch (int a) {
			cout << "Enter correct number: ";
		}
	}
	return stod(input);
}
long long int FoolCheck(bool non_negative = false) {
	string input;
	bool flag = true;
	int j = 0;
	while (flag) {
		try {
			cin >> input;
			if (input[0] == '-' && !non_negative) j = 1;
			if (j >= input.length()) throw ERR_INPUT_INT;
			for (int i = j; i < input.length(); i++) {
				if (int(input[i]) > 57 || int(input[i]) < 48)
					throw ERR_INPUT_INT;
			}
			flag = false;
		}
		catch (int a) {
			cout << "Enter correct number: ";
		}
	}
	return stoll(input);
}

ostream& operator<<(ostream& out, Point& point) {
	return out << setprecision(5) << point.GetX() << "\t" << point.GetY() << "\t" << point.GetLabel() << endl;
}
class Cloud {
private:
	pair<double, double> x, y; // spreading
	Point center; // center of this cloud
	vector<Point> cloud; // list of points in this cloud
	int cloud_count;
public:
	Cloud() { // Default constructor
		x.first = 0;
		x.second = 0;
		y.first = 0;
		y.second = 0;
		center = Point(0, 0);
		cloud.resize(0);
		cloud_count = 0;
	}
	Cloud(double x_center, double y_center, double x_mn, double x_mx, double y_mn, double y_mx, int count) {
		x.first = x_mn;
		x.second = x_mx;
		y.first = y_mn;
		y.second = y_mx;
		center = Point(x_center, y_center);
		cloud_count = count;
	}
	vector<Point>& GetCloud() {
		return cloud;
	}
	int GetCloudCount() {
		return cloud_count;
	}
	vector<double> CreateNorm(size_t numb, double minval, double maxval) { //norm gistogramm
		vector<double> normal_array(numb); // Create a vector of predefined size
		if (minval > maxval) {
			swap(minval, maxval);
		}
		if (minval < 0) {
			swap(minval, maxval);
		}
		random_device rd{};
		mt19937 gen{ rd() };
		if (abs(maxval - minval) / 12 != 0) {
			normal_distribution<double> rng_machine{ (maxval + minval) / 2, abs(maxval - minval) / 12 }; // set normal distribution of points
			for (size_t i = 0; i < numb; ++i) {
				normal_array[i] = rng_machine(gen);
			}
		}
		else {
			for (size_t i = 0; i < numb; ++i) {
				normal_array[i] = maxval;
			}
		}
		return normal_array;
	}
	void ShiftX(double dx) {
		for (auto& point : cloud) {
			point.MoveX(dx);
		}
		center.MoveX(dx);
	}
	void ShiftY(double dy) {
		for (auto& point : cloud) {
			point.MoveY(dy);
		}
		center.MoveY(dy);
	}
	double MoveXonAngle(Point point, double angle, double _x, double _y) {
		return _x + (point.GetX() - _x) * cos(angle) - (point.GetY() - _y) * sin(angle);
	}
	double MoveYonAngle(Point point, double angle, double _x, double _y) {
		return _y + (point.GetX() - _x) * sin(angle) + (point.GetY() - _y) * cos(angle);
	}
	void Turn(int angle, double _x = 0, double _y = 0) {
		for (auto& point : cloud) {
			double new_x = MoveXonAngle(point, angle, _x, _y);
			double new_y = MoveYonAngle(point, angle, _x, _y);
			point.SetX(new_x);
			point.SetY(new_y);
		}
		center.SetX(MoveXonAngle(center, angle, _x, _y));
		center.SetY(MoveXonAngle(center, angle, _x, _y));
	}
	void Filling(size_t count, Point& start, Point& end) {
		cloud.resize(count);
		center = Point((end.GetX() - start.GetX()) / 2, (end.GetY() - start.GetY()) / 2);
		vector <double> x_arr = CreateNorm(count, start.GetX(), end.GetX());
		vector <double> y_arr = CreateNorm(count, start.GetY(), end.GetY());
		for (size_t i = 0; i < count; ++i) {
			cloud[i].SetX(x_arr[i]);
			cloud[i].SetY(y_arr[i]);
			cloud[i].SetLabel(0);
			cloud_count++;
		}
	}
	~Cloud() {
		cloud.clear();
	}
};
class Field {
private:
	vector<vector<Point>> clouds; // all clouds in this field
	vector<Point> points; // all points in this field
	long long int clouds_count, points_count;
public:
	Field() {
		points.resize(0);
		clouds.resize(0);
		clouds_count = 0;
		points_count = 0;
	}
	long long int GetCloudsCount() {
		return clouds_count; // quantity of clouds
	}
	long long int GetPointsCount() {
		return points_count;
	}
	void SetCloudsCount(long long int _clouds_count) {
		clouds_count = _clouds_count;
	}
	void AddPointsCount(long long int _points_count) {
		points_count += _points_count;
	}
	vector<Point>& GetFieldPoints() {
		return points;
	}
	vector<vector<Point>>& GetClouds() {
		return clouds;
	}
	void AddCloud() {
		Cloud cloud;
		clouds.push_back(cloud.GetCloud());
		clouds_count++;
	}
	int GetFieldFromFile(string filename) { // read from file, file is set by user
		ifstream fin(filename);
		double x, y;
		long long int l;
		if (!fin.is_open()) {
			cout << "Couldn't open the file" << endl;
			return ERR_OPEN_FILE;
		}
		else {
			while (!fin.eof()) {
				fin >> x >> y >> l;
				points.push_back(Point(x, y, l));
				points_count++;
				if (l > clouds_count) clouds_count = l;
			}
			fin.close();
		}
		PrintField();
		return 0;
	}
	void PrintGNU() {
		fstream script("output.plt", ios::out | ios::trunc); // output to gnuplot
		if (!script.is_open()) {
			cout << "Couldn't open the file" << endl;
			exit(ERR_OPEN_FILE);
		}
		script << "plot 'data.txt' u 1:2:3 with points lc variable title 'work'" << endl;
		script.close();
	}
	void PrintField() { // print to file
		ofstream fout("data.txt");
		if (!fout.is_open()) {
			cout << "Couldn't open the file" << endl;
			exit(ERR_OPEN_FILE);
		}
		for (auto& cloud : clouds) { // print to file
			for (auto& point : cloud) {
				fout << point.GetX() << " " << point.GetY() << " " << point.GetLabel() << endl;
			}
		}
		fout.close();
		PrintGNU();
	}
	~Field() {
		points.clear();
		clouds.clear();
		points_count = 0;
		clouds_count = 0;
	}
};
class Wave {
private:
	vector<pair<Point, int>> points; // all points and their marks(not reached, in queue, visited(current cluster), written)
	queue<pair<Point*, int*>> work_points; // neighbours
	long long int label = 1; // current cluster label
	long long int not_visited; // quantity of not visited points
public:
	Wave() {
		points.resize(0);
		work_points;
		label = 1;
		not_visited = 0;
	}
	Wave(vector<vector<Point>>& clouds) {
		for (size_t j = 0; j < clouds.size(); ++j) {
			vector<Point> cloud = clouds[j];
			for (size_t i = 0; i < cloud.size(); ++i) {
				points.push_back(make_pair(cloud[i], 0));
			}
		}
		not_visited = points.size();
	}
	double Distance(Point& first, Point& second) {
		double x = first.GetX() - second.GetX();
		double y = first.GetY() - second.GetY();
		return sqrt((x * x + y * y));
	}
	void FindNotVisited() {
		for (size_t i = 0; i < points.size(); ++i) {
			if (points[i].second == 0) {
				work_points.push(make_pair(&points[i].first, &points[i].second));
				break;
			}
		}
	}
	void MarkCluster() { // mark found cluster
		for (auto& p : points) {
			if (p.second == 2) {
				p.second = 3;
				p.first.SetLabel(label);
				--not_visited;
			}
		}
		++label;
	}
	void Search(const double EPS) {
		while (not_visited > 0) {
			if (not_visited != 0 && work_points.empty()) { // cluster was formed, but there is points left in the field
				FindNotVisited();
			}
			while (!work_points.empty()) { // search algorithm
				pair<Point*, int*> p = work_points.front();
				work_points.pop();
				*p.second = 2; // mark point as visited(current cluster)
				for (size_t i = 0; i < points.size(); ++i) {
					pair<Point*, int*> to = make_pair(&points[i].first, &points[i].second);
					if (*to.second == 0 && Distance(*p.first, *to.first) <= EPS) {
						*to.second = 1; // mark point(in queue)
						work_points.push(to); // add to queue
					}
				}
			}
			MarkCluster();
		}
	}
	void Print() {
		ofstream fout("data.txt");
		if (!fout.is_open()) {
			cout << "Couldn't open the file" << endl;
			exit(ERR_OPEN_FILE);
		}
		vector<long long> colors;
		long long int color_count = 0;
		for (auto& point : points) {
			if (point.first.GetLabel() > color_count) {
				color_count = point.first.GetLabel();
			}
		}
		colors.resize(color_count + 1, 0);
		srand(time(NULL));
		for (size_t i = 0; i < color_count; ++i) {
			if (colors[i] == 0) {
				colors[i] = 65536 * rand() + 256 * rand() + rand();
			}
		}
		for (auto& point : points) { // print to file
			fout << point.first.GetX() << " " << point.first.GetY() << " " << colors[point.first.GetLabel()] << endl;
		}
		fout.close();
		fstream script("output.plt", ios::out | ios::trunc); // output to gnuplot
		if (!script.is_open()) {
			cout << "Couldn't open the file" << endl;
			exit(ERR_OPEN_FILE);
		}
		script << "plot 'data.txt' u 1:2:3 with points lc variable title 'wave'" << endl;
		script.close();
	}
};
class Dbscan {
private:
	vector<pair<Point, int>> points; // all points and their marks(not reached, in queue, visited(current cluster), written)
	queue<pair<Point*, int*>> work_points; // neighbours
	long long int label = 1; // current cluster label
	long long int not_visited; // quantity of not visited points
	long long int noise = 0; // label for lonely points
public:
	Dbscan() {
		points.resize(0);
		work_points;
		label = 1;
		not_visited = 0;
	}
	Dbscan(vector<vector<Point>>& clouds) {
		for (size_t j = 0; j < clouds.size(); ++j) {
			vector<Point> cloud = clouds[j];
			for (size_t i = 0; i < cloud.size(); ++i) {
				points.push_back(make_pair(cloud[i], 0));
			}
		}
		not_visited = points.size();
	}
	double Distance(Point& first, Point& second) {
		double x = first.GetX() - second.GetX();
		double y = first.GetY() - second.GetY();
		return sqrt((abs(x * x) + abs(y * y)));
	}
	void NewCluster() {
		for (size_t i = 0; i < points.size(); ++i) {
			if (points[i].second == 0) {
				points[i].second = 2;
				work_points.push(make_pair(&points[i].first, &points[i].second));
				not_visited--;
				break;
			}
		}
	}
	void MarkCluster(const long long int min_points, const long long int k) { // mark new cluster
		for (auto& p : points) {
			if (p.second == 2 && k >= min_points) {
				p.second = 3;
				p.first.SetLabel(label);
			}
			else if (p.second == 2 && k < min_points) { // lonely points(do not belong to any cluster), let them be, they're as happy as any of us
				p.second = 3;
				p.first.SetLabel(noise);
			}
		}
		if (k >= min_points) ++label; // if it was a cluster and not just noise(lonely point)
	}
	void Search(const double EPS, const long long int min_points) {
		while (not_visited > 0) {
			if (not_visited != 0 && work_points.empty()) { // cluster was formed, but there is points left
				NewCluster();
			}
			long long int k = 0; // how many points in this cluster(if enough to make a cluster)
			while (!work_points.empty()) { // algorithm
				pair<Point*, int*> p = work_points.front();
				work_points.pop();
				*p.second = 2;
				++k;
				for (size_t i = 0; i < points.size(); ++i) {
					pair<Point*, int*> to = make_pair(&points[i].first, &points[i].second);
					if (*to.second == 0 && Distance(*p.first, *to.first) <= EPS) {
						*to.second = 1;
						work_points.push(to);
						--not_visited;
					}
				}
			}
			MarkCluster(min_points, k);
		}
	}
	void Print() {
		ofstream fout("data.txt");
		if (!fout.is_open()) {
			cout << "Couldn't open the file" << endl;
			exit(ERR_OPEN_FILE);
		}
		vector<long long> colors;
		long long int color_count = 0;
		for (auto& point : points) {
			if (point.first.GetLabel() > color_count) {
				color_count = point.first.GetLabel();
			}
		}
		colors.resize(color_count + 1, 0);
		srand(time(NULL));
		for (size_t i = 0; i < color_count; ++i) {
			if (colors[i] == 0) {
				colors[i] = 65536 * rand() + 256 * rand() + rand();
			}
		}
		for (auto& point : points) { // print to file
			fout << point.first.GetX() << " " << point.first.GetY() << " " << colors[point.first.GetLabel()] << endl;
		}
		fout.close();
		fstream script("output.plt", ios::out | ios::trunc); // output to gnuplot
		if (!script.is_open()) {
			cout << "Couldn't open the file" << endl;
			exit(ERR_OPEN_FILE);
		}
		script << "plot 'data.txt' u 1:2:3 with points lc variable title 'dbscan'" << endl;
		script.close();
	}
};
class Hierarchy {
private:
	vector<Point> points; // all points in the field
	vector<vector<double>> distances; // distances between clusters
	vector<vector<Point>> clusters; // found clusters
public:
	Hierarchy() {
		points.resize(0);
		distances.resize(0);
		clusters.resize(0);
	}
	double Distance(Point& first, Point& second) {
		double x = first.GetX() - second.GetX();
		double y = first.GetY() - second.GetY();
		return sqrt((x * x + y * y));
	}
	Hierarchy(vector<vector<Point>>& clouds) {
		size_t label = 1, i = 0;
		size_t size_of_clusters = 0;
		for (auto& cloud : clouds) {
			size_of_clusters += cloud.size();
		}
		clusters.resize(size_of_clusters);
		for (auto& cloud : clouds) {
			for (auto& point : cloud) {
				point.SetLabel(label); // each point is a cluster by the start
				points.push_back(point);
				clusters[i].push_back(point);
				++label;
				++i;
			}
		}
		size_t points_size = points.size();
		for (size_t i = 0; i < points_size; ++i) {
			distances.resize(points_size);
			for (size_t j = 0; j < points_size; ++j) {
				distances[i].push_back(Distance(points[i], points[j]));
			}
		}
	}
	pair<long long int, long long int> FindMinDistance() {
		long long int mn = -1;
		long long int res_i = -1, res_j = -1;
		for (size_t i = 0; i < distances.size(); ++i) {
			for (size_t j = 0; j < distances[i].size(); ++j) {
				if ((distances[i][j] < mn || mn < 0) && i != j) { // set mn as -1, because there can't be negative distance
					mn = distances[i][j];
					res_i = i;
					res_j = j;
				}
			}
		}
		return make_pair(res_i, res_j);
	}
	void CreateNewDistances(pair<long long int, long long int> united) {
		vector<vector<double>> new_distances = distances;
		for (size_t i = 0; i < distances.size(); ++i) {
			if (i != united.first) { // formulas for distance
				distances[united.first][i] = (1 / 2) * distances[united.first][i] + (1 / 2) * distances[united.second][i] - (1 / 2) * abs(distances[united.first][i] - distances[united.second][i]);
				distances[i][united.first] = 1 / 2 * distances[united.first][i] + 1 / 2 * distances[united.second][i] - 1 / 2 * abs(distances[united.first][i] - distances[united.second][i]);
			}
		}
		for (size_t i = 0; i < distances.size(); ++i) { // delete one of the united clusters
			new_distances[i].erase(new_distances[i].begin() + united.second);
			if (i >= new_distances.size()) {
				break;
			}
		}
		new_distances.erase(new_distances.begin() + united.second);
		distances = new_distances;
	}
	void Search(const long long int numb_of_clusters) {
		while (distances.size() != numb_of_clusters) { // algorithm
			pair<long long int, long long int> united = FindMinDistance(); // find min distance between two clusters
			CreateNewDistances(united); // unite two min-distance clusters, update distances
			for (size_t i = 0; i < clusters[united.second].size(); ++i) {
				clusters[united.first].push_back(clusters[united.second][i]); // continue to unite two min-distance clusters, update clusters
			}
			clusters.erase(clusters.begin() + united.second); // delete one of the united clusters
		}
		long long int label = 1;
		for (auto& cluster : clusters) { // set labels in right order
			for (auto& point : cluster) {
				point.SetLabel(label);
			}
			++label;
		}
	}
	void Print() {
		ofstream fout("data.txt");
		if (!fout.is_open()) {
			cout << "Couldn't open the file" << endl;
			exit(ERR_OPEN_FILE);
		}
		vector<long long> colors;
		long long color_count = 0;
		for (auto& cloud : clusters) {
			for (auto& point : cloud) {
				if (point.GetLabel() > color_count) {
					color_count = point.GetLabel();
				}
			}
		}
		colors.resize(color_count + 1, 0);
		srand(time(NULL));
		for (size_t i = 0; i < color_count; ++i) {
			if (colors[i] == 0) {
				colors[i] = 65536 * rand() + 256 * rand() + rand();
			}
		}
		for (auto& points : clusters) { // print to file
			for (auto& point : points) {
				fout << point.GetX() << " " << point.GetY() << " " << colors[point.GetLabel()] << endl;
			}
		}
		fout.close();
		fstream script("output.plt", ios::out | ios::trunc); // output to gnuplot
		if (!script.is_open()) {
			cout << "Couldn't open the file" << endl;
			exit(ERR_OPEN_FILE);
		}
		script << "plot 'data.txt' u 1:2:3 with points lc variable title 'hierarchy'" << endl;
		script.close();
	}
};
class Interface {
private:
	string command;
	string work_com;
	string filename;
public:
	void Control() {
		Field field;
		Point point;
		Cloud cloud;
		CommandList();
		cout << "Enter command: " << endl;
		do {
			bool mistake = false;
			getline(cin, command);
			bool flag = true;
			if (command == "name") {
				do {
					cout << "Enter name of file: " << endl;
					getline(cin, filename);
					int n = filename.length();
					string format = ".txt";
					flag = true;
					if (n > 3) {
						for (int i = n - 4; i < n; i++) {
							if (filename[i] != format[i - n + 4]) {
								flag = false;
								break;
							}
						}
						if (!flag) {
							cout << "Not valid name of file" << endl;
							mistake = true;
						}
						else if (field.GetFieldFromFile(filename)) {
							mistake = true;
							flag = false;
						}
						else {
							cout << "The field have been initialized!" << endl;
							flag = true;
						}
					}
				} while (!flag);
			}
			else if (command == "help") {
				CommandList();
			}
			else if (command == "create") {
				long long int clouds_count = 0;
				Point first;
				Point second;
				string command;
				while (true) {
					cout << "Enter x_start, y_start, x_end, y_end: " << endl;
					first.SetX(FoolCheckDouble(true));
					first.SetY(FoolCheckDouble(true));
					second.SetX(FoolCheckDouble(true));
					second.SetY(FoolCheckDouble(true));
					cout << "Enter quantity of points: " << endl;
					long long int count = FoolCheck();
					cloud.Filling(count, first, second);
					field.GetClouds().push_back(cloud.GetCloud());
					if (cloud.GetCloud().size()) {
						++clouds_count;
					}
					field.SetCloudsCount(clouds_count);
					field.PrintField();
					if (!field.GetCloudsCount() || !cloud.GetCloud().size()) {
						cout << "You haven't initialized the field" << endl;
						mistake = true;
					}
					else {
						field.AddPointsCount(cloud.GetCloud().size());
						cin.ignore();
						WorkList();
						cout << "Enter work command: " << endl;
						do {
							getline(cin, work_com);
							if (work_com == "help") {
								WorkList();
							}
							else if (work_com == "turn") {
								cout << "Enter angle: " << endl;
								int angle = FoolCheck(true);
								cout << "Enter x: " << endl;
								long long int x = FoolCheck();
								cout << "Enter y: " << endl;
								long long int y = FoolCheck();
								cloud.Turn(angle, x, y);
								field.GetClouds().pop_back();
								field.GetClouds().push_back(cloud.GetCloud());
								field.PrintField();
								cout << "The cloud has been turned by " << angle << " degrees around (" << x << ", " << y << ")" << endl;
							}
							else if (work_com == "shift x") {
								cout << "Enter delta x: " << endl;
								long long int dx = FoolCheck();
								cloud.ShiftX(dx);
								field.GetClouds().pop_back();
								field.GetClouds().push_back(cloud.GetCloud());
								field.PrintField();
								cout << "The cloud has been shift by " << dx << endl;
							}
							else if (work_com == "shift y") {
								cout << "Enter delta y: " << endl;
								long long int dy = FoolCheck();
								cloud.ShiftY(dy);
								field.GetClouds().pop_back();
								field.GetClouds().push_back(cloud.GetCloud());
								field.PrintField();
								cout << "The cloud has been shift by " << dy << endl;
							}
							else if (work_com != "exit") {
								if (work_com != "") cout << "Not valid work command (enter help for list of commands)" << endl;
								else cout << "Enter work command: " << endl;
							}
						} while (work_com != "exit");
						break;
					}
				}
			}
			else if (command == "wave" || command == "hierarchy" || command == "dbscan") {
				if (!field.GetCloudsCount()) {
					cout << "You haven't initialized the field yet." << endl;
					mistake = true;
				}
				else {
					if (command == "wave") {
						cout << "Set EPS: " << endl;
						double EPS;
						EPS = FoolCheckDouble();
						Wave wave(field.GetClouds());
						wave.Search(EPS);
						wave.Print();
						cin.ignore();
					}
					if (command == "hierarchy") {
						long long int numb_of_clusters = -1;
						do {
							if (numb_of_clusters > field.GetPointsCount()) {
								cout << "There is not enough points for so many clusters. ";
							}
							if (numb_of_clusters == 0) {
								cout << "You have set zero clusters. ";
							}
							cout << "Set quantity of clusters: " << endl;
							numb_of_clusters = FoolCheck(true);
						} while (numb_of_clusters > field.GetPointsCount() || numb_of_clusters == 0);
						Hierarchy hierarchy(field.GetClouds());
						hierarchy.Search(numb_of_clusters);
						hierarchy.Print();
						cin.ignore();
					}
					if (command == "dbscan") {
						double EPS;
						long long int min_points;
						cout << "Set EPS and minimum points in cluster: " << endl;
						EPS = FoolCheckDouble();
						min_points = FoolCheck(true);
						Dbscan dbscan(field.GetClouds());
						dbscan.Search(EPS, min_points);
						dbscan.Print();
						cin.ignore();
					}
				}
			}
			else if (command != "exit") {
				cout << "Not valid command (enter help for list of commands)" << endl;
			}
			if (!mistake && command != "exit") cout << "Enter command: " << endl;
		} while (command != "exit");
	}
	void CommandList() {
		cout << "List of command: \n"
			"help - print list of commands\n"
			"name - enter name of file\n"
			"create - create clusters\n"
			"wave - run wave search\n"
			"dbscan - run dbscan search\n"
			"hierarchy - run hierarchy search\n"
			"exit" << endl;
	}
	void WorkList() {
		cout << "List of work commands: \n"
			"help - print list of work commands\n"
			"turn - turn cloud at the angle\n"
			"shift x - move cloud by x\n"
			"shift y - move cloud by y\n"
			"exit" << endl;
	}
};
int main() {
	setlocale(LC_ALL, "russian");
	Interface interface;
	interface.Control();
	return 0;
}