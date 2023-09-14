#include "stdafx.h"
#include<list>
#include<fstream>
#include<iostream>
#include<string>
#include<vector>
#include<new>
#include "vizi.h"
using namespace std;

class point
{
public:
	double x;
	double y;
	int n;
	bool log;

	point() {

	}
};

double min(double a, double b, double c, double d)
{
	double res = a;
	if (b < res) res = b;
	if (c < res) res = c;
	if (d < res) res = d;
	return res;
}

double max(vector<double> v)
{
	double res = v[0];
	for (int i = 1; i < v.size(); i++)
		if (v[i] < res)
			res = v[i];
	return res;
}

list<point> read(string name) {
	list<point> cont;
	FILE* file = fopen(("../input/" + name).c_str(), "r");
	// save n
	int step; // step scale
	fscanf(file, "%d", &step);
	if (step < 1) step = 1;
	//cont.push_back(p);
	while (!feof(file)) {
		point pnt;
		fscanf(file, "%lf %lf %d %d\n", &pnt.x, &pnt.y, &pnt.n, &pnt.log);
		cont.push_back(pnt);
	}
	fclose(file);
	list<point>::iterator p = cont.begin();
	while (p != cont.end())
	{
		(*p).n *= step;
		p++;
	}

	return cont;
}

void save(vector<double>*** tab, vector<int> sizes, string name)
{
	vector<vector<array<double, 2>>> grid(sizes[0], vector<array<double, 2>>(sizes[1]));
	for (int i = 0; i < sizes[0]; i++) {
		for (int j = 0; j < sizes[1]; j++) {
			grid[i][j][0] = (*tab)[i][j][0];
			grid[i][j][1] = (*tab)[i][j][1];
		}
	}
	vizi::SaveRegGridToVizi(name, grid);
	vector<vector<double>> nodes_row(sizes[0], vector<double>(sizes[1]));
	for (int i = 0; i < sizes[0]; i++) {
		for (int j = 0; j < sizes[1]; j++) {
			nodes_row[i][j] = i;
		}
	}
	vizi::SaveRegDataToVizi(name, nodes_row, "nodes_row", 0, 1);
	vector<vector<double>> nodes_column(sizes[0], vector<double>(sizes[1]));
	for (int i = 0; i < sizes[0]; i++) {
		for (int j = 0; j < sizes[1]; j++) {
			nodes_column[i][j] = j;
		}
	}
	vizi::SaveRegDataToVizi(name, nodes_column, "nodes_column", 0, 1);

	vector<vector<double>> cells_row(sizes[0] - 1, vector<double>(sizes[1] - 1));
	for (int i = 0; i < sizes[0] - 1; i++) {
		for (int j = 0; j < sizes[1] - 1; j++) {
			cells_row[i][j] = i;
		}
	}
	vizi::SaveRegDataToVizi(name, cells_row, "cells_row", 0, 0);
	vector<vector<double>> cells_column(sizes[0] - 1, vector<double>(sizes[1] - 1));
	for (int i = 0; i < sizes[0] - 1; i++) {
		for (int j = 0; j < sizes[1] - 1; j++) {
			cells_column[i][j] = j;
		}
	}
	vizi::SaveRegDataToVizi(name, cells_column, "cells_column", 0, 0);
}

vector<int> getmatsize(list<point> lst)
{
	//nodes pointers
	list<point>::iterator p1 = lst.begin();
	list<point>::iterator p2;
	list<point>::iterator p3;

	vector<int> res;

	// fist node
	while (!(*p1).log)
	{
		p1++;
	}
	p2 = p1;
	p2++;
	// 2nd node
	while (!(*p2).log)
	{
		p2++;
	}
	p3 = p2;
	p3++;
	// 3rd node
	while (!(*p3).log)
	{
		p3++;
	}

	// count sizes
	int size1 = 1;
	int size2 = 1;
	while (p1 != p2)
	{
		size1 += (*p1).n;
		p1++;
	}
	while (p2 != p3)
	{
		size2 += (*p2).n;
		p2++;
	}
	res.push_back(size2);
	res.push_back(size1);
	return res;
}

void interpolate(vector<double>*** tab, list<point> lst, vector<int> sizes)
{
	// s
	vector<double> st = { 0 }; // top
	vector<double> sb = { 0 }; // bottom
	vector<double> sl = { 0 }; // left
	vector<double> sr = { 0 }; // right
	/*~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ inner points ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~*/
	// Top & Bottom
	double xbp = (*tab)[0][0][0];
	double ybp = (*tab)[0][0][1];
	double xtp = (*tab)[sizes[0] - 1][0][0];
	double ytp = (*tab)[sizes[0] - 1][0][1];
	for (int j = 1; j < sizes[1]; j++)
	{
		// Bottom
		double sbt = sqrt(pow((*tab)[0][j][0] - xbp, 2) + pow((*tab)[0][j][1] - ybp, 2));
		sbt += sb[j - 1];
		xbp = (*tab)[0][j][0];
		ybp = (*tab)[0][j][1];
		sb.push_back(sbt);

		// Top
		double stp = sqrt(pow((*tab)[sizes[0] - 1][j][0] - xtp, 2) + pow((*tab)[sizes[0] - 1][j][1] - ytp, 2));
		stp += st[j - 1];
		xtp = (*tab)[0][j][0];
		ytp = (*tab)[0][j][1];
		st.push_back(stp);
	}
	// Norm
	for (int j = 1; j < sizes[1] - 1; j++)
	{
		st[j] /= st[sizes[1] - 1];
		sb[j] /= sb[sizes[1] - 1];
	}

	// Left & Right
	double xlp = (*tab)[0][0][0];
	double ylp = (*tab)[0][0][1];
	double xrp = (*tab)[0][sizes[1] - 1][0];
	double yrp = (*tab)[0][sizes[1] - 1][1];
	for (int i = 1; i < sizes[0]; i++)
	{
		// Left
		double slft = sqrt(pow((*tab)[i][0][0] - xlp, 2) + pow((*tab)[i][0][1] - ylp, 2));
		xlp = (*tab)[i][0][0];
		ylp = (*tab)[i][0][0];
		slft += sl[i - 1];
		sl.push_back(slft);

		// Right
		double srt = sqrt(pow((*tab)[i][sizes[1] - 1][0] - xrp, 2) + pow((*tab)[i][sizes[1] - 1][1] - yrp, 2));
		yrp = (*tab)[i][sizes[1] - 1][1];
		xrp = (*tab)[i][sizes[1] - 1][0];
		srt += sr[i - 1];
		sr.push_back(srt);
	}
	// Norm
	for (int i = 1; i < sizes[0] - 1; i++)
	{
		sl[i] /= sl[sizes[0] - 1];
		sr[i] /= sr[sizes[0] - 1];
	}

	// Making ksi and eta arrays
	double** ksi;
	double** eta;
	try {
		ksi = new double* [sizes[0]];
		eta = new double* [sizes[0]];
		for (int i = 0; i < sizes[0]; i++)
		{
			ksi[i] = new double[sizes[1]];
			eta[i] = new double[sizes[1]];
			for (int j = 0; j < sizes[1]; j++)
			{
				ksi[i][j] = 0;
				eta[i][j] = 0;
			}
		}
	}
	catch (bad_alloc xa) {
		return;
	}

	// Filling ksi & eta arrays
	for (int i = 1; i < sizes[0] - 1; i++)
		for (int j = 1; j < sizes[1] - 1; j++)
		{
			ksi[i][j] = sb[j] + (st[j] - sb[j]) * sl[i];
			ksi[i][j] /= 1 - (st[j] - sb[j]) * (sr[i] - sl[i]);
			eta[i][j] = sl[i] + (sr[i] - sl[i]) * sb[j];
			eta[i][j] /= 1 - (st[j] - sb[j]) * (sr[i] - sl[i]);
		}
	//print(ksi, sizes);

	// Points making
	for (int i = 1; i < sizes[0] - 1; i++)
		for (int j = 1; j < sizes[1] - 1; j++)
		{
			(*tab)[i][j][0] = (*tab)[0][j][0] * (1 - eta[i][j]) + (*tab)[sizes[0] - 1][j][0] * eta[i][j]
				+ ksi[i][j] * ((*tab)[i][sizes[1] - 1][0] - (*tab)[0][sizes[1] - 1][0] * (1 - sr[i]) - (*tab)[sizes[0] - 1][sizes[1] - 1][0] * sr[i])
				+ (1 - ksi[i][j]) * ((*tab)[i][0][0] - (*tab)[0][0][0] * (1 - sl[i]) - (*tab)[sizes[0] - 1][0][0] * sl[i])
				;
			(*tab)[i][j][1] = (*tab)[0][j][1] * (1 - eta[i][j]) + (*tab)[sizes[0] - 1][j][1] * eta[i][j]
				+ ksi[i][j] * ((*tab)[i][sizes[1] - 1][1] - (*tab)[0][sizes[1] - 1][1] * (1 - sr[i]) - (*tab)[sizes[0] - 1][sizes[1] - 1][1] * sr[i])
				+ (1 - ksi[i][j]) * ((*tab)[i][0][1] - (*tab)[0][0][1] * (1 - sl[i]) - (*tab)[sizes[0] - 1][0][1] * sl[i])
				;
		}
}

void distribute(vector<double>*** tab, list<point> lst, vector<int> sizes)
{
	// temporary list of middle points
	//list<point> subpoints;
	list<vector<double>> subpoints;
	int kn = 0; // number of key points in da list
	list<point>::iterator p = lst.begin();
	list<point>::iterator s = p;
	s++;
	lst.push_back(*p);
	while (s != lst.end())
	{
		// main points


		// ranges number
		int rn = (*p).n;
		// distance between the points
		double dx = ((*s).x - (*p).x) / rn;
		double dy = ((*s).y - (*p).y) / rn;
		// Add points
		for (int i = 1; i < rn; i++)
		{
			vector<double> v;
			double x = (*p).x + dx * i;
			double y = (*p).y + dy * i;
			v.push_back(x);
			v.push_back(y);
			subpoints.push_back(v);
		}
		if (!(*s).log)
		{
			vector<double> v;
			v.push_back((*s).x);
			v.push_back((*s).y);
			subpoints.push_back(v);
		}
		p++;
		s++;
	}
	// Writing the frame
	int ch = 0;
	list<vector<double>>::iterator t = subpoints.begin();
	// top
	for (int i = 1; i < sizes[1] - 1; i++)
	{
		ch++;
		(*tab)[0][i] = *t;
		t++;
	}
	// right
	for (int i = 1; i < sizes[0] - 1; i++)
	{
		ch++;
		(*tab)[i][sizes[1] - 1] = *t;
		t++;
	}
	// bottom
	for (int i = sizes[1] - 2; i > 0; i--)
	{
		ch++;
		(*tab)[sizes[0] - 1][i] = *t;
		t++;
	}
	//left
	for (int i = sizes[0] - 2; i > 0; i--)
	{
		ch++;
		(*tab)[i][0] = *t;
		t++;
	}

}

void godunov1(vector<double>*** tab, list<point> lst, vector<int> sizes)
{
	interpolate(tab, lst, sizes);
	int it = 10;
	//resudual arrays
	double** xres;
	double** yres;

	try {
		xres = new double* [sizes[0]];
		yres = new double* [sizes[0]];
		for (int i = 0; i < sizes[0]; i++)
		{
			xres[i] = new double[sizes[1]];
			yres[i] = new double[sizes[1]];
			for (int j = 0; j < sizes[1]; j++)
			{
				xres[i][j] = 0;
				yres[i][j] = 0;
			}
		}
	}
	catch (bad_alloc xa) {
		return;
	}


	while (it > 0)
	{
		// dksi & deta
		vector<double> deta;
		vector<double> dksi;
		for (int i = 0; i < sizes[1] - 1; i++)
		{
			dksi.push_back((*tab)[i + 1][0][0] - (*tab)[i][0][0]);
		}
		for (int j = 0; j < sizes[0] - 1; j++)
		{
			deta.push_back((*tab)[0][j + 1][1] - (*tab)[0][j][1]);
		}
		double m = max(dksi);
		for (int i = 0; i < dksi.size(); i++)
		{
			deta[i] /= m;
		}
		m = max(deta);
		for (int j = 0; j < deta.size(); j++)
		{
			dksi[j] /= m;
		}
		deta.push_back(0);
		dksi.push_back(0);
		reverse(dksi.begin(), dksi.end());
		reverse(deta.begin(), deta.end());

		// Go table
		for (int i = 1; i < sizes[0] - 1; i++)
			for (int j = 1; j < sizes[1] - 1; j++)
			{
				//derivatives
				double xksi = ((*tab)[i + 1][j][0] - (*tab)[i - 1][j][0]) / (2 * dksi[i]);								// x/ksi deprivative 1
				double yksi = ((*tab)[i + 1][j][1] - (*tab)[i - 1][j][1]) / (2 * dksi[i]); 								// y/ksi deprivative 1
				double xeta = ((*tab)[i][j + 1][0] - (*tab)[i][j - 1][0]) / (2 * deta[j]);								// x/eta deprivative 1
				double yeta = ((*tab)[i][j + 1][1] - (*tab)[i][j - 1][1]) / (2 * deta[j]);						   	    // y/eta deprivative 1
				double xksi2 = ((*tab)[i - 1][j][0] - 2 * (*tab)[i][j][0] + (*tab)[i + 1][j][0]) / (dksi[i] * dksi[i]);   // x/ksi deprivative 2
				double yksi2 = ((*tab)[i - 1][j][1] - 2 * (*tab)[i][j][1] + (*tab)[i + 1][j][1]) / (dksi[i] * dksi[i]);    // y/ksi deprivative 2
				double xeta2 = ((*tab)[i][j - 1][0] - 2 * (*tab)[i][j][0] + (*tab)[i][j + 1][0]) / (deta[j] * deta[j]);   // x/eta deprivative 2
				double yeta2 = ((*tab)[i][j - 1][1] - 2 * (*tab)[i][j][1] + (*tab)[i][j + 1][1]) / (deta[j] * deta[j]);   // y/eta deprivative 2
				// Mixed
				double xksip = ((*tab)[i + 1][j - 1][0] - (*tab)[i - 1][j - 1][0]) / 2 / dksi[i]; // previos
				double xksin = ((*tab)[i + 1][j + 1][0] - (*tab)[i - 1][j + 1][0]) / 2 / dksi[i]; // next
				double xksieta = (xksin - xksip) / 2 / deta[j];
				double yksip = ((*tab)[i + 1][j - 1][1] - (*tab)[i - 1][j - 1][1]) / 2 / dksi[i]; // previos
				double yksin = ((*tab)[i + 1][j + 1][1] - (*tab)[i - 1][j + 1][1]) / 2 / dksi[i]; // next
				double yksieta = (yksin - yksip) / 2 / deta[j];

				// alfa, beta, gamma, teta
				double a = xeta * xeta + yeta * yeta;
				double b = xeta * xksi + yeta * yksi;
				double g = xksi * xksi + yksi * yksi;
				double t = 1 / (2 * (a + g) + abs(b));

				xres[i][j] = (a * xksi2 - 2 * b * xksieta + g * xeta2) * t;
				yres[i][j] = (a * yksi2 - 2 * b * yksieta + g * yeta2) * t;
			}
		// Add resuduals
		for (int i = 1; i < sizes[0] - 1; i++)
			for (int j = 1; j < sizes[1] - 1; j++)
			{
				(*tab)[i][j][0] += xres[i][j];
				(*tab)[i][j][1] += yres[i][j];
			}
		it--;
	}
}

void godunov2(vector<double>*** tab, list<point> lst, vector<int> sizes)
{
	//interpolate(tab, lst, sizes);
	godunov1(tab, lst, sizes);

	vector<vector<double>> xres(sizes[0], vector<double>(sizes[1], 0.0));
	vector<vector<double>> yres(sizes[0], vector<double>(sizes[1], 0.0));
	vector<vector<double>> lam(sizes[0] - 1, vector<double>(sizes[1] - 1, 0.0)); // lambda
	vector<vector<double>> om(sizes[0] - 1, vector<double>(sizes[1] - 1, 0.0)); // w
	vector<vector<double>> A(sizes[0] - 1, vector<double>(sizes[1] - 1, 0.0));
	vector<vector<double>> B(sizes[0] - 1, vector<double>(sizes[1] - 1, 0.0));
	vector<vector<double>> C(sizes[0] - 1, vector<double>(sizes[1] - 1, 0.0));
	vector<vector<double>> Adx(sizes[0], vector<double>(sizes[1], 0.0)); // differentials
	vector<vector<double>> Bdx(sizes[0], vector<double>(sizes[1], 0.0));
	vector<vector<double>> Cdx(sizes[0], vector<double>(sizes[1], 0.0));
	vector<vector<double>> Ady(sizes[0], vector<double>(sizes[1], 0.0));
	vector<vector<double>> Bdy(sizes[0], vector<double>(sizes[1], 0.0));
	vector<vector<double>> Cdy(sizes[0], vector<double>(sizes[1], 0.0));
	vector<vector<double>> teta(sizes[0], vector<double>(sizes[1], 0.0));

	bool it = true;
	int stop = 0;
	while (it)
	{
		// lambda star, w star
		double ls = 0;
		double ws = 0;
		// p & q
		vector<double> p;
		vector<double> q;
		// alfa, beta
		vector<double> al;
		vector<double> b;
		for (int i = 0; i < sizes[0] - 1; i++)
			for (int j = 0; j < sizes[1] - 1; j++)
			{
				//derivatives p = plus, m = minus
				double xksim = ((*tab)[i + 1][j][0] - (*tab)[i][j][0]);
				double yksim = ((*tab)[i + 1][j][1] - (*tab)[i][j][1]);
				double xetam = ((*tab)[i][j + 1][0] - (*tab)[i][j][0]);
				double yetam = ((*tab)[i][j + 1][1] - (*tab)[i][j][1]);
				double xksip = ((*tab)[i + 1][j + 1][0] - (*tab)[i][j + 1][0]);
				double yksip = ((*tab)[i + 1][j + 1][1] - (*tab)[i][j + 1][1]);
				double xetap = ((*tab)[i + 1][j + 1][0] - (*tab)[i + 1][j][0]);
				double yetap = ((*tab)[i + 1][j + 1][1] - (*tab)[i + 1][j][1]);
				double xksi = (xksip + xksim) / 2;
				double xeta = (xetap + xetam) / 2;
				double yeta = (yetap + yetam) / 2;
				double yksi = (yksip + yksim) / 2;
				// lambda, w
				double l = log(sqrt((xeta * xeta + yeta * yeta) / (xksi * xksi + yksi * yksi)));
				double w = acos((xeta * xksi + yeta * yksi) / (sqrt(xksi * xksi + yksi * yksi)
					* sqrt(xeta * xeta + yeta * yeta)));

				lam[i][j] = l;
				om[i][j] = w;


				ls += l;
				ws += w;
			}
		// Find q
		for (int i = 0; i < sizes[0] - 1; i++)
		{
			double qk = 0;
			for (int j = 0; j < sizes[1] - 1; j++)
			{
				qk -= lam[i][j];
			}
			qk /= sizes[1] - 1;
			q.push_back(qk);
		}
		// Find p
		for (int j = 0; j < sizes[1] - 1; j++)
		{
			double pk = 0;
			for (int i = 0; i < sizes[0] - 1; i++)
			{
				pk += lam[i][j];
			}
			pk /= sizes[0] - 1;
			p.push_back(pk);
		}
		// Find alfa
		for (int i = 0; i < sizes[0] - 1; i++)
		{
			double alk = 0;
			for (int j = 0; j < sizes[1] - 1; j++)
			{
				alk -= om[i][j];
			}
			alk /= sizes[1] - 1;
			al.push_back(alk);
		}
		// Find beta
		for (int j = 0; j < sizes[1] - 1; j++)
		{
			double bk = 0;
			for (int i = 0; i < sizes[0] - 1; i++)
			{
				bk += om[i][j];
			}
			bk /= sizes[0] - 1;
			b.push_back(bk);
		}

		// Compliting lambda*, w*, alfa, beta, p & q
		ls /= 2 * (sizes[0] - 1) * (sizes[1] - 1);
		ws /= 2 * (sizes[0] - 1) * (sizes[1] - 1);
		for (int k = 0; k < p.size(); k++) p[k] -= ls;
		for (int k = 0; k < q.size(); k++) q[k] += ls;
		for (int k = 0; k < b.size(); k++) b[k] -= ws;
		for (int k = 0; k < al.size(); k++) al[k] += ws;

		// A, B, C
		for (int i = 0; i < sizes[0] - 1; i++)
			for (int j = 0; j < sizes[1] - 1; j++)
			{

				A[i][j] = exp(p[j] - q[i]) / sin(b[j] - al[i]);
				B[i][j] = 1 / tan(b[j] - al[i]);
				C[i][j] = exp(q[i] - p[j]) / sin(b[j] - al[i]);
			}
		// A diferential (minus included)
		for (int i = 1; i < sizes[0] - 1; i++)
			for (int j = 1; j < sizes[1] - 1; j++)
			{
				Adx[i][j] = -0.5 * (A[i][j] + A[i][j - 1]) * ((*tab)[i + 1][j][0] - (*tab)[i][j][0])
					+ 0.5 * (A[i - 1][j] + A[i - 1][j - 1]) * ((*tab)[i][j][0] - (*tab)[i - 1][j][0]);
				Cdx[i][j] = -0.5 * (C[i][j] + C[i - 1][j]) * ((*tab)[i][j + 1][0] - (*tab)[i][j][0])
					+ 0.5 * (C[i][j - 1] + C[i - 1][j - 1]) * ((*tab)[i][j][0] - (*tab)[i][j - 1][0]);
				Bdx[i][j] = 0.5 * B[i][j] * ((*tab)[i + 1][j + 1][0] - (*tab)[i][j][0]) - 0.5 * B[i - 1][j] * ((*tab)[i - 1][j + 1][0] - (*tab)[i][j][0])
					- 0.5 * B[i][j - 1] * ((*tab)[i + 1][j - 1][0] - (*tab)[i][j][0]) + 0.5 * B[i - 1][j - 1] * ((*tab)[i - 1][j - 1][0] - (*tab)[i][j][0]);

				Ady[i][j] = -0.5 * (A[i][j] + A[i][j - 1]) * ((*tab)[i + 1][j][1] - (*tab)[i][j][1])
					+ 0.5 * (A[i - 1][j] + A[i - 1][j - 1]) * ((*tab)[i][j][1] - (*tab)[i - 1][j][1]);
				Cdy[i][j] = -0.5 * (C[i][j] + C[i - 1][j]) * ((*tab)[i][j + 1][1] - (*tab)[i][j][1])
					+ 0.5 * (C[i][j - 1] + C[i - 1][j - 1]) * ((*tab)[i][j][1] - (*tab)[i][j - 1][1]);
				Bdy[i][j] = 0.5 * B[i][j] * ((*tab)[i + 1][j + 1][1] - (*tab)[i][j][1]) - 0.5 * B[i - 1][j] * ((*tab)[i - 1][j + 1][1] - (*tab)[i][j][1])
					- 0.5 * B[i][j - 1] * ((*tab)[i + 1][j - 1][1] - (*tab)[i][j][1]) + 0.5 * B[i - 1][j - 1] * ((*tab)[i - 1][j - 1][1] - (*tab)[i][j][1]);
			}
		// teta
		for (int i = 0; i < sizes[0] - 1; i++)
			for (int j = 0; j < sizes[1] - 1; j++)
			{
				teta[i][j] = 1 / (2 * (A[i][j] + C[i][j]) + abs(B[i][j]));
			}
		// Count resudual
		for (int i = 1; i < sizes[0] - 1; i++)
			for (int j = 1; j < sizes[1] - 1; j++)
			{
				if (b[j] - al[i] <= 1e-6 || b[j] - al[i] >= 3.1415)
				{
					xres[i][j] = 0;
					yres[i][j] = 0;
				}
				else
				{
					double minteta = min(teta[i][j], teta[i - 1][j - 1], teta[i][j - 1], teta[i - 1][j]);
					xres[i][j] = minteta * (-Adx[i][j] - Cdx[i][j] - Bdx[i][j]);
					yres[i][j] = minteta * (-Ady[i][j] - Cdy[i][j] - Bdy[i][j]);
				}
			}
		// Stop condition
		stop++;
		if (stop > 30) it = false;
		/*		for (int i = 1; i < sizes[0] - 1; i++)
					for (int j = 1; j < sizes[1] - 1; j++)
					{
						double c1 = abs(((*tab)[i][j][0] - (*tab)[i - 1][j][0]));
						double c2 = abs(((*tab)[i + 1][j][0] - (*tab)[i][j][0]));
						double c3 = abs(((*tab)[i][j + 1][0] - (*tab)[i][j][0]));
						double c4 = abs(((*tab)[i][j][0] - (*tab)[i][j - 1][0]));
						double minside = min(c1, c2, c3, c4);
						if (abs(xres[i][j]) < 0.0001 * minside) it = false;
					}
				*/
				// Add resuduals
		for (int i = 1; i < sizes[0] - 1; i++)
			for (int j = 1; j < sizes[1] - 1; j++)
			{
				(*tab)[i][j][0] += xres[i][j];
				(*tab)[i][j][1] += yres[i][j];
			}


	}
}

void frame(vector<double>*** tab, list<point> u, vector<int> sizes)
{

	// putting key points into tha table
	// 1) define these points 
	list<point>::iterator p = u.begin();
	vector<point> keypoints;
	while (p != u.end())
	{
		if ((*p).log)
			keypoints.push_back(*p);
		p++;
	}
	int k = 0;
	// 2) put 'em into the table
	for (int j = 0; j < sizes[1]; j += sizes[1] - 1)
	{
		vector<double> v;
		v.push_back(keypoints[k].x);
		v.push_back(keypoints[k].y);
		(*tab)[0][j] = v;
		k++;
	}
	for (int j = (sizes[1] - 1); j >= 0; j -= (sizes[1] - 1))
	{
		vector<double> v;
		v.push_back(keypoints[k].x);
		v.push_back(keypoints[k].y);
		(*tab)[sizes[0] - 1][j] = v;
		k++;
	}

}

vector<string> getfiles()
{
	vector<string> res;
	system("dir ..\\input /B > files.txt");
	FILE* file = fopen("files.txt", "r");
	while (!feof(file))
	{
		char str[256];
		fscanf(file, "%s\n", &str);
		res.push_back(string(str));
	}
	fclose(file);
	system("del files.txt");
	return res;
}

string notxt(string s)
{
	string res = "";
	int i = 0;
	while (s[i] != '.')
	{
		res += s[i];
		i++;
	}
	return res;
}
int main()
{
	vector<string> files = getfiles();
	for (auto f = files.begin(); f != files.end(); f++)
	{
		list<point> u;
		u = read(*f);
		vector<int> sizes = getmatsize(u);
		// The table
		vector<double>** table;
		try {
			table = new vector<double>*[sizes[0]];
			for (int i = 0; i < sizes[0]; i++)
			{
				table[i] = new vector<double>[sizes[1]];
				for (int j = 0; j < sizes[1]; j++)
				{
					table[i][j].resize(2, 0);
				}
			}
		}
		catch (bad_alloc xa) {
			return -1;
		}

		frame(&table, u, sizes);
		distribute(&table, u, sizes);

		cout << *f << endl;

		// The first method
		interpolate(&table, u, sizes);
		string s = "../output/interpolation_" + notxt(*f) + ".dat";
		save(&table, sizes, s);
		cout << "\tinterpolation" << endl;

		// The second method
		godunov1(&table, u, sizes);
		s = "../output/godunov1_" + notxt(*f) + ".dat";
		save(&table, sizes, s);
		cout << "\tgodunov1" << endl;
			godunov2(&table, u, sizes);
			s = "../output/godunov2_" + notxt(*f) + ".dat";
			save(&table, sizes, s);
			cout << "\tgodunov2" << endl;

		delete[] table;
	}

	return 0;
}
