#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include "unistd.h"
#include "polytope.h"
#include "glpk.h"

using namespace std;

#define MAX_LOOP_LIMIT 500

long total_generate;
long total_sample;

polyvest::polytope* ReadFile(string fileName)
{

	polyvest::polytope *p;
	
	string line;

	ifstream fileReader(fileName);
	if (fileReader.is_open())
	{
		vector<double> arr;
	
		// parse the file and store all numbers
		while (getline(fileReader, line))
		{
			if (line.find('#') != string::npos)
				line.erase(line.find('#'));
				
			istringstream iss(line);
			
			double num;
			while ((iss >> num))
				arr.push_back(num);
		}
		fileReader.close();
		
		if (arr.size() < 2)
		{
			cout << "Error: Lack file content! Please Check Input File." << endl;
			exit(-1);
		}
			
		unsigned m = arr[0];
		unsigned n = arr[1] - 2;
		
		if (arr.size() != m * (n + 2) + 2)
		{
			cout << "Error: Incorrect number of elements! Please Check Input File" << endl;
			exit(-1);
		}
			
		p = new polyvest::polytope(m, n);
		
		unsigned index = 2;
		
//		for (unsigned i = 0; i < n; i++)
//		{
//			p->var_flag[i] = arr[index];
//			index++;
//		}
		
		for (unsigned i = 0; i < m; i++) 
		{
			//P->vecOP[i] = (int)arr[index];	// 1: <=; 0: =;
			index++;		
		
			for (unsigned j = 0; j < n; j++)
			{
				p->matA(arr[index], i, j);
				index++;
			}
			
			p->vecb(arr[index], i);
			index++;
		}
		
		return p;
		
	}
	else
	{
		cout << "Error: " << fileName << " cannot open!" << endl;
		exit(-1);
	}
	
}

void print_help()
{
	cout << "Usage: ./ApproxLatCount <infile> [outfile]\n";
//	cout << "Options:\n";
//	cout << "  -h or -help\t\tPrint help.\n";
//	cout << "  -a or -approx\t\tEmploy approximate method (by default).\n";
//	cout << "    -seed=..\t\t  Set random seed.\n";
//	cout << "    -epsilon=..\t\t  Set accuracy parameter (by default epsilon = 0.1).\n";
//	cout << "    -delta=..\t\t  Set probability parameter (by default delta = 0.05).\n";
//	cout << "    -maxs=..\t\t  Set the max number of samples (by default maxs = 100000).\n";
//	cout << "    -vc or -vinci\t  Employ volume computation tool Vinci.\n";
//	cout << "    -pv or -polyvest\t  Employ volume estimation tool PolyVest.\n";
//	cout << "  -e or -exact\t\tEmploy exact method.\n";
}

void SampleGenWR(polyvest::polytope *P, const int maxs, vector<arma::vec> &samples)
{
	int n = P->n;
		
	// copy from P then enlarge it
	polyvest::polytope *bigP = P->Clone();
	bigP->Enlarge(0.5);
	
	// apply ellipsoid
	polyvest::polytope *bigQ = bigP->Clone();
	bigQ->AffineTrans();
	bigQ->PrepForWalk();
	if (samples.size() > 0)
		bigQ->x = bigQ->GetTransPoint(samples[rand() % samples.size()]);
	else if (!bigQ->isInside(bigQ->x))
		bigQ->x.zeros();

	long generate_count = 0;
	long sample_count = 0;
	arma::vec point(n);

	while ((int)samples.size() < maxs)
	{		
		do 
		{
			// at least n steps
			for (int i = 0; i < n; i++) 
				bigQ->Walk();
			point = bigQ->GetInvPoint(bigQ->x);

			for (int i = 0; i < n; i++)
				point(i) = round(point(i));
					
			generate_count++;
					
		} while (!P->isInside(point));

		samples.push_back(point);
		
		sample_count++;
	}

	total_generate += generate_count;
	total_sample += sample_count;
	cout << "Samples / Generated: " << sample_count << " / " << generate_count << endl;

	delete bigQ;
	delete bigP;

}

void SampleGen(polyvest::polytope *P, const int maxs, vector<arma::vec> &samples)
{
	int n = P->n;

	// copy from P then enlarge it
	polyvest::polytope *bigP = P->Clone();
	bigP->Enlarge(0.5);			
	bigP->PrepForWalk();

	if (samples.size() > 0)
		bigP->x = samples[rand() % samples.size()];
	else if (!bigP->isInside(bigP->x))
	{
		polyvest::polytope *bigQ = bigP->Clone();
		bigQ->AffineTrans();
		bigQ->x.zeros();
		bigP->x = bigQ->GetInvPoint(bigQ->x);
		delete bigQ;
	}

	long generate_count = 0;
	long sample_count = 0;
	arma::vec point(n);

	while ((int)samples.size() < maxs)
	{
		do 
		{
			// at least n steps
			for (int i = 0; i < n; i++)
				bigP->Walk();
			point = bigP->x;

			for (int i = 0; i < n; i++)
				point(i) = round(point(i));
					
			generate_count++;
				
		} while (!P->isInside(point));

		samples.push_back(point);

		sample_count++;
	}
	
	total_generate += generate_count;
	total_sample += sample_count;
	cout << "Samples / Generated: " << sample_count << " / " << generate_count << endl;
	
	delete bigP;
}

polyvest::polytope* GenRect(polyvest::polytope* P, bool callILP)
{
	int n = P->n;

	polyvest::polytope *Rect = new polyvest::polytope(2 * n, n);

	double count = 1;
	
	if (callILP)
	{
		P->GetBounds(true);
		for (int i = 0; i < n; i++)
		{
			count *= round(P->xmax[i]) - round(P->xmin[i]) + 1;
		
			Rect->matA(1, i, i);
			Rect->vecb(round(P->xmax[i]), i);
			Rect->matA(-1, i + n, i);
			Rect->vecb(-round(P->xmin[i]), i + n);
		}
	} else
	{
		P->GetBounds(false);
		for (int i = 0; i < n; i++)
		{
			count *= floor(P->xmax[i]) - ceil(P->xmin[i]) + 1;
		
			Rect->matA(1, i, i);
			Rect->vecb(floor(P->xmax[i]), i);
			Rect->matA(-1, i + n, i);
			Rect->vecb(-ceil(P->xmin[i]), i + n);
		}
	}
	
	Rect->latcount = count;
	
	return Rect;
}

bool FindB(arma::vec &a, double &b, vector<arma::vec> &samples, const int lb, const int ub)
{
	double originb = b;

	const int mid = (lb + ub) / 2 - 1;
	int lidx = 0;
	int bidx = 0;
	vector<double> dist;
	vector<arma::vec> newsamples;
	
	for (long unsigned i = 0; i < samples.size(); i++)
	{
		arma::mat matval = a.t() * samples[i];
		double val = matval(0, 0);
		dist.push_back(val);
		if (val < b) bidx++;
	}
	
	if (bidx >= (int)dist.size()) return true;
	
	sort(dist.begin(), dist.end());
	
	if (dist[mid] > dist[bidx])
	{
		// find the index which is closest to mid
		lidx = mid;
		while (dist[lidx] == dist[lidx - 1] && lidx > 0) lidx--;
		lidx = mid * 2 - lidx;
		while (dist[lidx] == dist[lidx - 1] && lidx > 0) lidx--;
		lidx--;
	
		if (lidx < lb || lidx > ub) 
			return false;
		else 
			b = dist[lidx];
	}
		
	for (long unsigned i = 0; i < samples.size(); i++)
	{
		arma::mat matval = a.t() * samples[i];
		double val = matval(0, 0);
		if (val <= b)
			newsamples.push_back(samples[i]);
	}
	
	samples.resize(0);
	samples = newsamples;
	
	assert(b >= originb);
	
	return true;
}

// P \cap ax<=b
void AddConstraint(polyvest::polytope* &P, arma::vec a, double b)
{
	int m = P->m;
	int n = P->n;
	
	// from bottom find the ith row of A which is same as vec a
	for (int i = m - 1; i >= 0; i--)
	{
		bool same = true;
		for (int j = 0; j < n; j++)
			if (a(j) != P->matA(i, j))
			{
				same = false;
				break;
			}
		if (same) 
		{
			if (b < P->vecb(i))
				P->vecb(b, i);
			return;
		}
	}
	
	// add new constraint newa * x <= newb
	polyvest::polytope *tP = new polyvest::polytope(m + 1, n);
	
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			tP->matA(P->matA(i, j), i, j);
		tP->vecb(P->vecb(i), i);
	}
	for (int j = 0; j < n; j++)
		tP->matA(a(j), m, j);
	tP->vecb(b, m);
	
	delete P;
	P = tP;
}

double GetMean(vector<double> &x)
{
	double sum = accumulate(x.begin(), x.end(), 0.0);
	return sum / x.size();
}

double GetVariance(vector<double> &x)
{
	double mean = GetMean(x);
	double variance = 0.0;
	for (long unsigned i = 0; i < x.size(); i++)
		variance += pow(x[i] - mean, 2);
	variance /= x.size();
	return variance;
}

// IN:	ost, ratios, vars, RectSize, delta
// OUT:	approx C(P), err
double FinalApprox(int msg_lvl, ostream &ost, vector<double> &ratios, vector<double> &vars, double &err, const double RectSize, const double delta)
{

	double ratio = 1.0;
	double var = 1.0;
	double var_r = 1.0;

	if (msg_lvl == 2) 
		ost << "Ratios:\t1.0\n";

	for (long unsigned i = 0; i < ratios.size(); i++)
	{
		ratio *= ratios[i];
		
		double rsqr = pow(ratios[i], 2);
		var *= vars[i] + rsqr;
		var_r *= rsqr;
		
		if (msg_lvl == 2) 
			ost << " x\tE " << ratios[i] << ", Var " << vars[i] << endl;
	}
	var -= var_r;
	err = sqrt(var / delta) / ratio;

	if (msg_lvl == 2)
	{
		ost << " =\tE " << ratio << ", V " << var << ", err " << err << endl;
		ost << "RectSize: " << RectSize << endl;
		ost << "Est C(P) = Ratio * RectSize = " << ratio * RectSize << endl;
		ost << "Total Samples / Generated: " << total_sample << " / " << total_generate << endl;
	} else if (msg_lvl == 1)
	{
		ost << ratios.size() << " " << ratio << " " << var << " " << err << " " << RectSize << " " 
			<< ratio * RectSize << " " << total_sample << " " << total_generate;
	}
	
	return ratio * RectSize;
}

// IN: 	p, samples, rlist, ratio, gcount, loops
// OUT:	ratio, var
void SplitSamples(polyvest::polytope *p, vector<arma::vec> &samples, vector<double> &rlist, double &ratio, double &var, const int gcount, const int loops)
{

	int count = 0;
	int ssize = samples.size();	
	int gsize = ssize / gcount;
	vector<arma::vec> newsamples;
			
	for (int i = 0; i < gcount; i++)
	{
		int tmpc = 0;
		for (int j = 0; j < gsize; j++)
		{
			if (p->isInside(samples[i * gsize + j]))
			{
				newsamples.push_back(samples[i * gsize + j]);
				count++, tmpc++;
			}
		}
		rlist.push_back((double)tmpc / gsize);
	}
	samples.resize(0);
	samples = newsamples;

	ratio = (double)(ratio * ssize * loops + count) / (ssize * loops + ssize);
	var= GetVariance(rlist) / rlist.size();

}


int main(int argc, char **argv) 
{

	bool enableRounding = true;

	double 	epsilon = 0.2;
	double 	delta = 0.1;
	double  upper_coef = 0.6;
	double  lower_coef = 0.4;
	long unsigned maxs = round(2 * 1.0 / (delta * epsilon * epsilon));
	int 	gcount = 10;
	
	auto timepoint = chrono::system_clock::now();
    auto fine = chrono::time_point_cast<chrono::milliseconds>(timepoint);
	srand((int)fine.time_since_epoch().count());

	if (argc < 2)
	{
		print_help();
		exit(-1);
	}

	string inFileName = "";
	string outFileName = "";

	for (int i = 1; i < argc; i++)
	{
		string arg_str = argv[i];
		int len = arg_str.length();

		if (arg_str[0] != '-' || len < 2)
		{
			if (inFileName == "")
				inFileName = arg_str;
			else if (outFileName == "")
				outFileName = arg_str;
			else
			{
				cout << "Error: Multiple input or output files.\n";
				exit(-1);
			}
		}
		
		if (arg_str == "-h" || arg_str == "-help")
		{
			print_help();
			exit(1);
		}
	}
	
	
	auto t1 = chrono::high_resolution_clock::now();
	
	polyvest::polytope *P;
	vector<polyvest::polytope *> PL;
	vector<arma::vec> samples;
	vector<arma::vec> samples_backup;
	vector<double> ratios;
	vector<double> vars;
	vector< vector<double> > rll;
	double latcount;
	
	P = ReadFile(inFileName);
//	P->Print();
	cout << "Read file \'" << inFileName << "\' finished.\n";
	
	int m = P->m;
	int n = P->n;
	
	bool callILP = true;
	if (n >= 20) callILP = false;
	PL.push_back(GenRect(P, callILP));
	
	int Pidx = 0;
	int HPidx = 0;
	total_generate = 0;
	total_sample = 0;
	
	cout << "\nROUNDS: 1" << endl;
	
	while (Pidx < MAX_LOOP_LIMIT && HPidx < m)
	{
		polyvest::polytope *Pi = PL[Pidx];
		polyvest::polytope *newP;
	
		if (enableRounding)
			SampleGenWR(Pi, maxs, samples);
		else
			SampleGen(Pi, maxs, samples);
			
		samples_backup.resize(0);
		samples_backup = samples;	

		int count = 0;
		for (long unsigned i = 0; i < samples.size(); i++)
			if (P->isInside(samples[i])) 
				count++;
		
		if (count > lower_coef * maxs)
			newP = P, HPidx = m;
		else
			newP = Pi->Clone();

		//newP->Print();

		while (samples.size() > upper_coef * maxs && HPidx < m)
		{			
			arma::vec a = P->A.row(HPidx).t();			
			arma::vec newa = a;
			double b = P->vecb(HPidx);
			double newb = b;
			bool accept = FindB(newa, newb, samples, lower_coef * maxs, upper_coef * maxs);
			//cout << "TEST A " << HPidx << " " << newb << endl;
			
			if (accept && newb == b)
				HPidx++;

			while (!accept)
			{
				newa = a, newb = b;
				for (int i = 0; i < n; i++) 
					newa(i) *= 1 + ((double)(rand() % 10000) / 10000 - 0.5) * 0.01;
				accept = FindB(newa, newb, samples, lower_coef * maxs, upper_coef * maxs);
				//cout << "TEST C " << HPidx << " " << newb << endl;
			}

			AddConstraint(newP, newa, newb);
			
		}
		
		double ratio, var;
		vector<double> rl;
		SplitSamples(newP, samples_backup, rl, ratio, var, gcount, 0);
		rll.push_back(rl);
		ratios.push_back(ratio);
		vars.push_back(var);

		PL.push_back(newP), Pidx++;
	}
	
	double err = 0;
	latcount = FinalApprox(2, cout, ratios, vars, err, PL[0]->LatticeCount(), delta);
	
	int loops = 2;
	while (err > epsilon)
	{
		cout << "\nROUNDS: " << loops << endl;

		samples.resize(0);

		for (long unsigned Pidx = 0; Pidx < PL.size() - 1; Pidx++)
		{
		
			polyvest::polytope *Pi = PL[Pidx];
			polyvest::polytope *Pj = PL[Pidx + 1];
	
			if (enableRounding)
				SampleGenWR(Pi, maxs, samples);
			else
				SampleGen(Pi, maxs, samples);

			assert(maxs == samples.size());
			
			SplitSamples(Pj, samples, rll[Pidx], ratios[Pidx], vars[Pidx], gcount, loops);

		}
		
		latcount = FinalApprox(2, cout, ratios, vars, err, PL[0]->LatticeCount(), delta);
		
		loops++;
	
	}

	auto t2 = std::chrono::high_resolution_clock::now();

	double timecost = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / 1000000000;
	cout << "Runime: " << timecost << endl;

	if (outFileName != "")
	{
		ofstream ofile(outFileName, ios::app);
		latcount = FinalApprox(1, ofile, ratios, vars, err, PL[0]->LatticeCount(), delta);
		ofile << " " << timecost << endl;
		ofile.close();
	}
	
	//ofstream ofile("tmp.out");
	//ofile << latcount << endl;
	//ofile.close();
/*
	double ExactP = P->ExactCount();
	cout << "Exa C(P):\t" << ExactP << " " << latcount / ExactP << endl;	


	if (latcount / ExactP < 0.97)
	{
		double CountPi = PL[0]->ExactCount();
		for (long unsigned i = 0; i < PL.size() - 1; i++)
		{
			double CountPj = PL[i + 1]->ExactCount();
			cout << "CountP: " << CountPi << " " << CountPj << " " << CountPj / CountPi << endl;
			CountPi = CountPj;
		}
	}
*/
	for (long unsigned i = 0; i < PL.size() - 1; i++)
		delete PL[i];
	delete P;

	return 1;

}





