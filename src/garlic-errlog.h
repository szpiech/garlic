#ifndef __GARLIC_ERRLOG_H__
#define __GARLIC_ERRLOG_H__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

class errlog {
public:
	errlog();
	~errlog();
	errlog(string file);

	string logfile;
	string errfile;

	void init(string file);

	void err(string str);
	void errn(string str);

	void err(string str, int val, bool nl = true);
	void err(string str, double val, bool nl = true);
	void err(string str, bool val, bool nl = true);
	void err(string str, string val, bool nl = true);
	void err(string str, char val, bool nl = true);

	void errv(string str, vector<int> &val, bool nl = true);
	void errv(string str, vector<double> &val, bool nl = true);

	void log(string str);
	void logn(string str);

	void log(string str, int val, bool nl = true);
	void log(string str, double val, bool nl = true);
	void log(string str, bool val, bool nl = true);
	void log(string str, string val, bool nl = true);
	void log(string str, char val, bool nl = true);

	void logv(string str, vector<int> &val, bool nl = true);
	void logv(string str, vector<double> &val, bool nl = true);

private:
	ofstream *errstream;
	ofstream *logstream;

	void out(ofstream *out, string str);
	void outn(ofstream *out, string str);
	
	void out(ofstream *out, string str, int val, bool nl = true);
	void out(ofstream *out, string str, double val, bool nl = true);
	void out(ofstream *out, string str, bool val, bool nl = true);
	void out(ofstream *out, string str, string val, bool nl = true);
	void out(ofstream *out, string str, char val, bool nl = true);

	void outv(ofstream *out, string str, vector<int> &val, bool nl = true);
	void outv(ofstream *out, string str, vector<double> &val, bool nl = true);
};

extern errlog LOG;

#endif