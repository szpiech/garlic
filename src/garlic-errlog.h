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
	void err(double val);
	void errn(double val);
	void err(int val);
	void errn(int val);
	void err(char val);
	void errn(char val);

	void err(string str, int val, bool nl = true);
	void err(string str, double val, bool nl = true);
	void err(string str, bool val, bool nl = true);
	void err(string str, string val, bool nl = true);
	void err(string str, char val, bool nl = true);

	void erra(string str, string *val, int size, bool nl = true);
	void erra(string str, int *val, int size, bool nl = true);
	void erra(string str, double *val, int size, bool nl = true);
	void erra(string str, char *val, int size, bool nl = true);

	void errv(string str, vector<int> &val, bool nl = true);
	void errv(string str, vector<double> &val, bool nl = true);

	void log(string str);
	void logn(string str);
	void log(double val);
	void logn(double val);
	void log(int val);
	void logn(int val);
	void log(char val);
	void logn(char val);

	void log(string str, int val, bool nl = true);
	void log(string str, double val, bool nl = true);
	void log(string str, bool val, bool nl = true);
	void log(string str, string val, bool nl = true);
	void log(string str, char val, bool nl = true);

	void loga(string str, string *val, int size, bool nl = true);
	void loga(string str, int *val, int size, bool nl = true);
	void loga(string str, double *val, int size, bool nl = true);
	void loga(string str, char *val, int size, bool nl = true);

	void logv(string str, vector<int> &val, bool nl = true);
	void logv(string str, vector<double> &val, bool nl = true);

private:
	ofstream *errstream;
	ofstream *logstream;

	void out(ostream *out, string str);
	void outn(ostream *out, string str);

	void out(ostream *out, double val);
	void outn(ostream *out, double val);
	void out(ostream *out, int val);
	void outn(ostream *out, int val);
	void out(ostream *out, char val);
	void outn(ostream *out, char val);

	void out(ostream *out, string str, int val, bool nl = true);
	void out(ostream *out, string str, double val, bool nl = true);
	void out(ostream *out, string str, bool val, bool nl = true);
	void out(ostream *out, string str, string val, bool nl = true);
	void out(ostream *out, string str, char val, bool nl = true);

	void outv(ostream *out, string str, vector<int> &val, bool nl = true);
	void outv(ostream *out, string str, vector<double> &val, bool nl = true);

	void outa(ostream *out, string str, string *val, int size, bool nl = true);
	void outa(ostream *out, string str, int *val, int size, bool nl = true);
	void outa(ostream *out, string str, double *val, int size, bool nl = true);
	void outa(ostream *out, string str, char *val, int size, bool nl = true);


};

extern errlog LOG;

#endif