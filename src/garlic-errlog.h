#ifndef __GARLIC_ERRLOG_H__
#define __GARLIC_ERRLOG_H__

#include <iostream>
#include <fstream>
#include <sstream>
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

	//init() only records the output basename.  Log output accumulates in memory
	//until commit(), so an invocation that fails argument validation does not
	//create or truncate <out>.log / <out>.error.  <out>.error is created only
	//if something is actually written to it.
	void init(string file);
	void commit();

	//Verbosity.  --quiet silences progress and informational chatter (errors
	//still go to stderr); --verbose forces the progress bar on even when
	//stderr is not a terminal.
	void setVerbosity(bool quiet, bool verbose);
	bool isQuiet();
	bool isVerbose();

	void err(string str);
	void errn(string str);
	void err(double val);
	void errn(double val);
	void err(int val);
	void errn(int val);
	void err(char val);
	void errn(char val);

	void err(string str, int val, bool nl = true);
	//A physical position is a pos_t (int64), which matches neither the int
	//nor the double overload; without this, printing one is ambiguous, and
	//routing it through double would render large values in scientific
	//notation.
	void err(string str, long long val, bool nl = true);
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
	void log(string str, long long val, bool nl = true);
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
	ostream *errstream;
	ostream *logstream;
	ostringstream logbuf;
	ostringstream errbuf;
	ofstream *errfilestream;
	ofstream *logfilestream;
	bool committed;
	bool quiet;
	bool verbose;
	ostream *errOut();

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

//Reports the exception currently being handled.  Call it from inside a catch
//block; it re-throws into its own typed handlers to recover the type.
//
//Every `catch (...)` in garlic used to discard the exception AND any
//diagnostic.  That was harmless for garlic's own failures, which log a cause
//before throwing an int, but silent for anything the standard library raises:
//there are ~130 vector::at() calls inside these try blocks, and a bad_alloc on
//a large callset, so a run could exit 2 having printed nothing at all.
void logCurrentException(const string &stage);

extern errlog LOG;

#endif