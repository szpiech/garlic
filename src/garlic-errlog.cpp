#include "garlic-errlog.h"

errlog::errlog()
{
	logfile = "_none";
	errfile = "_none";
	errstream = NULL;
	logstream = NULL;
}

errlog::errlog(string file)
{
	this->init(file);
}

errlog::~errlog()
{
	if(errstream) errstream->close();
	if(logstream) logstream->close();
	delete errstream;
	delete logstream;
}

void errlog::init(string file)
{
	logfile = file;
	logfile += ".log";

	logstream = new ofstream;
	logstream->open(logfile.c_str());

	if (logstream->fail())
	{
		cerr << "ERROR: Could not open " << logfile << " for logging.\n";
		throw 0;
	}

	errfile = file;
	errfile += ".error";

	errstream = new ofstream;
	errstream->open(errfile.c_str());

	if (errstream->fail())
	{
		cerr << "ERROR: Could not open " << errfile << " for logging.\n";
		throw 0;
	}

	return;
}

void errlog::errn(string str)
{
	this->outn(errstream, str);
	return;
}

void errlog::err(string str)
{
	this->out(errstream, str);
	return;
}

void errlog::err(string str, int val, bool nl)
{
	this->out(errstream, str, val, nl);
	return;
}

void errlog::errv(string str, vector<int> &val, bool nl)
{
	this->outv(errstream, str, val, nl);
	return;
}

void errlog::err(string str, double val, bool nl)
{
	this->out(errstream, str, val, nl);
	return;
}

void errlog::errv(string str, vector<double> &val, bool nl)
{
	this->outv(errstream, str, val, nl);
	return;
}

void errlog::err(string str, bool val, bool nl)
{
	this->out(errstream, str, val, nl);
	return;
}

void errlog::err(string str, string val, bool nl)
{
	this->out(errstream, str, val, nl);
	return;
}

void errlog::err(string str, char val, bool nl)
{
	this->out(errstream, str, val, nl);
	return;
}

void errlog::logn(string str)
{
	this->outn(logstream, str);
	return;
}

void errlog::log(string str)
{
	this->out(logstream, str);
	return;
}

void errlog::log(string str, int val, bool nl)
{
	this->out(logstream, str, val, nl);
	return;
}

void errlog::logv(string str, vector<int> &val, bool nl)
{
	this->outv(logstream, str, val, nl);
	return;
}

void errlog::log(string str, double val, bool nl)
{
	this->out(logstream, str, val, nl);
	return;
}

void errlog::logv(string str, vector<double> &val, bool nl)
{
	this->outv(logstream, str, val, nl);
	return;
}

void errlog::log(string str, bool val, bool nl)
{
	this->out(logstream, str, val, nl);
	return;
}

void errlog::log(string str, string val, bool nl)
{
	this->out(logstream, str, val, nl);
	return;
}

void errlog::log(string str, char val, bool nl)
{
	this->out(logstream, str, val, nl);
	return;
}

void errlog::outn(ofstream *out, string str)
{
	if (out)
	{
		*(out) << str;
		out->flush();
	}
	return;
}

void errlog::out(ofstream *out, string str)
{
	if (out)
	{
		*(out) << str << "\n";
		out->flush();
	}
	return;
}

void errlog::out(ofstream *out, string str, int val, bool nl)
{
	if (out)
	{
		*(out) << str << " " << val;
		if (nl) *(out) << "\n";
		out->flush();
	}
	return;
}

void errlog::outv(ofstream *out, string str, vector<int> &val, bool nl)
{
	if (out)
	{
		*(out) << str;
		for (unsigned int i = 0; i < val.size(); i++) *(out) << " " << val[i];
		if (nl) *(out) << "\n";
		out->flush();
	}
	return;
}

void errlog::out(ofstream *out, string str, double val, bool nl)
{
	if (out)
	{
		*(out) << str << " " << val;
		if (nl) *(out) << "\n";
		out->flush();
	}
	return;
}

void errlog::outv(ofstream *out, string str, vector<double> &val, bool nl)
{
	if (out)
	{
		*(out) << str;
		for (unsigned int i = 0; i < val.size(); i++) *(out) << " " << val[i];
		if (nl) *(out) << "\n";
		out->flush();
	}
	return;
}

void errlog::out(ofstream *out, string str, bool val, bool nl)
{
	if (out)
	{
		string b = val ? "TRUE" : "FALSE";
		*(out) << str << " " << b;
		if (nl) *(out) << "\n";
		out->flush();
	}
	return;
}

void errlog::out(ofstream *out, string str, string val, bool nl)
{
	if (out)
	{
		*(out) << str << " " << val;
		if (nl) *(out) << "\n";
		out->flush();
	}
	return;
}

void errlog::out(ofstream *out, string str, char val, bool nl)
{
	if (out)
	{
		*(out) << str << " " << val;
		if (nl) *(out) << "\n";
		out->flush();
	}
	return;
}

errlog LOG;
