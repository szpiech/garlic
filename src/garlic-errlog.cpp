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
	if (errstream) errstream->close();
	if (logstream) logstream->close();
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
	this->outn(&cerr, str);
	this->outn(errstream, str);
	return;
}

void errlog::err(string str)
{
	this->out(&cerr, str);
	this->out(errstream, str);
	return;
}

void errlog::err(double val)
{
	this->out(&cerr, val);
	this->out(errstream, val);
	return;
}

void errlog::errn(double val)
{
	this->out(&cerr, val);
	this->out(errstream, val);
	return;
}

void errlog::err(int val)
{
	this->out(&cerr, val);
	this->out(errstream, val);
	return;
}

void errlog::errn(int val)
{
	this->out(&cerr, val);
	this->out(errstream, val);
	return;
}

void errlog::err(char val)
{
	this->out(&cerr, val);
	this->out(errstream, val);
	return;
}

void errlog::errn(char val)
{
	this->out(&cerr, val);
	this->out(errstream, val);
	return;
}

void errlog::err(string str, int val, bool nl)
{
	this->out(&cerr, str, val, nl);
	this->out(errstream, str, val, nl);
	return;
}

void errlog::errv(string str, vector<int> &val, bool nl)
{
	this->outv(&cerr, str, val, nl);
	this->outv(errstream, str, val, nl);
	return;
}

void errlog::err(string str, double val, bool nl)
{
	this->out(&cerr, str, val, nl);
	this->out(errstream, str, val, nl);
	return;
}

void errlog::errv(string str, vector<double> &val, bool nl)
{
	this->outv(&cerr, str, val, nl);
	this->outv(errstream, str, val, nl);
	return;
}

void errlog::err(string str, bool val, bool nl)
{
	this->out(&cerr, str, val, nl);
	this->out(errstream, str, val, nl);
	return;
}

void errlog::err(string str, string val, bool nl)
{
	this->out(&cerr, str, val, nl);
	this->out(errstream, str, val, nl);
	return;
}

void errlog::err(string str, char val, bool nl)
{
	this->out(&cerr, str, val, nl);
	this->out(errstream, str, val, nl);
	return;
}

void errlog::erra(string str, string *val, int size, bool nl)
{
	this->outa(&cerr, str, val, size, nl);
	this->outa(errstream, str, val, size, nl);
	return;
}

void errlog::erra(string str, int *val, int size, bool nl)
{
	this->outa(&cerr, str, val, size, nl);
	this->outa(errstream, str, val, size, nl);
	return;
}

void errlog::erra(string str, double *val, int size, bool nl)
{
	this->outa(&cerr, str, val, size, nl);
	this->outa(errstream, str, val, size, nl);
	return;
}

void errlog::erra(string str, char *val, int size, bool nl)
{
	this->outa(&cerr, str, val, size, nl);
	this->outa(errstream, str, val, size, nl);
	return;
}

void errlog::logn(string str)
{
	this->outn(&cout, str);
	this->outn(logstream, str);
	return;
}

void errlog::log(string str)
{
	this->out(&cout, str);
	this->out(logstream, str);
	return;
}

void errlog::log(double val)
{
	this->out(&cout, val);
	this->out(logstream, val);
	return;
}

void errlog::logn(double val)
{
	this->out(&cout, val);
	this->out(logstream, val);
	return;
}

void errlog::log(int val)
{
	this->out(&cout, val);
	this->out(logstream, val);
	return;
}

void errlog::logn(int val)
{
	this->out(&cout, val);
	this->out(logstream, val);
	return;
}

void errlog::log(char val)
{
	this->out(&cout, val);
	this->out(logstream, val);
	return;
}

void errlog::logn(char val)
{
	this->out(&cout, val);
	this->out(logstream, val);
	return;
}

void errlog::log(string str, int val, bool nl)
{
	this->out(&cout, str, val, nl);
	this->out(logstream, str, val, nl);
	return;
}

void errlog::logv(string str, vector<int> &val, bool nl)
{
	this->outv(&cout, str, val, nl);
	this->outv(logstream, str, val, nl);
	return;
}

void errlog::log(string str, double val, bool nl)
{
	this->out(&cout, str, val, nl);
	this->out(logstream, str, val, nl);
	return;
}

void errlog::logv(string str, vector<double> &val, bool nl)
{
	this->outv(&cout, str, val, nl);
	this->outv(logstream, str, val, nl);
	return;
}

void errlog::log(string str, bool val, bool nl)
{
	this->out(&cout, str, val, nl);
	this->out(logstream, str, val, nl);
	return;
}

void errlog::log(string str, string val, bool nl)
{
	this->out(&cout, str, val, nl);
	this->out(logstream, str, val, nl);
	return;
}

void errlog::log(string str, char val, bool nl)
{
	this->out(&cout, str, val, nl);
	this->out(logstream, str, val, nl);
	return;
}

void errlog::loga(string str, string *val, int size, bool nl)
{
	this->outa(&cout, str, val, size, nl);
	this->outa(logstream, str, val, size, nl);
	return;
}

void errlog::loga(string str, int *val, int size, bool nl)
{
	this->outa(&cout, str, val, size, nl);
	this->outa(logstream, str, val, size, nl);
	return;
}

void errlog::loga(string str, double *val, int size, bool nl)
{
	this->outa(&cout, str, val, size, nl);
	this->outa(logstream, str, val, size, nl);
	return;
}

void errlog::loga(string str, char *val, int size, bool nl)
{
	this->outa(&cout, str, val, size, nl);
	this->outa(logstream, str, val, size, nl);
	return;
}

void errlog::outn(ostream *out, string str)
{
	if (out)
	{
		*(out) << str;
		out->flush();
	}
	return;
}

void errlog::out(ostream *out, string str)
{
	if (out)
	{
		*(out) << str << endl;
		out->flush();
	}
	return;
}

void errlog::out(ostream *out, string str, int val, bool nl)
{
	if (out)
	{
		*(out) << str << " " << val;
		if (nl) *(out) << endl;
		out->flush();
	}
	return;
}

void errlog::outv(ostream *out, string str, vector<int> &val, bool nl)
{
	if (out)
	{
		*(out) << str;
		for (unsigned int i = 0; i < val.size(); i++) *(out) << " " << val[i];
		if (nl) *(out) << endl;
		out->flush();
	}
	return;
}

void errlog::out(ostream *out, string str, double val, bool nl)
{
	if (out)
	{
		*(out) << str << " " << val;
		if (nl) *(out) << endl;
		out->flush();
	}
	return;
}

void errlog::outv(ostream *out, string str, vector<double> &val, bool nl)
{
	if (out)
	{
		*(out) << str;
		for (unsigned int i = 0; i < val.size(); i++) *(out) << " " << val[i];
		if (nl) *(out) << endl;
		out->flush();
	}
	return;
}

void errlog::out(ostream *out, string str, bool val, bool nl)
{
	if (out)
	{
		string b = val ? "TRUE" : "FALSE";
		*(out) << str << " " << b;
		if (nl) *(out) << endl;
		out->flush();
	}
	return;
}

void errlog::out(ostream *out, string str, string val, bool nl)
{
	if (out)
	{
		*(out) << str << " " << val;
		if (nl) *(out) << endl;
		out->flush();
	}
	return;
}

void errlog::out(ostream *out, string str, char val, bool nl)
{
	if (out)
	{
		*(out) << str << " " << val;
		if (nl) *(out) << endl;
		out->flush();
	}
	return;
}

void errlog::outa(ostream *out, string str, string *val, int size, bool nl)
{
	if (out && size > 0)
	{
		*(out) << str;
		for (int i = 0; i < size; i++) *(out) << " " << val[i];
		if (nl) *(out) << endl;
		out->flush();
	}
	return;
}

void errlog::outa(ostream *out, string str, int *val, int size, bool nl)
{
	if (out && size > 0)
	{
		*(out) << str;
		for (int i = 0; i < size; i++) *(out) << " " << val[i];
		if (nl) *(out) << endl;
		out->flush();
	}
	return;
}

void errlog::outa(ostream *out, string str, double *val, int size, bool nl)
{
	if (out && size > 0)
	{
		*(out) << str;
		for (int i = 0; i < size; i++) *(out) << " " << val[i];
		if (nl) *(out) << endl;
		out->flush();
	}
	return;
}

void errlog::outa(ostream *out, string str, char *val, int size, bool nl)
{
	if (out && size > 0)
	{
		*(out) << str;
		for (int i = 0; i < size; i++) *(out) << " " << val[i];
		if (nl) *(out) << endl;
		out->flush();
	}
	return;
}

void errlog::out(ostream *out, double val)
{
	if (out)
	{
		*(out) << val << endl;
		out->flush();
	}
	return;
}

void errlog::outn(ostream *out, double val)
{
	if (out)
	{
		*(out) << val;
		out->flush();
	}
	return;
}

void errlog::out(ostream *out, int val)
{
	if (out)
	{
		*(out) << val << endl;
		out->flush();
	}
	return;
}

void errlog::outn(ostream *out, int val)
{
	if (out)
	{
		*(out) << val;
		out->flush();
	}
	return;
}

void errlog::out(ostream *out, char val)
{
	if (out)
	{
		*(out) << val << endl;
		out->flush();
	}
	return;
}

void errlog::outn(ostream *out, char val)
{
	if (out)
	{
		*(out) << val;
		out->flush();
	}
	return;
}

errlog LOG;
