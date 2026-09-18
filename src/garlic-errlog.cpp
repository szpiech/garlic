#include "garlic-errlog.h"
#include <exception>
#include <stdexcept>

errlog::errlog()
{
	logfile = "_none";
	errfile = "_none";
	errstream = &errbuf;
	logstream = &logbuf;
	errfilestream = NULL;
	logfilestream = NULL;
	committed = false;
	quiet = false;
	verbose = false;
}

void errlog::setVerbosity(bool q, bool v) { quiet = q; verbose = v; }
bool errlog::isQuiet() { return quiet; }
bool errlog::isVerbose() { return verbose; }

errlog::errlog(string file)
{
	this->init(file);
}

errlog::~errlog()
{
	if (logfilestream) { logfilestream->close(); delete logfilestream; }
	if (errfilestream) { errfilestream->close(); delete errfilestream; }
}

void errlog::init(string file)
{
	logfile = file;
	logfile += ".log";
	errfile = file;
	errfile += ".error";
	//Nothing is opened here on purpose; see commit().
	return;
}

void errlog::commit()
{
	if (committed) return;

	logfilestream = new ofstream;
	//ios::binary: see writeROHData in garlic-roh.cpp.
	logfilestream->open(logfile.c_str(), ios::binary);
	if (logfilestream->fail())
	{
		cerr << "ERROR: Could not open " << logfile << " for logging.\n";
		throw 0;
	}
	(*logfilestream) << logbuf.str();
	logfilestream->flush();
	logstream = logfilestream;

	committed = true;

	//If anything was written to the error stream before now, materialise the
	//file; otherwise leave it absent until the first error.
	if (!errbuf.str().empty()) errOut();

	return;
}

//Opens <out>.error on first use after commit() and replays whatever had been
//buffered.  Before commit() this just returns the in-memory buffer.
ostream *errlog::errOut()
{
	if (committed && errstream == &errbuf)
	{
		errfilestream = new ofstream;
		errfilestream->open(errfile.c_str(), ios::binary);
		if (errfilestream->fail())
		{
			cerr << "ERROR: Could not open " << errfile << " for logging.\n";
			throw 0;
		}
		(*errfilestream) << errbuf.str();
		errstream = errfilestream;
	}
	return errstream;
}

void errlog::errn(string str)
{
	this->outn(&cerr, str);
	this->outn(errOut(), str);
	return;
}

void errlog::err(string str)
{
	this->out(&cerr, str);
	this->out(errOut(), str);
	return;
}

void errlog::err(double val)
{
	this->out(&cerr, val);
	this->out(errOut(), val);
	return;
}

void errlog::errn(double val)
{
	this->out(&cerr, val);
	this->out(errOut(), val);
	return;
}

void errlog::err(int val)
{
	this->out(&cerr, val);
	this->out(errOut(), val);
	return;
}

void errlog::errn(int val)
{
	this->out(&cerr, val);
	this->out(errOut(), val);
	return;
}

void errlog::err(char val)
{
	this->out(&cerr, val);
	this->out(errOut(), val);
	return;
}

void errlog::errn(char val)
{
	this->out(&cerr, val);
	this->out(errOut(), val);
	return;
}

void errlog::err(string str, int val, bool nl)
{
	this->out(&cerr, str, val, nl);
	this->out(errOut(), str, val, nl);
	return;
}

void errlog::err(string str, long val, bool nl)
{
	//See the header: int64_t is long on LP64 and long long on LLP64, so
	//both need an exact match.  Same string routing as the long long form.
	stringstream ss;
	ss << val;
	this->err(str, ss.str(), nl);
	return;
}

void errlog::err(string str, long long val, bool nl)
{
	//Formatted through the string overload rather than adding another out():
	//the int overload would narrow a 64-bit position and the double overload
	//would render it in scientific notation.
	stringstream ss;
	ss << val;
	this->err(str, ss.str(), nl);
	return;
}

void errlog::errv(string str, vector<int> &val, bool nl)
{
	this->outv(&cerr, str, val, nl);
	this->outv(errOut(), str, val, nl);
	return;
}

void errlog::err(string str, double val, bool nl)
{
	this->out(&cerr, str, val, nl);
	this->out(errOut(), str, val, nl);
	return;
}

void errlog::errv(string str, vector<double> &val, bool nl)
{
	this->outv(&cerr, str, val, nl);
	this->outv(errOut(), str, val, nl);
	return;
}

void errlog::err(string str, bool val, bool nl)
{
	this->out(&cerr, str, val, nl);
	this->out(errOut(), str, val, nl);
	return;
}

void errlog::err(string str, string val, bool nl)
{
	this->out(&cerr, str, val, nl);
	this->out(errOut(), str, val, nl);
	return;
}

void errlog::err(string str, char val, bool nl)
{
	this->out(&cerr, str, val, nl);
	this->out(errOut(), str, val, nl);
	return;
}

void errlog::erra(string str, string *val, int size, bool nl)
{
	this->outa(&cerr, str, val, size, nl);
	this->outa(errOut(), str, val, size, nl);
	return;
}

void errlog::erra(string str, int *val, int size, bool nl)
{
	this->outa(&cerr, str, val, size, nl);
	this->outa(errOut(), str, val, size, nl);
	return;
}

void errlog::erra(string str, double *val, int size, bool nl)
{
	this->outa(&cerr, str, val, size, nl);
	this->outa(errOut(), str, val, size, nl);
	return;
}

void errlog::erra(string str, char *val, int size, bool nl)
{
	this->outa(&cerr, str, val, size, nl);
	this->outa(errOut(), str, val, size, nl);
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

void errlog::log(string str, long val, bool nl)
{
	//See the header: int64_t is long on LP64 and long long on LLP64, so
	//both need an exact match.  Same string routing as the long long form.
	stringstream ss;
	ss << val;
	this->log(str, ss.str(), nl);
	return;
}

void errlog::log(string str, long long val, bool nl)
{
	//See err(string, long long, bool): routed through the string overload so a
	//64-bit position is neither narrowed nor printed in scientific notation.
	stringstream ss;
	ss << val;
	this->log(str, ss.str(), nl);
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

void logCurrentException(const string &stage)
{
    try
    {
        throw;                      //re-throw whatever is being handled
    }
    catch (const std::exception &e)
    {
        LOG.err("ERROR: " + stage + " failed:", string(e.what()));
    }
    catch (int code)
    {
        //garlic's own signalled failures.  The site that threw has already
        //reported the cause, so only the stage is added here.
        LOG.err("ERROR: " + stage + " failed, code", code);
    }
    catch (...)
    {
        LOG.err("ERROR: " + stage + " failed with an unrecognised exception.");
    }
    return;
}
