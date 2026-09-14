/* param_t -- a class for basic command line argument parsing
   Copyright (C) 2014  Zachary A Szpiech

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
*/
#include "param_t.h"
#include <iterator>
#include <cerrno>

using namespace std;

const string ARG_HELP = "--help";
const string ARG_HELP_SHORT = "-h";

//Several flags -- list flags in particular -- state their real default in the
//description, because the auto-generated one can only show a single sentinel
//element ("0.000000", "_ALL").  Appending a second Default line in that case
//made --help print both.
static string defaultSuffix(const string &description, const string &buffer)
{
    if (description.find("Default:") != string::npos) return "";
    return "\n\tDefault: " + buffer;
}

bool param_t::addFlag(string flag, bool value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        if (value) buffer = "true";
        else buffer = "false";
        argb[flag] = value;
        help[flag] = "<bool>: " + description + defaultSuffix(description, buffer);
        labels[flag] = label;
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addFlag(string flag, double value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        char charBuffer[100];
        //%g rather than %e: the default for --auto-winsize-threshold used to
        //print as "5.000000e-01" in --help and in the generated manual.
        sprintf(charBuffer, "%g", value);
        buffer = charBuffer;
        argd[flag] = value;
        help[flag] = "<double>: " + description + defaultSuffix(description, buffer);
        labels[flag] = label;
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addFlag(string flag, int value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        char charBuffer[100];
        sprintf(charBuffer, "%d", value);
        buffer = charBuffer;
        argi[flag] = value;
        help[flag] = "<int>: " + description + defaultSuffix(description, buffer);
        labels[flag] = label;
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addFlag(string flag, char value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        char charBuffer[100];
        sprintf(charBuffer, "%c", value);
        buffer = charBuffer;
        argch[flag] = value;
        help[flag] = "<char>: " + description + defaultSuffix(description, buffer);
        labels[flag] = label;
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addFlag(string flag, string value, string label, string description)
{
    if (!flagExists(flag))
    {
        args[flag] = value;
        help[flag] = "<string>: " + description + defaultSuffix(description, value);
        labels[flag] = label;
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addFlag(string flag, const char value[], string label, string description)
{
    return this->addFlag(flag, string(value), label, description);
}

bool param_t::addListFlag(string flag, double value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        char charBuffer[100];
        sprintf(charBuffer, "%f", value);
        buffer = charBuffer;
        listargd[flag].push_back(value);
        help[flag] = "<double1> ... <doubleN>: " + description + defaultSuffix(description, buffer);
        labels[flag] = label;
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}
bool param_t::addListFlag(string flag, char value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        char charBuffer[100];
        sprintf(charBuffer, "%c", value);
        buffer = charBuffer;
        listargch[flag].push_back(value);
        help[flag] = "<char1> ... <charN>: " + description + defaultSuffix(description, buffer);
        labels[flag] = label;
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addListFlag(string flag, int value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        char charBuffer[100];
        sprintf(charBuffer, "%d", value);
        buffer = charBuffer;
        listargi[flag].push_back(value);
        help[flag] = "<int1> ... <intN>: " + description + defaultSuffix(description, buffer);
        labels[flag] = label;
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool  param_t::addListFlag(string flag, string value, string label, string description)
{
    if (!flagExists(flag))
    {
        listargs[flag].push_back(value);
        help[flag] = "<string1> ... <stringN>: " + description + defaultSuffix(description, value);
        labels[flag] = label;
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addListFlag(string flag, const char value[], string label, string description)
{
    return this->addListFlag(flag, string(value), label, description);
}

void param_t::printHelp()
{
    map<string, string>::iterator it;

    cerr << preamble << endl;

    cerr << "----------Command Line Arguments----------\n\n";

    for (it = help.begin(); it != help.end(); it++)
    {
        if (labels[(*it).first].compare("SILENT") != 0)
        {
            cerr << (*it).first << " " << (*it).second << "\n\n";
        }
    }

    return;
}

bool param_t::goodDouble(string str)
{
    //Hand-rolled validation rejected scientific notation, so --size-bounds 1e6
    //-- a natural way to write a ROH size threshold -- was reported as 'not a
    //valid double' and then as an unrecognised flag.  Defer to strtod, which
    //accepts exactly what atof will subsequently parse.
    if (str.empty()) return 0;
    const char *s = str.c_str();
    char *end = NULL;
    errno = 0;
    strtod(s, &end);
    if (end == s) return 0;          //nothing consumed
    if (*end != '\0') return 0;      //trailing junk
    return 1;
}

bool param_t::goodInt(string str)
{
    string::iterator it;
    //int dashCount = 0;
    for (it = str.begin(); it != str.end(); it++)
    {
        if (!isdigit(*it) && *it != '-') return 0;
        if (*it == '-' && it != str.begin()) return 0;
        //if (dashCount > 1) return 0;
    }
    return 1;
}

bool param_t::goodChar(string str)
{
    if (str.length() > 1) return 0;
    return 1;
}

//help[flag] is built by addFlag as "<TYPE>: DESCRIPTION\n\tDefault: VALUE".
//Split it back into those three parts.
static void splitHelpEntry(const string &entry, string &type, string &desc, string &def)
{
    type.clear(); desc.clear(); def.clear();
    size_t lt = entry.find('<'), gt = entry.find(">:");
    size_t descStart = 0;
    if (lt == 0 && gt != string::npos)
    {
        type = entry.substr(1, gt - 1);
        descStart = gt + 2;
        while (descStart < entry.size() && entry[descStart] == ' ') descStart++;
    }
    const string marker = "\n\tDefault: ";
    size_t d = entry.find(marker, descStart);
    if (d == string::npos) { desc = entry.substr(descStart); return; }
    desc = entry.substr(descStart, d - descStart);
    def  = entry.substr(d + marker.size());
    size_t again = def.find("Default:");
    if (again != string::npos) def = def.substr(0, again);
}

//The HELP_ strings wrap with "\n\t" for terminal output; documentation
//formats re-flow, so collapse those to single spaces.
static string unwrap(const string &s)
{
    string o;
    for (size_t i = 0; i < s.size(); i++)
    {
        if (s[i] == '\n' || s[i] == '\t') { if (!o.empty() && o[o.size()-1] != ' ') o += ' '; }
        else o += s[i];
    }
    while (!o.empty() && o[o.size()-1] == ' ') o.erase(o.size()-1);
    return o;
}

//A bare "--flag" in a description is typeset as an en-dash followed by the
//name in LaTeX text mode, so cross-references to other flags came out as
//"-weighted" in the manual.  Wrap each one in \\texttt{} instead, which also
//makes them look like the flags they are.  Runs after texEscape, so the input
//is already escaped.
static string texWrapFlags(const string &s)
{
    string o;
    size_t i = 0;
    while (i < s.size())
    {
        bool atStart = (i == 0) || s[i-1] == ' ' || s[i-1] == '(' || s[i-1] == ',';
        if (atStart && s.compare(i, 2, "--") == 0 && i + 2 < s.size() && isalpha(s[i+2]))
        {
            size_t j = i + 2;
            while (j < s.size() && (isalnum(s[j]) || s[j] == '-')) j++;
            o += "\\texttt{--" + s.substr(i + 2, j - (i + 2)) + "}";
            i = j;
            continue;
        }
        o += s[i];
        i++;
    }
    return o;
}

static string texEscape(const string &s)
{
    string o;
    for (size_t i = 0; i < s.size(); i++)
    {
        char c = s[i];
        switch (c)
        {
        case '\\': o += "\\textbackslash{}"; break;
        case '{': case '}': case '$': case '&': case '#': case '%': case '_':
            o += '\\'; o += c; break;
        case '^': o += "\\textasciicircum{}"; break;
        case '~': o += "\\textasciitilde{}"; break;
        case '<': o += "\\textless{}"; break;
        case '>': o += "\\textgreater{}"; break;
        default: o += c;
        }
    }
    return o;
}

bool param_t::writeHelpDoc(ostream &out, string format)
{
    if (format != "txt" && format != "tex")
    {
        cerr << "ERROR: unknown documentation format '" << format << "'. Use txt or tex.\n";
        return false;
    }

    if (format == "tex") out << "% Generated by 'make docs' from the HELP_ strings in garlic-cli.cpp.\n"
                             //Descriptions contain unbreakable tokens such as
                             //slope*log(density), which overflowed the text block
                             //(10 overfull hboxes) without \sloppy.
                             << "% Do not edit by hand.\n\\begingroup\\sloppy\n\\begin{description}\n";
    else out << "This section is generated by 'make docs' from the program's own help\n"
             << "strings. Do not edit by hand.\n\n";

    map<string, string>::iterator it;
    for (it = help.begin(); it != help.end(); it++)
    {
        if (labels[it->first].compare("SILENT") == 0) continue;
        string type, desc, def;
        splitHelpEntry(it->second, type, desc, def);
        if (format == "txt")
        {
            out << it->first << " <" << type << ">: " << unwrap(desc) << "\n";
            if (!def.empty()) out << "\tDefault: " << unwrap(def) << "\n";
            out << "\n";
        }
        else
        {
            out << "\\item[\\texttt{" << texEscape(it->first) << "}";
            if (!type.empty()) out << " \\textnormal{\\textless{}" << texEscape(type) << "\\textgreater{}}";
            out << "] " << texWrapFlags(texEscape(unwrap(desc)));
            if (!def.empty()) out << " \\\\ \\textit{Default:} \\texttt{" << texEscape(unwrap(def)) << "}";
            out << "\n";
        }
    }
    if (format == "tex") out << "\\end{description}\n\\endgroup\n";
    return true;
}

int param_t::parseCommandLine(int argc, char *argv[])
{
    int badFlags = 0;

    for (int i = 1; i < argc; i++)
    {
        if (isSet.count(argv[i]) > 0)
        {
            cerr << "ERROR: Duplicate " << argv[i] << " found.\n";
            badFlags++;
            break;
        }
        else if (argb.count(argv[i]) > 0)
        {
            //Set, do not toggle: toggling silently inverts any flag whose
            //default is true.
            argb[argv[i]] = true;
            isSet[argv[i]] = true;
        }
        else if (argi.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else if (!goodInt(string(argv[i + 1])))
            {
                cerr << "ERROR: " << argv[i + 1] << " is not a valid integer.\n";
                badFlags++;
                break;
            }
            else
            {
                argi[argv[i]] = atoi(argv[i + 1]);
                isSet[argv[i]] = true;
                i++;
            }
        }
        else if (listargi.count(argv[i]) > 0)
        {
            //List flags never recorded isSet, so isFlagSet() reported them unset
            //and the duplicate-flag check never fired for them.
            isSet[argv[i]] = true;
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else
            {
                listargi[argv[i]].clear();//clear the default value
                int flagIndex = i;//remember where the flag is in argv
                while (i + 1 < argc) //go until the end of the list
                {
                    if (goodInt(string(argv[i + 1]))) //make sure the next value is OK
                    {
                        listargi[argv[flagIndex]].push_back(atoi(argv[i + 1]));
                        i++;
                    }
                    //if it is a bad int...
                    else if (!goodInt(string(argv[i + 1])) && !flagExists(string(argv[i + 1])))
                    {
                        cerr << "ERROR: " << argv[i + 1] << " is not a valid integer.\n";
                        badFlags++;
                        break;
                    }
                    else //if the next value is another flag..
                    {
                        if (listargi[argv[flagIndex]].size() == 0)
                        {
                            cerr << "ERROR: No arguments found for " << argv[flagIndex] << ".\n";
                            badFlags++;
                        }
                        break;
                    }
                }
            }
        }
        else if (argd.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else if (!goodDouble(string(argv[i + 1])))
            {
                cerr << "ERROR: " << argv[i + 1] << " is not a valid double.\n";
                badFlags++;
                break;
            }
            else
            {
                argd[argv[i]] = atof(argv[i + 1]);
                isSet[argv[i]] = true;
                i++;
            }
        }
        else if (listargd.count(argv[i]) > 0)
        {
            //List flags never recorded isSet, so isFlagSet() reported them unset
            //and the duplicate-flag check never fired for them.
            isSet[argv[i]] = true;
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else
            {
                listargd[argv[i]].clear();//clear the default value
                int flagIndex = i;//remember where the flag is in argv
                while (i + 1 < argc) //go until the end of the list
                {
                    if (goodDouble(string(argv[i + 1]))) //make sure the next value is OK
                    {
                        listargd[argv[flagIndex]].push_back(atof(argv[i + 1]));
                        i++;
                    }
                    //if it is a bad int...
                    else if (!goodDouble(string(argv[i + 1])) && !flagExists(string(argv[i + 1])))
                    {
                        cerr << "ERROR: " << argv[i + 1] << " is not a valid double.\n";
                        badFlags++;
                        break;
                    }
                    else //if the next value is another flag..
                    {
                        if (listargd[argv[flagIndex]].size() == 0)
                        {
                            cerr << "ERROR: No arguments found for " << argv[flagIndex] << ".\n";
                            badFlags++;
                        }
                        break;
                    }
                }
            }
        }
        else if (args.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else
            {
                args[argv[i]] = argv[i + 1];
                isSet[argv[i]] = true;
                i++;
            }
        }
        else if (listargs.count(argv[i]) > 0)
        {
            //List flags never recorded isSet, so isFlagSet() reported them unset
            //and the duplicate-flag check never fired for them.
            isSet[argv[i]] = true;
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else
            {
                listargs[argv[i]].clear();//clear the default value
                int flagIndex = i;//remember where the flag is in argv
                while (i + 1 < argc) //go until the end of the list
                {
                    if (argv[i + 1][0] != '-') //make sure the next value isn't another flag
                    {
                        listargs[argv[flagIndex]].push_back(string(argv[i + 1]));
                        i++;
                    }
                    else //if the next value is another flag...
                    {
                        if (listargs[argv[flagIndex]].size() == 0)
                        {
                            cerr << "ERROR: No arguments found for " << argv[flagIndex] << ".\n";
                            badFlags++;
                        }
                        break;
                    }
                }
            }
        }
        else if (argch.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else if (!goodChar(string(argv[i + 1])))
            {
                cerr << "ERROR: " << argv[i + 1] << " is not a valid character.\n";
                badFlags++;
                break;
            }
            else
            {
                argch[argv[i]] = argv[i + 1][0];
                isSet[argv[i]] = true;
                i++;
            }
        }
        else if (listargch.count(argv[i]) > 0)
        {
            //List flags never recorded isSet, so isFlagSet() reported them unset
            //and the duplicate-flag check never fired for them.
            isSet[argv[i]] = true;
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else
            {
                listargch[argv[i]].clear();//clear the default value
                int flagIndex = i;//remember where the flag is in argv
                while (i + 1 < argc) //go until the end of the list
                {
                    if (goodChar(string(argv[i + 1]))) //make sure the next value is OK
                    {
                        listargch[argv[flagIndex]].push_back(atoi(argv[i + 1]));
                        i++;
                    }
                    //if it is a bad int...
                    else if (!goodChar(string(argv[i + 1])) && !flagExists(string(argv[i + 1])))
                    {
                        cerr << "ERROR: " << argv[i + 1] << " is not a valid character.\n";
                        badFlags++;
                        break;
                    }
                    else //if the next value is another flag..
                    {
                        if (listargch[argv[flagIndex]].size() == 0)
                        {
                            cerr << "ERROR: No arguments found for " << argv[flagIndex] << ".\n";
                            badFlags++;
                        }
                        break;
                    }
                }
            }
        }
        else //if (argv[i][0] == '-')
        {
            cerr << "ERROR: command line flag " << argv[i] << " not recognized.\n";
            badFlags++;
        }
    }

    if (getBoolFlag(ARG_HELP) || getBoolFlag(ARG_HELP_SHORT))
    {
        this->printHelp();
        return PARAM_HELP;
    }

    if (badFlags) return PARAM_ERROR;

    return PARAM_OK;
}

bool param_t::flagExists(string flag)
{
    return (help.count(flag) > 0);
}

param_t::param_t()
{
    this->addFlag(ARG_HELP, false, "__help", "Prints this help dialog.");
    this->addFlag(ARG_HELP_SHORT, false, "__help", "Prints this help dialog.");
}

bool param_t::isFlagSet(string flag)
{
    return (isSet.count(flag) > 0);
}

static string jsonEscape(const string &s)
{
    string o;
    for (unsigned int i = 0; i < s.size(); i++)
    {
        char c = s[i];
        if (c == '"' || c == '\\') { o += '\\'; o += c; }
        else if (c == '\n') o += "\\n";
        else if (c == '\t') o += "\\t";
        else o += c;
    }
    return o;
}

vector<string> param_t::setFlags()
{
    vector<string> out;
    map<string, bool>::iterator it;
    for (it = isSet.begin(); it != isSet.end(); it++)
        if (it->second) out.push_back(it->first);
    return out;
}

bool param_t::setIntFlag(string flag, int value)
{
    if (argi.count(flag) == 0) return false;
    argi[flag] = value;
    //Mark it supplied, so it appears in the params record's "set" list and a
    //replay uses this value rather than drawing a fresh one.
    isSet[flag] = true;
    return true;
}

void param_t::writeFlagsJSON(ostream &out, string indent)
{
    bool first = true;
    map<string, bool>::iterator bi;
    for (bi = argb.begin(); bi != argb.end(); bi++)
    {
        if (bi->first == ARG_HELP || bi->first == ARG_HELP_SHORT) continue;
        if (!first) out << ",\n";
        out << indent << "\"" << bi->first << "\": " << (bi->second ? "true" : "false");
        first = false;
    }
    map<string, int>::iterator ii;
    for (ii = argi.begin(); ii != argi.end(); ii++)
    { if (!first) out << ",\n"; out << indent << "\"" << ii->first << "\": " << ii->second; first = false; }
    map<string, double>::iterator di;
    for (di = argd.begin(); di != argd.end(); di++)
    { if (!first) out << ",\n"; out << indent << "\"" << di->first << "\": " << di->second; first = false; }
    map<string, string>::iterator si;
    for (si = args.begin(); si != args.end(); si++)
    { if (!first) out << ",\n"; out << indent << "\"" << si->first << "\": \"" << jsonEscape(si->second) << "\""; first = false; }
    map<string, char>::iterator ci;
    for (ci = argch.begin(); ci != argch.end(); ci++)
    { if (!first) out << ",\n"; out << indent << "\"" << ci->first << "\": \"" << ci->second << "\""; first = false; }

    map<string, vector<int> >::iterator li;
    for (li = listargi.begin(); li != listargi.end(); li++)
    {
        if (!first) out << ",\n";
        out << indent << "\"" << li->first << "\": [";
        for (unsigned int k = 0; k < li->second.size(); k++) { if (k) out << ", "; out << li->second[k]; }
        out << "]"; first = false;
    }
    map<string, vector<double> >::iterator ld;
    for (ld = listargd.begin(); ld != listargd.end(); ld++)
    {
        if (!first) out << ",\n";
        out << indent << "\"" << ld->first << "\": [";
        for (unsigned int k = 0; k < ld->second.size(); k++) { if (k) out << ", "; out << ld->second[k]; }
        out << "]"; first = false;
    }
    map<string, vector<string> >::iterator ls;
    for (ls = listargs.begin(); ls != listargs.end(); ls++)
    {
        if (!first) out << ",\n";
        out << indent << "\"" << ls->first << "\": [";
        for (unsigned int k = 0; k < ls->second.size(); k++) { if (k) out << ", "; out << "\"" << jsonEscape(ls->second[k]) << "\""; }
        out << "]"; first = false;
    }
    if (!first) out << "\n";
    return;
}

//Minimal reader for the object this program writes: a flat set of
//"--flag": value members, where value is a number, string, boolean or array of
//those.  Nested objects (the "resolved" block) are skipped.
bool param_t::loadFlagsJSON(string file)
{
    ifstream fin(file.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Could not open " << file << " for reading.\n";
        return false;
    }
    string text((istreambuf_iterator<char>(fin)), istreambuf_iterator<char>());
    fin.close();

    //Read the "set" allow-list first, if the file has one.
    vector<string> allow;
    bool haveAllow = false;
    {
        size_t a = text.find("\"set\"");
        if (a != string::npos)
        {
            size_t lb = text.find('[', a);
            size_t rb = (lb == string::npos) ? string::npos : text.find(']', lb);
            if (lb != string::npos && rb != string::npos)
            {
                haveAllow = true;
                string body = text.substr(lb + 1, rb - lb - 1);
                string cur;
                for (unsigned int k = 0; k <= body.size(); k++)
                {
                    if (k == body.size() || body[k] == ',')
                    {
                        while (!cur.empty() && (cur[0] == ' ' || cur[0] == '"' || cur[0] == '\n' || cur[0] == '\r' || cur[0] == '\t')) cur.erase(cur.begin());
                        while (!cur.empty() && (cur[cur.size()-1] == ' ' || cur[cur.size()-1] == '"' || cur[cur.size()-1] == '\n' || cur[cur.size()-1] == '\r' || cur[cur.size()-1] == '\t')) cur.erase(cur.size()-1);
                        if (!cur.empty()) allow.push_back(cur);
                        cur.clear();
                    }
                    else cur += body[k];
                }
            }
        }
    }

    size_t i = 0;
    const size_t n = text.size();
    int applied = 0;

    while (i < n)
    {
        while (i < n && text[i] != '"') i++;
        if (i >= n) break;
        size_t ks = ++i;
        while (i < n && text[i] != '"') i++;
        if (i >= n) break;
        string key = text.substr(ks, i - ks);
        i++;
        while (i < n && (text[i] == ' ' || text[i] == '\t' || text[i] == '\n' || text[i] == '\r')) i++;
        if (i >= n || text[i] != ':') continue;
        i++;
        while (i < n && (text[i] == ' ' || text[i] == '\t' || text[i] == '\n' || text[i] == '\r')) i++;
        if (i >= n) break;

        //Collect the raw value.
        string val;
        if (text[i] == '{')
        {
            //Descend into the flags object; skip any other nested object
            //(currently "resolved", which is a record, not input).
            if (key == "flags") { i++; continue; }
            int depth = 0;
            while (i < n) { if (text[i] == '{') depth++; else if (text[i] == '}') { depth--; if (!depth) { i++; break; } } i++; }
            continue;
        }
        else if (text[i] == '[')
        {
            if (key == "set") { while (i < n && text[i] != ']') i++; if (i < n) i++; continue; }
            size_t vs = i++;
            while (i < n && text[i] != ']') i++;
            val = text.substr(vs + 1, i - vs - 1);
            if (i < n) i++;
        }
        else if (text[i] == '"')
        {
            size_t vs = ++i;
            while (i < n && text[i] != '"') i++;
            val = text.substr(vs, i - vs);
            if (i < n) i++;
        }
        else
        {
            size_t vs = i;
            while (i < n && text[i] != ',' && text[i] != '}' && text[i] != '\n') i++;
            val = text.substr(vs, i - vs);
            while (!val.empty() && (val[val.size()-1] == ' ' || val[val.size()-1] == '\r')) val.erase(val.size()-1);
        }

        if (key.size() < 2 || key[0] != '-') continue;   //not a flag
        if (isSet.count(key) > 0) continue;              //command line wins
        if (haveAllow)
        {
            bool ok = false;
            for (unsigned int k = 0; k < allow.size(); k++) if (allow[k] == key) { ok = true; break; }
            if (!ok) continue;                            //recorded, but not user-supplied
        }

        bool known = true;
        if (argb.count(key) > 0)            argb[key] = (val == "true" || val == "1");
        else if (argi.count(key) > 0)       argi[key] = atoi(val.c_str());
        else if (argd.count(key) > 0)       argd[key] = atof(val.c_str());
        else if (args.count(key) > 0)       args[key] = val;
        else if (argch.count(key) > 0)      argch[key] = val.empty() ? ' ' : val[0];
        else if (listargi.count(key) > 0 || listargd.count(key) > 0 || listargs.count(key) > 0)
        {
            vector<string> toks; string cur;
            for (unsigned int k = 0; k <= val.size(); k++)
            {
                if (k == val.size() || val[k] == ',')
                {
                    while (!cur.empty() && (cur[0] == ' ' || cur[0] == '"')) cur.erase(cur.begin());
                    while (!cur.empty() && (cur[cur.size()-1] == ' ' || cur[cur.size()-1] == '"')) cur.erase(cur.size()-1);
                    if (!cur.empty()) toks.push_back(cur);
                    cur.clear();
                }
                else cur += val[k];
            }
            if (toks.empty()) { known = false; }
            else if (listargi.count(key) > 0)
            { listargi[key].clear(); for (unsigned int k = 0; k < toks.size(); k++) listargi[key].push_back(atoi(toks[k].c_str())); }
            else if (listargd.count(key) > 0)
            { listargd[key].clear(); for (unsigned int k = 0; k < toks.size(); k++) listargd[key].push_back(atof(toks[k].c_str())); }
            else
            { listargs[key].clear(); for (unsigned int k = 0; k < toks.size(); k++) listargs[key].push_back(toks[k]); }
        }
        else known = false;

        if (known) { isSet[key] = true; applied++; }
    }

    cerr << "Loaded " << applied << " parameters from " << file << ".\n";
    return true;
}

bool param_t::getBoolFlag(string flag)
{
    if (argb.count(flag) > 0) return argb[flag];

    cerr << "ERROR: There are no bool flags named " << flag << "\n";
    throw 0;
}

double param_t::getDoubleFlag(string flag)
{
    if (argd.count(flag) > 0) return argd[flag];

    cerr << "ERROR: There are no double flags named " << flag << "\n";
    throw 0;
}

int param_t::getIntFlag(string flag)
{
    if (argi.count(flag) > 0) return argi[flag];

    cerr << "ERROR: There are no int flags named " << flag << "\n";
    throw 0;
}

char param_t::getCharFlag(string flag)
{
    if (argch.count(flag) > 0) return argch[flag];

    cerr << "ERROR: There are no char flags named " << flag << "\n";
    throw 0;
}

string param_t::getStringFlag(string flag)
{
    if (args.count(flag) > 0) return args[flag];

    cerr << "ERROR: There are no string flags named " << flag << "\n";
    throw 0;
}

vector<string> param_t::getStringListFlag(string flag)
{
    if (listargs.count(flag) > 0) return listargs[flag];

    cerr << "ERROR: There are no string list flags named " << flag << "\n";
    throw 0;
}

vector<int> param_t::getIntListFlag(string flag)
{
    if (listargi.count(flag) > 0) return listargi[flag];

    cerr << "ERROR: There are no int list flags named " << flag << "\n";
    throw 0;
}

vector<double> param_t::getDoubleListFlag(string flag)
{
    if (listargd.count(flag) > 0) return listargd[flag];

    cerr << "ERROR: There are no double list flags named " << flag << "\n";
    throw 0;
}

vector<char> param_t::getCharListFlag(string flag)
{
    if (listargch.count(flag) > 0) return listargch[flag];

    cerr << "ERROR: There are no int list flags named " << flag << "\n";
    throw 0;
}

void param_t::setPreamble(string str)
{
    preamble = str;
    return;
}
