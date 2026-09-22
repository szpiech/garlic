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
#ifndef __PARAM_T_H__
#define __PARAM_T_H__

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cctype>
#include <cstdlib>
#include <cstdio>

using namespace std;

extern const string ARG_HELP;
extern const string ARG_HELP_SHORT;

//parseCommandLine outcomes.  It used to return 0 for BOTH --help and a bad
//command line, so the caller could not tell success from a usage error and
//main exited 0 on a typo'd flag.
const int PARAM_OK    =  1;
const int PARAM_HELP  =  0;
const int PARAM_ERROR = -1;

class param_t
{
public:

  bool addFlag(string flag, bool value, string label, string description);
  bool addFlag(string flag, double value, string label, string description);
  bool addFlag(string flag, int value, string label, string description);
  bool addFlag(string flag, char value, string label, string description);
  bool addFlag(string flag, string value, string label, string description);
  bool addFlag(string flag, const char value[], string label, string description);

  bool addListFlag(string flag, string value, string label, string description);
  bool addListFlag(string flag, const char value[], string label, string description);
  bool addListFlag(string flag, int value, string label, string description);
  bool addListFlag(string flag, double value, string label, string description);
  bool addListFlag(string flag, char value, string label, string description);

  //Help output is grouped rather than alphabetical: 77 flags in one flat list
  //is a reference you can search, not one you can read, and a newcomer cannot
  //tell from it which flags they must set and which have defaults that will do.
  //
  //The taxonomy is declared ONCE by the caller (see setHelpCategories in
  //garlic-cli.cpp) rather than as an argument on each addFlag, so the whole
  //grouping can be read and reordered in one place instead of being recovered
  //by grepping seventy call sites.  The cost of that choice is that a new flag
  //has to be named in two places, so it is enforced rather than trusted: a
  //non-SILENT flag in no category, a category naming a flag that does not
  //exist, and a flag named twice are all errors.  --help says so and
  //--dump-docs fails, which makes the docs job in CI the backstop.
  //
  //Leaving this unset keeps the old flat alphabetical output.
  void setHelpCategories(const vector< pair< string, vector<string> > > &cats);

  void printHelp();

  //Emit the flag reference in a documentation format, so README and the
  //manual do not have to restate by hand what the HELP_ strings already say.
  //format is "txt" (README block) or "tex" (LaTeX description list).
  //Driven by the same help/labels maps printHelp uses, so a flag cannot be
  //added to the program without appearing here.
  bool writeHelpDoc(ostream &out, string format);

  int parseCommandLine(int argc, char *argv[]);

  //Was the flag given on the command line?  Lets callers stop inferring
  //"unset" from magic default values like --lod-cutoff -999999.
  bool isFlagSet(string flag);

  //Dump every flag and its effective value as JSON object members (no braces),
  //and read such a file back.  loadFlagsJSON only fills flags that were NOT
  //given on the command line, so an explicit argument always wins over the
  //file.
  //Used to record the seed that was actually drawn, so the flags block of
  //<out>.params.json is a command line that reproduces the run.
  bool setIntFlag(string flag, int value);

  //Flags explicitly supplied (on the command line, or loaded from a params
  //file).  The params record dumps every flag's effective value for the
  //reader, but a replay must apply only these -- several flags carry sentinel
  //defaults that fail their own validation if presented as user input.
  vector<string> setFlags();

  void writeFlagsJSON(ostream &out, string indent);
  bool loadFlagsJSON(string file);

  bool getBoolFlag(string flag);
  double getDoubleFlag(string flag);
  int getIntFlag(string flag);
  char getCharFlag(string flag);
  string getStringFlag(string flag);

  vector<string> getStringListFlag(string flag);
  vector<int> getIntListFlag(string flag);
  vector<double> getDoubleListFlag(string flag);
  vector<char> getCharListFlag(string flag);

  void setPreamble(string str);

  param_t();


private:

  map<string, bool> argb;
  map<string, double> argd;
  map<string, int> argi;
  map<string, char> argch;
  map<string, string> args;

  map<string, vector< string > > listargs;
  map<string, vector< int > > listargi;
  map<string, vector< double > > listargd;
  map<string, vector< char > > listargch;

  map<string, string> help;
  map<string, bool> isSet;
  map<string, string> labels;

  //Declared order of the help groups, each with the flags it holds.
  vector< pair< string, vector<string> > > helpCats;

  //helpCats resolved against the registered flags: SILENT ones dropped and
  //each group sorted by flag name, so the grouping is editorial and the order
  //within a group is not.  Every complaint -- an uncategorised flag, an
  //unknown one, a duplicate -- goes into `problems`; an empty helpCats yields
  //one unnamed group holding everything, which is the old flat output.
  vector< pair< string, vector<string> > > groupedHelp(vector<string> &problems) const;

  bool goodDouble(string str);
  bool goodInt(string str);
  bool goodChar(string str);

  bool flagExists(string flag);

  string preamble;
};

#endif
