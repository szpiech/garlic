#ifndef __GARLIC_PLATFORM_H__
#define __GARLIC_PLATFORM_H__

//The three places where POSIX and the Windows CRT genuinely differ, kept
//together rather than as #ifdefs at the call sites.

#include <sys/stat.h>  //stat/_stat, S_IFMT, S_IFDIR

#ifdef _WIN32
  #include <direct.h>    //_mkdir
  #include <io.h>        //_isatty
#else
  #include <unistd.h>    //isatty, STDERR_FILENO
#endif

//mkdir is not portable.  POSIX takes (path, mode); the Windows CRT takes
//(path) only -- MinGW declares `int mkdir(const char *)` and deprecates it in
//favour of _mkdir, so passing a mode is "too many arguments to function".
//That is exactly how the first MinGW build of --outdir failed.  The mode has
//no meaning on Windows anyway: NTFS ACLs are not a umask.
static inline int garlicMkdir(const char *path)
{
#ifdef _WIN32
    return _mkdir(path);
#else
    return mkdir(path, 0777);
#endif
}

//isatty and STDERR_FILENO come from <unistd.h> on POSIX.  MinGW does ship a
//<unistd.h>, but the supported spellings there are _isatty and a literal
//descriptor number, so use those rather than depend on its compatibility
//layer being present -- this cannot be verified without a Windows toolchain,
//and the explicit form costs nothing.
static inline bool garlicStderrIsTTY()
{
#ifdef _WIN32
    return _isatty(2) != 0;
#else
    return isatty(STDERR_FILENO) != 0;
#endif
}

//Does this path exist AND is it a directory?
//
//mkdir returning EEXIST is not enough to conclude the output directory is
//usable: a regular FILE of that name gives EEXIST too.  garlic used to accept
//that and then abort when the writers could not open their outputs -- exit 134
//on POSIX for `--outdir <an existing file>`, and on Windows it made an
//uncreatable path look like success.  Tested with the S_IFMT mask rather than
//S_ISDIR, which POSIX guarantees but the Windows CRT does not.
static inline bool garlicIsDir(const char *path)
{
#ifdef _WIN32
    struct _stat st;
    if (_stat(path, &st) != 0) return false;
#else
    struct stat st;
    if (stat(path, &st) != 0) return false;
#endif
    return (st.st_mode & S_IFMT) == S_IFDIR;
}

//Windows accepts either separator in a path; POSIX treats a backslash as an
//ordinary filename character, so it must NOT be split on there.
static inline bool garlicIsPathSep(char c)
{
#ifdef _WIN32
    return c == '/' || c == '\\';
#else
    return c == '/';
#endif
}

#endif
