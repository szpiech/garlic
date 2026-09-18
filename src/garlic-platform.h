#ifndef __GARLIC_PLATFORM_H__
#define __GARLIC_PLATFORM_H__

//The three places where POSIX and the Windows CRT genuinely differ, kept
//together rather than as #ifdefs at the call sites.

#ifdef _WIN32
  #include <direct.h>    //_mkdir
  #include <io.h>        //_isatty
#else
  #include <sys/stat.h>  //mkdir
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
