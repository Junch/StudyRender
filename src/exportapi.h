#pragma once

// Define EXPORTAPI for any platform
#ifdef _WIN32
#if defined (_MSC_VER)
#pragma warning(disable: 4251)
#endif
# ifdef WIN_EXPORT
#   define EXPORTAPI  __declspec(dllexport)
# else
#   define EXPORTAPI  __declspec(dllimport)
# endif
#else
# define EXPORTAPI
#endif
