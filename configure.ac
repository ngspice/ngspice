# configure.ac script for ngspice
#
# This file is part of ngspice.
#
# (Process this file with autoconf to produce a configure script.)
#
# This configure file has been updated to mostly follow the guidence in the
# GNU Autoconf manual for version 2.59


# Initialisation
# ---------------

# Initialisation of configure
AC_INIT([ngspice], [26], [http://ngspice.sourceforge.net/bugrep.html])

# Unique file in the source directory
AC_CONFIG_SRCDIR([src/ngspice.c])

# Ensure a recent version of Autoconf is in use
# Older versions may not work with this script and this will report the
# problem to the user.  (2.52 is a guess and might need some adjustment)
AC_PREREQ([2.59])

# Revison stamp the generated ./configure script
AC_REVISION([$Revision$])

# Create a configuration header
AC_CONFIG_HEADER([src/include/ngspice/config.h])

# Initialise automake stuff
AM_INIT_AUTOMAKE([-Wall -Werror nostdinc])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

ext_CFLAGS="${CFLAGS+yes}"

AC_PROG_CC

m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

AC_SUBST([AM_CPPFLAGS], ['-I. -I$(srcdir) -I$(top_builddir)/src/include'])

AC_CONFIG_MACRO_DIR([m4])

# Package Options
# ---------------

# --enable-ansi: try to force --ansi option to the compiler
AC_ARG_ENABLE([ansi],
    [AS_HELP_STRING([--enable-ansi], [Force --ansi option for compilation])])

# --disable-debug: remove -g and -Wall  option to the compiler
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--disable-debug], [Remove -g option for compilation (Added by default)])])

# --enable-gprof: add -pg option to the compiler (GCC)
AC_ARG_ENABLE([gprof],
   [AS_HELP_STRING([--enable-gprof], [Produce gprof profiling data in 'gmon.out' (GCC only).])])

# --enable-checker: add --with-checker-debug  option to the compiler
AC_ARG_ENABLE([checker],
    [AS_HELP_STRING([--enable-checkergcc], [Option for compilation with checkergcc])])

# --enable-gc: Enable Boehm-Weiser Conservative Garbage Collector
AC_ARG_ENABLE([checker],
    [AS_HELP_STRING([--enable-gc], [Use Boehm-Weiser Conservative Garbage Collector. Disabled by default.])])

# --enable-nobypass: define NOBYPASS for the code
AC_ARG_ENABLE([nobypass],
    [AS_HELP_STRING([--enable-nobypass], [Don't bypass recalculations of slowly changing variables])])

# --enable-capbypass: define CAPBYPASS for the code
AC_ARG_ENABLE([capbypass],
    [AS_HELP_STRING([--enable-capbypass], [Bypass calculation of cbd/cbs in the mosfets if the vbs/vbd voltages are unchanged])])

# --enable-nodelimiting: define NODELIMITING for the code
AC_ARG_ENABLE([nodelimiting],
    [AS_HELP_STRING([--enable-nodelimiting], [Experimental damping scheme.])])

# --enable-predictor: define PREDICTOR for the code
AC_ARG_ENABLE([predictor],
    [AS_HELP_STRING([--enable-predictor], [Enable a predictor method for convergence.])])

# --enable-newpred: define NEWPRED for the code
AC_ARG_ENABLE([newpred],
    [AS_HELP_STRING([--enable-newpred], [Enable NEWPRED whatever it is(?)])])

# --enable-newtrunc: define NEWTRUNC for the code
AC_ARG_ENABLE([newtrunc],
    [AS_HELP_STRING([--enable-newtrunc], [Enable, how we want extrapolate capacitances.])])

# --enable-sense2: define WANT_SENSE2 for the code
AC_ARG_ENABLE([sense2],
    [AS_HELP_STRING([--enable-sense2], [Use spice2 sensitivity analysis.])])

# --enable-ftedebug: enable frontend debug macros
AC_ARG_ENABLE([ftedebug],
    [AS_HELP_STRING([--enable-ftedebug], [Enable ngspice frontend debug])])

# --enable-cpdebug: enable frontend debug macros
AC_ARG_ENABLE([cpdebug],
    [AS_HELP_STRING([--enable-cpdebug], [Enable ngspice shell debug])])

# --enable-sensdebug: define SENSDEBUG for the code
AC_ARG_ENABLE([sensdebug],
    [AS_HELP_STRING([--enable-sensdebug], [Debug sensitivity code *SENSDEBUG*])])

# --enable-asdebug: define ASDEBUG for the code
AC_ARG_ENABLE([asdebug],
    [AS_HELP_STRING([--enable-asdebug], [Debug sensitivity code *ASDEBUG*])])

# --enable-stepdebug: define STEPDEBUG for the code
AC_ARG_ENABLE([stepdebug],
    [AS_HELP_STRING([--enable-stepdebug], [Unknown debug option])])

# --enable-pzdebug: define PZDEBUG for the code
AC_ARG_ENABLE([pzdebug],
    [AS_HELP_STRING([--enable-pzdebug], [Debug pole/zero code])])

# --enable-blktmsdebug: define D_DBG_BLOCKTIMES for the code
AC_ARG_ENABLE([blktmsdebug],
    [AS_HELP_STRING([--enable-blktmsdebug], [Debug distortion code *BLOCKTIMES*])])

# --enable-smltmsdebug: define D_DBG_SMALLTIMES for the code
AC_ARG_ENABLE([smltmsdebug],
    [AS_HELP_STRING([--enable-smltmsdebug], [Debug distortion code *SMALLTIMES*])])

# --enable-smoketest: a smoketest
AC_ARG_ENABLE([smoketest],
    [AS_HELP_STRING([--enable-smoketest], [Enable smoketest compile])])

# --enable-experimental: define EXPERIMENTAL_CODE for the code
AC_ARG_ENABLE([experimental],
    [AS_HELP_STRING([--enable-experimental], [Enable some experimental code])])

# --enable-expdevices: Enable the compilation of experimental devices
AC_ARG_ENABLE([expdevices],
    [AS_HELP_STRING([--enable-expdevices], [Enable experimental devices (they do not compile)])])

# --enable-xspice: define XSPICE in the code. This is for xspice support
AC_ARG_ENABLE([xspice],
    [AS_HELP_STRING([--enable-xspice], [Enable XSpice enhancements, (experimental)])])

# --enable-cider: define CIDER in the code. This is for CIDER support
AC_ARG_ENABLE([cider],
    [AS_HELP_STRING([--enable-cider], [Enable CIDER enhancements, (experimental)])])

# --enable-adms: define ADMS in the code. This is for the adms Verilog-A compiler support
AC_ARG_ENABLE([adms],
    [AS_HELP_STRING([--enable-adms], [Enable ADMS code models, (experimental)])])

# --enable-pss: enable PSS Analysis
AC_ARG_ENABLE([pss],
    [AS_HELP_STRING([--enable-pss], [Enable PSS Analysis, (experimental)])])

# --enable-ndev: define NDEV in the code. An interface for external device i.e. numerical device
AC_ARG_ENABLE([ndev],
    [AS_HELP_STRING([--enable-ndev], [Enable NDEV interface, (experimental)])])

# --enable-cluster: define CLUSTER in the code. This is for cluster support
AC_ARG_ENABLE([cluster],
    [AS_HELP_STRING([--enable-cluster], [Enable cluster support, (experimental)])])

# --enable-help: try to force --ansi option to the compiler
AC_ARG_ENABLE([help],
    [AS_HELP_STRING([--enable-help], [Force building nghelp (deprecated)])])

# --enable-xgraph: Compile the xgraph plotting program. Default is "no".
# Xgraph is a plotting package for X11 once very popular.
AC_ARG_ENABLE([xgraph],
    [AS_HELP_STRING([--enable-xgraph], [Enable xgraph compilation.])])

# --with-readline:  Includes GNU readline support into CLI.  Default is "no".
AC_ARG_WITH([readline],
    [AS_HELP_STRING([--with-readline[=yes/no]], [Enable GNU readline support for CLI. Default=no.])])

# --with-editline:  Includes BSD Editline support into CLI.  Default is "no".
AC_ARG_WITH([editline],
    [AS_HELP_STRING([--with-editline[=yes/no]], [Enable BSD editline support for CLI. Default=no.])])

# --with-tcl: define TCL_MODULE in the code. This is for tcl support
AC_ARG_WITH([tcl],
    [AS_HELP_STRING([--with-tcl=tcldir], [Compiles  the tcl module instead,  experimental, see README.tcl])])

# --with-ngshared: define SHARED_MODULE in the code. This is for ngspice as shared library (or dll)
AC_ARG_WITH([ngshared],
    [AS_HELP_STRING([--with-ngshared], [Compiles ngspice as shared library (dll)])])

# readline and editline cannot both be enabled
if test "x$with_editline" = xyes; then
    if test "x$with_readline" = xyes; then
        AC_MSG_ERROR([Readline and editline cannot both be enabled])
    fi
fi


# Enable maintainer commands only if requested
AM_MAINTAINER_MODE


# Compiler checks
# ---------------

# Work on compiler options according to system:
# Set default CFLAGS - only use -Wall if we have gcc

# the above AC_PROG_CC may set CFLAGS to "-O2 -g"
if test "x$ext_CFLAGS" != xyes; then

    if test "x$enable_debug" = xno; then
        AC_MSG_WARN([Removing debugging option!])
        if test "x$GCC" = xyes; then
            CFLAGS="-O2 -s"
        else
            CFLAGS="-O2"
        fi
    else
        AC_DEFINE([NGDEBUG], [1], [Compile with debug info])
        if test "x$GCC" = xyes; then
            CFLAGS="-g -O1"
        else
            CFLAGS="-g"
        fi
    fi

    if test "x$GCC" = xyes; then
        CFLAGS="$CFLAGS -Wall -Wextra -Wmissing-prototypes -Wstrict-prototypes -Wnested-externs -Wold-style-definition -Wredundant-decls -Wconversion"
    fi

fi

# the profiling option
if test "x$enable_gprof" = xyes; then
    # Probably other compilers support -pg or something similar;
    # feel free to extend this to include them.
    if test "x$GCC" = xyes; then
        if test "x$enable_shared" = xyes ; then
            AC_MSG_ERROR([Can't have --enable-gprof without --disable-shared (we recommend also using --enable-all-static).])
        fi
        if test "x$enable_static" != xyes ; then
            AC_MSG_WARN([We recommend --enable-static with --enable-gprof.])
        fi
        AC_MSG_NOTICE([Enabling gprof profiling data (to gmon.out).])
        CFLAGS="$CFLAGS -pg"
        LDFLAGS="$LDFLAGS -pg"
    else
        AC_MSG_ERROR([We only support --enable-gprof with GCC right now.])
    fi
fi

# Not sure that this will work....
if test "x$with_checkergcc" = xyes; then
    CC="checkergcc"
fi


# Checks for ANSI-C header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([strings.h])

if test "x$enable_ansi" = xyes; then
    AC_PROG_CC()
    AC_DIAGNOSE([obsolete], [AM_PROG_CC_STDC:
	your code should no longer depend upon `am_cv_prog_cc_stdc', but upon
	`ac_cv_prog_cc_stdc'.  Remove this warning and the assignment when
	you adjust the code.  You can also remove the above call to
	AC_PROG_CC if you already called it elsewhere.])
    am_cv_prog_cc_stdc=$ac_cv_prog_cc_stdc

    if test "x$am_cv_prog_cc_stdc" != xyes ; then
        AC_MSG_WARN([Failed to find Ansi flag!])
    fi
fi

if test "x$enable_smoketest" = xyes; then
# CFLAGS="$CFLAGS -Werror"
    CFLAGS="$CFLAGS -pedantic -W -Wmissing-prototypes"
    CFLAGS="$CFLAGS -Wstrict-prototypes -Wtraditional"
    CFLAGS="$CFLAGS -Wconversion -Wshadow -Wpointer-arith"
    CFLAGS="$CFLAGS -Wcast-qual -Wcast-align -Wwrite-strings"
    CFLAGS="$CFLAGS -Waggregate-return -fshort-enums -fno-common"
    CFLAGS="$CFLAGS -Wnested-externs -Dinline= -g -O4"
fi


# Check system we're on , and tune accordingly
AC_CANONICAL_HOST

# Expand the prefix variable
if test "x$prefix" = xNONE ; then
#  with the MS WINDOWS executable we suggest C:/Spice/bin as *.exe-path
    if test "x$with_wingui" = xyes; then
        dprefix="C:/Spice"
        prefix=$dprefix
    else
        dprefix=$ac_default_prefix
    fi
else
    dprefix=$prefix
fi

##########################################################################
#
# ngspice as shared library (ngspice.so or ngspice.dll)
#
##########################################################################

if test "x$with_ngshared" != x  &&  test "x$with_ngshared" != xno ; then
    has_shared_module=true
    AC_DEFINE([SHARED_MODULE], [1], [shared ngspice module])
    case $host_os in
        *mingw* )
            has_shwin=true
            has_no_help=true
            AC_DEFINE([HAS_SHWIN], [1], [ngspice.dll Windows])
            ;;
        *cygwin* )
            has_shcyg=true
            AC_CHECK_LIB([pthread], [pthread_create])
            ;;
        *)
            CFLAGS="$CFLAGS -fvisibility=hidden"
            AC_CHECK_LIB([pthread], [pthread_create])
            ;;
    esac

    with_x=no
    with_windows=no

    AC_SUBST([STATIC], [-shared])

    TCL_PACKAGE_PATH=""
    TCL_LIB_SPEC=""


##########################################################################
#
# tcl libraries test
#
##########################################################################

#
# The tclSpice options
#
elif test "x$with_tcl" != x  &&  test "x$with_tcl" != xno ; then
    has_tcl_module=true
    AC_DEFINE([TCL_MODULE], [1], [Tcl Module])
    case $host_os in
        *mingw* )
            has_tclwin=true
            has_no_help=true
            AC_DEFINE([HAS_TCLWIN], [1], [Tcl Windows])
            ;;
        *cygwin* )
            has_tclcyg=true
            ;;
        *)
            ;;
    esac
 
    with_x=no
 
    AC_MSG_CHECKING([for tclConfig.sh])
    tcl_config_sh=""
    if test "x$with_tcl" != xyes ; then
        for dir in \
            $with_tcl
        do
            if test -r "$dir/tclConfig.sh" ; then
                tcl_config_sh="$dir/tclConfig.sh"
                break
            elif test -r "$dir/lib/tclConfig.sh" ; then
                tcl_config_sh="$dir/lib/tclConfig.sh"
                break
            elif test -r "$dir/unix/tclConfig.sh" ; then
                tcl_config_sh="$dir/unix/tclConfig.sh"
                break
            fi
        done
    else
        for dir in \
            $prefix \
            $exec_prefix
        do
            if test -r "$dir/tclConfig.sh" ; then
                tcl_config_sh="$dir/tclConfig.sh"
                break
            elif test -r "$dir/lib/tclConfig.sh" ; then
                tcl_config_sh="$dir/lib/tclConfig.sh"
                break
            elif test -r "$dir/unix/tclConfig.sh" ; then
                tcl_config_sh="$dir/unix/tclConfig.sh"
                break
            fi
        done
 
 
       if test "x$tcl_config_sh" = x ; then
           for dir in \
               `ls -dr /usr/local/tcltk/tcl[[7-9]].[[0-9]]* 2>/dev/null` \
               `ls -dr /usr/local/tcl/tcl[[7-9]].[[0-9]]* 2>/dev/null` \
               `ls -dr /usr/local/tcl[[7-9]].[[0-9]]* 2>/dev/null` \
               `ls -dr /usr/share/tcltk/tcl[[7-9]].[[0-9]]* 2>/dev/null` \
               `ls -dr /usr/share/tcl/tcl[[7-9]].[[0-9]]* 2>/dev/null` \
               `ls -dr /usr/share/tcl[[7-9]].[[0-9]]* 2>/dev/null` \
               `ls -dr /usr/lib/tcl/tcl[[7-9]].[[0-9]]* 2>/dev/null` \
               `ls -dr /usr/lib/tcl[[7-9]].[[0-9]]* 2>/dev/null` \
               /usr/local/tcl \
               /usr/local \
               /usr \
               /mingw
          do
              if test -r "$dir/tclConfig.sh" ; then
                  echo "Automatic research of tclConfig.sh Found one, accepting it and stopping research."
                  tcl_config_sh="$dir/tclConfig.sh"
                  break
              elif test -r "$dir/lib/tclConfig.sh" ; then
                  echo "Automatic research of tclConfig.sh Found one, accepting it and stopping research."
                  tcl_config_sh="$dir/lib/tclConfig.sh"
                  break
              fi
          done
       fi
    fi
 
    AC_MSG_RESULT([${tcl_config_sh}])
 
    if test "x$tcl_config_sh" = x ; then
        echo "can't find Tcl configuration script \"tclConfig.sh\""
        find /usr/ -name 'tclConfig.sh' -exec echo "Should you add --with-tcl={} to ./configure arguments?" \; | sed -re "s/\/tclConfig.sh/ /"
        exit 1
    fi
 
 .  $tcl_config_sh
 
    CFLAGS="$CFLAGS $TCL_INCLUDE_SPEC"
    CPPFLAGS="$CPPFLAGS $TCL_INCLUDE_SPEC"
 
    AC_CHECK_HEADERS([tcl.h blt.h], [], [AC_MSG_ERROR([Couldn't find Tcl/BLT headers])])
 
    if test ! -x "$TCL_EXEC_PREFIX/bin/tclsh$TCL_VERSION" ; then
        AC_MSG_ERROR([Couldn't find $TCL_EXEC_PREFIX/bin/tclsh$TCL_VERSION])
    fi

    AC_MSG_CHECKING([for TCL module BLT])

    rm -f conftest.tcl
    cat > conftest.tcl << EOF
package require BLT;
exit;
EOF

    if ($TCL_EXEC_PREFIX/bin/tclsh$TCL_VERSION conftest.tcl; exit) 2>/dev/null
    then
        AC_MSG_RESULT([Found])
    else
        AC_MSG_ERROR([Couldn't find BLT])
    fi
    rm -f conftest.tcl
    
    AC_CHECK_LIB([pthread], [pthread_create])
    
    if test -n "$libdir" ; then
        AC_SUBST([libdir])
        echo "libraries in $libdir"
    fi
    
    if test -z "$libdir" ; then
        for TCL_PACKAGE_PATH_ELEMENT in $TCL_PACKAGE_PATH ; do
            if test -a "$TCL_PACKAGE_PATH_ELEMENT" ; then
                libdir=$TCL_PACKAGE_PATH_ELEMENT
                AC_SUBST([libdir])
                echo "library in $libdir"
                break
            fi
        done
    fi

    AC_SUBST([STATIC], [-shared])

else
    TCL_PACKAGE_PATH=""
    TCL_LIB_SPEC=""
    AC_SUBST([STATIC], [-static])
fi


AM_CONDITIONAL([SHARED_MODULE], [test "x$has_shared_module" = xtrue])
AM_CONDITIONAL([SHWIN], [test "x$has_shwin" = xtrue])
AM_CONDITIONAL([SHCYG], [test "x$has_shcyg" = xtrue])
AM_CONDITIONAL([TCL_MODULE], [test "x$has_tcl_module" = xtrue])
AM_CONDITIONAL([TCLWIN], [test "x$has_tclwin" = xtrue])
AM_CONDITIONAL([TCLCYG], [test "x$has_tclcyg" = xtrue])

AC_SUBST([TCL_PACKAGE_PATH])
AC_SUBST([TCL_LIB_SPEC])
##################################################################
#
#
#End of tcl libraries test
#
#################################################################



# Checks for programs
AC_LIBTOOL_WIN32_DLL
AC_LIBTOOL_DLOPEN
AM_PROG_CC_C_O
## AM_PROG_LIBTOOL
AC_PROG_LIBTOOL


# --with-windows : the user wants to use generate the MS WINDOWS GUI executable
# old command version, deprecated, use --with-wingui
AC_ARG_WITH([windows],
    [AS_HELP_STRING([--with-windows], [MS WINDOWS GUI executable, deprecated])],
    [if test "x$with_windows" = xyes; then
        AC_DEFINE([HAS_WINGUI], [1], [define to specify MS Windows executable with simple GUI])
    fi],
    [with_windows=no])

# --with-wingui : the user wants to use generate the MS WINDOWS GUI executable
AC_MSG_CHECKING([whether windows code with gui is enabled])
AC_ARG_WITH([wingui],
    [AS_HELP_STRING([--with-wingui], [MS WINDOWS GUI executable])],
    [if test "x$with_wingui" = xyes; then
        AC_DEFINE([HAS_WINGUI], [1], [define to specify MS Windows executable with simple GUI])
    fi],
    [if test "x$with_windows" = xyes; then
        with_wingui=yes
    else
        with_wingui=no
    fi])
AC_MSG_RESULT([$with_wingui])
AM_CONDITIONAL([WINGUI], [test "x$with_wingui" = xyes])

case $with_wingui in
    yes )
        AC_DEFINE([X_DISPLAY_MISSING])
        AC_MSG_RESULT([No X display!])
        has_no_x=true
        if test "x$enable_help" = xyes; then
            has_no_help=false
        else
            has_no_help=true
        fi
        CFLAGS="$CFLAGS -mwindows"
        ;;
    * )

# Check for /proc (virtual process information file system)
AC_CHECK_HEADERS([/proc/meminfo])

# Checks for X11 header files and libraries - X11 support can be disabled
# by passing the '--without-x' option to configure:

# Try to locate the X-Windows include files and libraries
#   will set no_x="yes" or no_x=""
AC_PATH_X

# will set X_CFLAGS, X_LIBS, and might define X_DISPLAY_MISSING
AC_PATH_XTRA

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ALLOCA
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_CLOSEDIR_VOID
AC_FUNC_SELECT_ARGTYPES
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_VPRINTF

#  Check if the user wants to use the appropriate rpath commands to compile in
#  the shared library path for the X libraries.  This is off by default since
#  that is how it has been historically.  Using --enable-rpath will turn this on

_use_rpath=no
AC_MSG_CHECKING([if hardcoding of the X11 runtime library path is desired])
AC_ARG_ENABLE([rpath],
    [AS_HELP_STRING([--enable-rpath], [Enable hardcoding the X11 runtime library path])],
    [   if test "x$enable_rpath" = xno ; then
            AC_MSG_RESULT([no])
            _use_rpath=no
        else
            AC_MSG_RESULT([yes])
            _use_rpath=yes
        fi
    ],
    [   AC_MSG_RESULT([no])
        _use_rpath=no
    ]
)

if test "x$_use_rpath" = xyes ; then
    # Try to figure out if we need -Rpath for finding X11 libs
    # at runtime.  Why autoconf doesn't already do this, I don't
    # know...
    xlib_path=""
    for p in $X_LDFLAGS ; do
        case $p in
            -L*)
                xlib_path="$xlib_path $p"
                ;;
        esac
    done

    _save_LIBS=$LIBS
    LIBS="$LIBS $X_LDFLAGS"
    rpath=""
    for fl in "-Wl,-rpath " "-Wl,--rpath " "-rpath " "--rpath " "-Wl,-R" "-R" "-R " ; do
        xlib_rpath=`echo $xlib_path | sed "s/-L/$fl/g"`
        LIBS="$_save_LIBS $X_LIBS $xlib_rpath"
        AC_MSG_CHECKING([if the compiler accepts ${fl}path for runtime libraries])
        AC_LINK_IFELSE([AC_LANG_PROGRAM()],
            [AC_MSG_RESULT([yes])
            rpath=$fl],
            [AC_MSG_RESULT([no])])
        test -n "$rpath" && break
    done
    if test -n "$rpath"; then
        X_LDFLAGS="$X_LDFLAGS $xlib_rpath"
    fi
    LIBS=$_save_LIBS
fi

# Checks for X libraries - if X11 wasn't found then don't make following
# tests and compile without X11 support - otherwise, check if the following
# libraries are present (error if they are not)
# In CYGWIN library ordering has to be changed. Is this compatible to LINUX?
# XShmAttach is a struct in CYGWIN, not a function
# SJB: 13th march 2005
# Library order is giving linker warnings on MacOSX
# It's not clear to me which order is required for Cygwin (see comment above)
# and neither can I find any authoritative answer for the correct link order
# for MacOSX or Linux, but
#            -lXaw -lXmu -lXt -lXext -lX11
# seems to be the popular choice.
# (The previous order was -lX11 -lXt -lXext -lXmu -lXaw)

if test "x$no_x" != xyes ; then
    AC_CHECK_LIB([Xaw], [main], [X_LIBS="$X_LIBS -lXaw"],
        [AC_MSG_ERROR([Couldn't find Xaw library])], [$X_LIBS $X_EXTRA_LIBS])
    AC_CHECK_LIB([Xmu], [main], [X_LIBS="$X_LIBS -lXmu"],
        [AC_MSG_ERROR([Couldn't find Xmu library])], [$X_LIBS $X_EXTRA_LIBS])
    X_LIBS="$X_LIBS -lXt"
    AC_CHECK_LIB([Xext], [XShmAttach], [X_LIBS="$X_LIBS -lXext"],
        [AC_MSG_ERROR([Couldn't find Xext library])], [$X_LIBS $X_EXTRA_LIBS])
    X_LIBS="$X_LIBS -lX11"
    has_no_x=false
else
    AC_DEFINE([X_DISPLAY_MISSING])
    AC_MSG_RESULT([No X display!])
    has_no_x=true
fi

if test "x$enable_help" = xyes  &&  test "x$no_x" != xyes ; then
    has_no_help=false
else
    has_no_help=true
fi ;;

esac
# end of case --with-windows

AM_CONDITIONAL([NO_X], [test "x$has_no_x" = xtrue])
AM_CONDITIONAL([NO_HELP], [test "x$has_no_help" = xtrue])

# Check for a few typdefs:
AC_TYPE_PID_T
AC_TYPE_SIGNAL

# sighandler_t is the type of the signal handler on GNU variants,
# sig_t is the type of a signal handler on 4.4BSD's,
# other systems use __sighandler_t.
AC_CHECK_TYPES([sighandler_t, sig_t, __sighandler_t], [], [],
    [#include <sys/types.h>
     #include <signal.h>])

# Check for a few libraries and headers:
AC_HEADER_DIRENT
AC_CHECK_HEADERS([unistd.h ctype.h pwd.h fcntl.h sys/ioctl.h stropts.h])
AC_HEADER_SYS_WAIT
AC_HEADER_STAT

AC_CHECK_HEADERS([arpa/inet.h netdb.h netinet/in.h stddef.h sys/file.h sys/param.h sys/socket.h sys/time.h sys/timeb.h sys/io.h])


# Check time and resources headers and functions:
AC_HEADER_TIME
AC_STRUCT_TM
AC_STRUCT_TIMEZONE
AC_CHECK_FUNCS([localtime])

AC_CHECK_FUNCS([ftime gettimeofday])
# Do not use time or getrusage function for CPU time measurement under OpenMP
if test "x$enable_openmp" != xyes; then
    AC_CHECK_FUNCS([time getrusage])
fi
AC_CHECK_FUNCS([utimes])
AC_CHECK_FUNCS([getrlimit ulimit], [break])

AC_CHECK_FUNCS([endpwent gethostbyname memset select socket strdup strerror strncasecmp strstr strtol])

# Look for termios first (posix)
AC_CHECK_HEADERS([termios.h termio.h sgtty.h], [break])
AC_CHECK_FUNCS([isatty tcgetattr tcsetattr])

# Check for a few functions:
AC_FUNC_FORK([])
AC_CHECK_FUNCS([access bcopy bzero qsort dup2 popen])
AC_CHECK_FUNCS([strchr index], [break])
AC_CHECK_FUNCS([strrchr rindex], [break])
AC_CHECK_FUNCS([getcwd getwd], [break])


AC_MSG_RESULT([Checking mathematical features of the system:])
# Look for math library:
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_HEADERS([float.h limits.h values.h ieeefp.h])
AC_CHECK_HEADERS([fftw3.h])
AC_CHECK_LIB([fftw3], [fftw_plan_dft_1d],
        [AC_DEFINE([HAVE_LIBFFTW3], [], [Have fft routines in libfftw3])
         LIBS="$LIBS -lfftw3"])

# Check for a few mathematical functions:
AC_CHECK_FUNCS([erfc logb scalb scalbn asinh acosh atanh finite])
# According POSIX we should look for macros first
AC_CHECK_DECLS([isinf], [], [], [[#include <math.h>]])
AC_CHECK_DECLS([isnan], [], [], [[#include <math.h>]])
# But may be we have still functions
if test "x$ac_cv_have_decl_isinf" != xyes; then
    AC_CHECK_FUNC([isinf])
fi
if test "x$ac_cv_have_decl_isnan" != xyes; then
    AC_CHECK_FUNC([isnan])
fi

# If user enables garbage collection, look for garbage collector
if test "x$TCL_PACKAGE_PATH" = x; then
    if test "x$enable_gc" = xyes; then
        AC_MSG_RESULT([Checking for the presence of the Garbage Collector:])
        AC_CHECK_LIB([gc], [GC_malloc],
            [AC_DEFINE([HAVE_LIBGC], [], [Define if we want garbage collection enabled])
             LIBS="$LIBS -lgc"])
    fi
fi
# check, if we have sigsetjmp and siglongjmp.
# A trivial AC_CHECK_FUNCS(sigsetjmp) won't do because sigsetjmp() might be a
# macro declared in <setjmp.h>. (joze)
AC_MSG_CHECKING([for sigsetjmp])
AC_TRY_LINK([#include <setjmp.h>],[jmp_buf env; sigsetjmp(env, 1);],
    [AC_MSG_RESULT([yes])
        AC_DEFINE([HAVE_SIGSETJMP], [], [Define if we have sigsetjmp().])],
    [AC_MSG_RESULT([no])])


# Check for the snprintf function:
AC_CHECK_FUNCS([snprintf])

# Check for the dirname function:
AC_CHECK_FUNCS([dirname], [],
    [AC_CHECK_LIB([gen], [dirname],
        [AC_DEFINE([HAVE_DIRNAME], [1], [Have dirname in libgen])
         LIBS="$LIBS -lgen"])])

AC_CHECK_HEADERS([getopt.h])
AC_CHECK_FUNC([getopt_long], [getopt_long=true])
if test "x$getopt_long" = xtrue; then
    AC_DEFINE([HAVE_GETOPT_LONG], [1], [Have fcn getopt_long()])
fi

AC_DEFINE_UNQUOTED([NGSPICEBINDIR], ["`echo $dprefix/bin`"], [Define the directory for executables])
AC_DEFINE_UNQUOTED([NGSPICEDATADIR], ["`echo $dprefix/share/ngspice`"], [Define the directory for architecture independent data files])
AC_DEFINE_UNQUOTED([NGSPICEBUILDDATE], ["`date`"], [Define the build date])

if test "x$with_wingui" = xyes; then
    AC_MSG_RESULT([WINDOWS GUI code enabled])
    AC_CHECK_FUNCS([memmove])
fi

# Recapitulate settings:
AC_MSG_RESULT([Settings which were chosen:])

if test "x$enable_sense2" = xyes; then
    AC_DEFINE([WANT_SENSE2], [], [Define if we want spice2 sensitivity analysis])
    AC_MSG_RESULT([Spice2 sensitivity analysis enabled])
    AC_MSG_WARN([This feature is UNSUPPORTED])
fi
if test "x$enable_nobypass" = xyes; then
    AC_DEFINE([NOBYPASS], [], [Define if we want NOBYPASS])
    AC_MSG_RESULT([NOBYPASS option enabled])
fi
if test "x$enable_capbypass" = xyes; then
    AC_DEFINE([CAPBYPASS], [], [Define if we want to bypass cbd/cbs calculation for non varying vbs/vbd])
    AC_MSG_RESULT([CAPBYPASS option enabled])
fi
if test "x$enable_nodelimiting" = xyes; then
    AC_DEFINE([NODELIMITING], [], [Experimental code never implemented to damp Newton iterations])
    AC_MSG_RESULT([NODELIMITING option enabled])
fi
if test "x$enable_predictor" = xyes; then
    AC_DEFINE([PREDICTOR], [], [Define if we want predictor algorithm])
    AC_MSG_RESULT([PREDICTOR algorithm enabled])
fi
if test "x$enable_newpred" = xyes; then
    AC_DEFINE([NEWPRED], [], [Define if you want to discover :)])
    AC_MSG_RESULT([NEWPRED enabled])
fi
if test "x$enable_newtrunc" = xyes; then
    AC_DEFINE([NEWTRUNC], [], [Do not trigger unwanted traps by default])
    AC_MSG_RESULT([New truncation error calculation enabled])
fi
if test "x$enable_experimental" = xyes; then
    AC_DEFINE([EXPERIMENTAL_CODE], [], [Define if we want some experimental code])
    AC_MSG_RESULT([EXPERIMENTAL_CODE enabled])
fi
if test "x$enable_cpdebug" = xyes; then
    AC_DEFINE([CPDEBUG], [], [Define if you want to debug ngspice shell])
    AC_MSG_RESULT([WARNING: Shell debug is enabled])
fi
if test "x$enable_ftedebug" = xyes; then
    AC_DEFINE([FTEDEBUG], [], [Define if you want to debug frontend])
    AC_MSG_RESULT([WARNING: Frontend debug is enabled])
fi
if test "x$enable_sensdebug" = xyes; then
    AC_DEFINE([SENSDEBUG], [], [Define if we want debug sensititvity analysis])
    AC_MSG_RESULT([WARNING: Sensitivity code debug *SENSDEBUG* is enabled])
fi
if test "x$enable_asdebug" = xyes; then
    AC_DEFINE([ASDEBUG], [], [Define if we want debug sensitivity analysis])
    AC_MSG_RESULT([WARNING: Sensitivity code debug *ASDEBUG* is enabled])
fi
if test "x$enable_stepdebug" = xyes; then
    AC_DEFINE([STEPDEBUG], [], [Define if we want stepdebug])
    AC_MSG_RESULT([WARNING: STEPDEBUG debug is enabled])
fi
if test "x$enable_pzdebug" = xyes; then
    AC_DEFINE([PZDEBUG], [], [Define if you want to debug pole-zero analysis])
    AC_MSG_RESULT([WARNING: Pole/Zero analysis debug is enabled])
fi
if test "x$enable_pss" = xyes; then
    AC_DEFINE([WITH_PSS], [], [Define if you want PSS analysis])
    AC_MSG_RESULT([WARNING: PSS analysis enabled])
fi
if test "x$enable_blktmsdebug" = xyes; then
    AC_DEFINE([D_DBG_BLOCKTIMES], [], [Define if we want debug distortion analysis (BLOCKTIMES)])
    AC_MSG_RESULT([WARNING: Distortion analysis debug *D_DBG_BLOCKTIMES* is enabled])
fi
if test "x$enable_smltmsdebug" = xyes; then
    AC_DEFINE([D_DBG_SMALLTIMES], [], [Define if we want debug distortion analysis (SMALLTIMES)])
    AC_MSG_RESULT([WARNING: Distortion analysis debug *D_DBG_SMALLTIMES* is enabled])
fi
if test "x$enable_xgraph" = xyes; then
    AC_MSG_RESULT([Xgraph compilation enabled.])
    AC_CONFIG_SUBDIRS([xgraph])
    XGRAPHDIR="xgraph"
    NOTXGRAPH=""
else
    XGRAPHDIR=""
    NOTXGRAPH="xgraph"
fi

AC_SUBST([XGRAPHDIR])
AC_SUBST([NOTXGRAPH])

AC_CHECK_PROGS([YACC], ['bison -y' byacc yacc])

#################  XSPICE ##################################################
# Add new code models to the build by pointing to them here.

if test "x$enable_xspice" = xyes; then
    if test -z "${YACC}" ; then
        AC_MSG_ERROR([XSPICE build requires bison, byacc or yacc parser generator])
    fi
    AC_MSG_RESULT([X-Spice features included])
    AC_DEFINE([XSPICE], [1], [The xspice enhancements])
# Define variables for LEX
    AC_CHECK_PROGS([LEX], [flex lex])
    if test -z "$LEX" ; then
        AC_MSG_ERROR([Flex is required for building XSPICE])
    fi
    case $host_os in
        *mingw* )
            AC_DEFINE([IPC_DEBUG_VIA_STDIO], [1], [Client-Server only via stdio.])
            XSPICEDLLIBS=""
            ;;
        *cygwin* )
            AC_DEFINE([IPC_UNIX_SOCKETS], [1], [Client-Server via socket.])
            XSPICEDLLIBS=""
            ;;
        *freebsd* )
            AC_DEFINE([IPC_UNIX_SOCKETS], [1], [Client-Server via socket.])
            XSPICEDLLIBS=""
            ;;
        *openbsd* )
            AC_DEFINE([IPC_UNIX_SOCKETS], [1], [Client-Server via socket.])
            XSPICEDLLIBS=""
            ;;
        * )
            AC_DEFINE([IPC_UNIX_SOCKETS], [1], [Client-Server via socket.])
            XSPICEDLLIBS="-ldl"
            ;;
    esac
    XSPICEINIT=""
    AC_CHECK_FUNCS([modf])
    AC_CHECK_HEADERS([libintl.h malloc.h])
else
    XSPICEINIT="*"

    if test -z "${YACC}" ; then
        AC_MSG_WARN([No bison, byacc, yacc found: Uses prebuilt default parsers])
    fi
fi
AC_SUBST([XSPICEINIT])
AC_SUBST([XSPICEDLLIBS])

AM_CONDITIONAL([XSPICE_WANTED], [test "x$enable_xspice" = xyes])


# Add CIDER enhancements to ngspice.
if test "x$enable_cider" = xyes; then
    AC_MSG_RESULT([CIDER features enabled])
    AC_DEFINE([CIDER], [1], [The CIDER enhancements])
fi

AM_CONDITIONAL([CIDER_WANTED], [test "x$enable_cider" = xyes])
AM_CONDITIONAL([NUMDEV_WANTED], [test "x$enable_cider" = xyes])


AM_CONDITIONAL([PSS_WANTED], [test "x$enable_pss" = xyes])

AM_CONDITIONAL([SENSE2_WANTED], [test "x$enable_sense2" = xyes])

# adms option
if test "x$enable_adms" = xyes ; then
    AC_MSG_RESULT([**********************************
* ADMS  support is  experimental *
**********************************])
    AC_CHECK_PROGS([ADMSXML], [admsXml admsXml.exe], [no])

    if test "x$ADMSXML" = xno; then
        AC_MSG_ERROR([If you want Verilog-A models you should install admsXml])
    fi
    AC_DEFINE([ADMS], [1], [Support for Verilog-A models])

    VLADEVDIR=" adms/ekv \
                adms/hicum0 \
                adms/hicum2 \
                adms/mextram \
                adms/psp102 "

# The makefiles for adms (to be added to AC_CONFIG_FILES by ./autogen.sh --adms)
#VLAMKF    src/spicelib/devices/adms/ekv/Makefile
#VLAMKF    src/spicelib/devices/adms/hicum0/Makefile
#VLAMKF    src/spicelib/devices/adms/hicum2/Makefile
#VLAMKF    src/spicelib/devices/adms/mextram/Makefile
#VLAMKF    src/spicelib/devices/adms/psp102/Makefile

    NOTVLADEVDIR=""

    VLADEV=" spicelib/devices/adms/ekv/libekv.la \
             spicelib/devices/adms/hicum0/libhicum0.la \
             spicelib/devices/adms/hicum2/libhicum2.la \
             spicelib/devices/adms/mextram/libbjt504t.la \
             spicelib/devices/adms/psp102/libpsp102.la "

else

    VLADEVDIR=""
    NOTVLADEVDIR="adms"

fi

AC_SUBST([ADMSXML])
AC_SUBST([VLADEVDIR])
AC_SUBST([VLADEV])
AC_SUBST([NOTVLADEVDIR])

# NDEV option
if test "x$enable_ndev" = xyes; then
    AC_MSG_RESULT([NDEV features enabled])
    AC_DEFINE([NDEV], [], [The NDEV interface])
fi
AM_CONDITIONAL([NDEV_WANTED], [test "x$enable_ndev" = xyes])


# Cluster option
if test "x$enable_cluster" = xyes; then
    AC_MSG_RESULT([Cluster version is being compiled])
    AC_DEFINE([CLUSTER], [], [Spice cluster support])
    LIBS="$LIBS -lpthread"
fi

if test "x$enable_expdevices" = xyes; then
    AC_DEFINE([EXP_DEV], [], [Define if we want to enable experimental devices])
    AC_MSG_RESULT([WARNING: Experimental devices enabled])
fi

#  ----  Option to include GNU readline support in ngspice CLI  ----
#  ----  Default: disabled.                                     ----
#  ----  Hope to see in the future readline replacement.        ----

if test "x$with_readline" != xyes; then
    AC_MSG_RESULT([GNU readline disabled.])
else
    if test "x$with_tcl" = x  ||  test "x$with_tcl" = xno ; then
        AC_MSG_RESULT([Checking for readline:])
        AC_CHECK_HEADERS([readline/readline.h readline/history.h],
            [AC_DEFINE([HAVE_GNUREADLINE], [], [Define if we have GNU readline])],
            [AC_MSG_ERROR([Couldn't find GNU readline headers.])])
        AC_SEARCH_LIBS([tputs], [ncurses termcap],
            [AC_DEFINE([HAVE_TERMCAP], [], [Define if we have ncurses or termcap])],
            [AC_MSG_ERROR([Found neither ncurses or termcap])])
        AC_CHECK_LIB([readline], [readline],
            [LIBS="$LIBS -lreadline"],
            [AC_MSG_ERROR([Couldn't find readline libraries.])])
    fi
fi


#  ----  Option to include BSD editline support in ngspice CLI  ----
#  ----  Default: disabled.                                     ----

if test "x$with_editline" != xyes; then
    AC_MSG_RESULT([BSD editline disabled.])
else
    AC_MSG_RESULT([Checking for editline:])
    AC_CHECK_HEADERS([editline/readline.h],
        [AC_DEFINE([HAVE_BSDEDITLINE], [1], [Define to enable BSD editline])],
        [AC_MSG_ERROR([Couldn't find BSD editline headers.])])
    AC_SEARCH_LIBS([tputs], [ncurses termcap],
        [AC_DEFINE([HAVE_TERMCAP], [], [Define if we have ncurses or termcap])],
        [AC_MSG_ERROR([Found neither ncurses or termcap])])
    AC_CHECK_LIB([edit], [readline],
        [LIBS="$LIBS -ledit"],
        [AC_MSG_ERROR([Couldn't find editline libraries.])],
        [-lncurses])
fi

# Use AC_CHECK_HEADERS so the HAVE_*_H symbol gets defined
AC_CHECK_HEADERS([ncurses/termcap.h termcap.h])

# --enable-openmp: Use OpenMP on multi-core processors
AC_ARG_ENABLE([openmp],
    [AS_HELP_STRING([--enable-openmp], [Use OpenMP parallel processing])])

# Add OpenMP to ngspice.
: ${enable_openmp:=no}
m4_ifdef([AC_OPENMP], [AC_OPENMP])

if test "x$enable_openmp" = xyes; then
    AC_DEFINE([USE_OMP], [1], [OpenMP parallel processing])
    CFLAGS="$CFLAGS $OPENMP_CFLAGS"
    AC_MSG_RESULT([OpenMP feature enabled])
fi

# Output Files
# ------------

AC_CONFIG_FILES([Makefile
                 man/Makefile
                 man/man1/Makefile
                 src/Makefile
                 src/spicelib/Makefile
                 src/spicelib/analysis/Makefile
                 src/spicelib/devices/Makefile
                 src/spicelib/devices/asrc/Makefile
                 src/spicelib/devices/bjt/Makefile
                 src/spicelib/devices/bsim1/Makefile
                 src/spicelib/devices/bsim2/Makefile
                 src/spicelib/devices/bsim3/Makefile
                 src/spicelib/devices/bsim3v0/Makefile
                 src/spicelib/devices/bsim3v1/Makefile
                 src/spicelib/devices/bsim3v32/Makefile
                 src/spicelib/devices/bsim4/Makefile
                 src/spicelib/devices/bsim4v5/Makefile
                 src/spicelib/devices/bsim4v6/Makefile
                 src/spicelib/devices/bsim4v7/Makefile
                 src/spicelib/devices/bsim3soi_pd/Makefile
                 src/spicelib/devices/bsim3soi_fd/Makefile
                 src/spicelib/devices/bsim3soi_dd/Makefile
                 src/spicelib/devices/bsimsoi/Makefile
                 src/spicelib/devices/cap/Makefile
                 src/spicelib/devices/cccs/Makefile
                 src/spicelib/devices/ccvs/Makefile
                 src/spicelib/devices/csw/Makefile
                 src/spicelib/devices/cpl/Makefile
                 src/spicelib/devices/dio/Makefile
                 src/spicelib/devices/ind/Makefile
                 src/spicelib/devices/isrc/Makefile
                 src/spicelib/devices/hfet1/Makefile
                 src/spicelib/devices/hfet2/Makefile
                 src/spicelib/devices/hisim2/Makefile
                 src/spicelib/devices/hisimhv1/Makefile
                 src/spicelib/devices/jfet/Makefile
                 src/spicelib/devices/jfet2/Makefile
                 src/spicelib/devices/ltra/Makefile
                 src/spicelib/devices/mes/Makefile
                 src/spicelib/devices/mesa/Makefile
                 src/spicelib/devices/mos1/Makefile
                 src/spicelib/devices/mos2/Makefile
                 src/spicelib/devices/mos3/Makefile
                 src/spicelib/devices/mos6/Makefile
                 src/spicelib/devices/mos9/Makefile
                 src/spicelib/devices/ndev/Makefile
                 src/spicelib/devices/res/Makefile
                 src/spicelib/devices/soi3/Makefile
                 src/spicelib/devices/sw/Makefile
                 src/spicelib/devices/tra/Makefile
                 src/spicelib/devices/txl/Makefile
                 src/spicelib/devices/urc/Makefile
                 src/spicelib/devices/vbic/Makefile
                 src/spicelib/devices/vccs/Makefile
                 src/spicelib/devices/vcvs/Makefile
                 src/spicelib/devices/vsrc/Makefile
                 src/spicelib/devices/nbjt/Makefile
                 src/spicelib/devices/nbjt2/Makefile
                 src/spicelib/devices/numd/Makefile
                 src/spicelib/devices/numd2/Makefile
                 src/spicelib/devices/numos/Makefile
                 src/spicelib/parser/Makefile
                 src/ciderlib/Makefile
                 src/ciderlib/input/Makefile
                 src/ciderlib/support/Makefile
                 src/ciderlib/oned/Makefile
                 src/ciderlib/twod/Makefile
                 src/frontend/Makefile
                 src/frontend/numparam/Makefile
                 src/frontend/help/Makefile
                 src/frontend/parser/Makefile
                 src/frontend/plotting/Makefile
                 src/frontend/trannoise/Makefile
                 src/frontend/wdisp/Makefile
                 src/include/ngspice/Makefile
                 src/maths/Makefile
                 src/maths/cmaths/Makefile
                 src/maths/fft/Makefile
                 src/maths/misc/Makefile
                 src/maths/ni/Makefile
                 src/maths/deriv/Makefile
                 src/maths/poly/Makefile
                 src/maths/sparse/Makefile
                 src/misc/Makefile
                 src/xspice/Makefile
                 src/xspice/cm/Makefile
                 src/xspice/cmpp/Makefile
                 src/xspice/icm/makedefs
                 src/xspice/icm/GNUmakefile
                 src/xspice/mif/Makefile
                 src/xspice/evt/Makefile
                 src/xspice/enh/Makefile
                 src/xspice/ipc/Makefile
                 src/xspice/idn/Makefile
                 src/unsupported/Makefile
                 tests/Makefile
                 tests/bsim1/Makefile
                 tests/bsim2/Makefile
                 tests/bsim3/Makefile
                 tests/bsim3soidd/Makefile
                 tests/bsim3soifd/Makefile
                 tests/bsim3soipd/Makefile
                 tests/bsim4/Makefile
                 tests/bsimsoi/Makefile
                 tests/filters/Makefile
                 tests/general/Makefile
                 tests/hfet/Makefile
                 tests/hisim/Makefile
                 tests/hisimhv1/Makefile
                 tests/jfet/Makefile
                 tests/mes/Makefile
                 tests/mesa/Makefile
                 tests/mos6/Makefile
                 tests/polezero/Makefile
                 tests/regression/Makefile
                 tests/regression/subckt-processing/Makefile
                 tests/regression/lib-processing/Makefile
                 tests/regression/parser/Makefile
                 tests/sensitivity/Makefile
                 tests/transient/Makefile
                 tests/transmission/Makefile
                 tests/xspice/Makefile
                 tests/xspice/digital/Makefile
                 tests/xspice/digital/spinit
                 tests/resistance/Makefile
                 tests/vbic/Makefile])
AC_OUTPUT
