AC_DEFUN([AC_CHECK_SCILAB],
[
  AH_TEMPLATE(HAVE_SCILAB,
   [Defined to 1 if Scilab is present on the system])

  AC_ARG_ENABLE(scilab,
  [AS_HELP_STRING([--enable-scilab],[turn on/off scilab support])],
  [case "${enableval}" in
  	 yes) usescilab=YES ;;
   	 no)  usescilab=NO ;;
   	 *) AC_MSG_ERROR([bad value ${enableval} for --enable-scilab]) ;;
   esac],[usescilab=NO])

  AC_ARG_WITH(scilab_prefix,
		AC_HELP_STRING([--with-scilab-prefix=DIR],[Set the path to Scilab]),
		[with_scilab_prefix=$withval],
		[with_scilab_prefix='yes']
		)

  AC_ARG_WITH(scilab_version,
		AC_HELP_STRING([--with-scilab-version="major.minor.micro"],[Set the required Scilab version]),
		[with_scilab_version=$withval],
		[with_scilab_version='yes']
		)

  AC_ARG_WITH(scilab_toolbox_dir,
		AC_HELP_STRING([--with-scilab-toolbox-dir=DIR],[Set the path to the toolbox installation directory]),
		[with_scilab_toolbox_dir=$withval],
		[with_scilab_toolbox_dir='yes']
		)

  if test "x$with_scilab_version" != "xyes"
  then
    REQUIRED_SCILAB_MAJOR=`expr match "$with_scilab_version" " *\([0-9]*\)"`
    REQUIRED_SCILAB_MINOR=`expr match "$with_scilab_version" " *[0-9]*.\([0-9]*\)"
    REQUIRED_SCILAB_MICRO=`expr match "$with_scilab_version" " *[0-9]*.[0-9]*.\([0-9]*\)"
  else
    REQUIRED_SCILAB_MAJOR=5
    REQUIRED_SCILAB_MINOR=1
    REQUIRED_SCILAB_MICRO=1
  fi

  dnl check for Scilab

  if test "x$with_scilab_prefix" != "xyes"
  then
    if test -x "$with_scilab_prefix/bin/scilab"
    then 
      AC_MSG_RESULT([Scilab binary program was found in $with_scilab_prefix])
    else
      AC_MSG_ERROR([Scilab binary program was not found in $with_scilab_prefix/bin])
    fi
    SCILAB_EXE="$with_scilab_prefix/bin/scilab"
    AC_DEFINE(HAVE_SCILAB)
  else
    AC_CHECK_PROG([has_scilab],[scilab],[yes],[no], , )
    if test x$has_scilab = xno; then
      AC_MSG_ERROR([[Scilab binary program was found in your PATH], your PATH is $PATH])
    fi
    SCILAB_EXE="scilab"
    AC_DEFINE(HAVE_SCILAB)
  fi

  dnl test if SCI is defined. If not, test if scilab is in path,
  dnl then get SCI from scilab -nwni.
  if test -z "$SCI"; then
    cmd='F=mopen("getpath.incl","w");
         mfprintf(F,SCI);
         mclose(F);exit;'
    echo "$cmd" > getpath.sci
    $SCILAB_EXE -nw -f getpath.sci >/dev/null
    SCILAB_DIR=`cat getpath.incl`
    rm -f getpath.sci getpath.incl 
  else
    SCILAB_DIR="$SCI"
  fi

  dnl get scilab version
  cmd='F=mopen("version.incl","w");
       ver=getversion();
       mfprintf(F,ver);
       mclose(F);exit;'
  echo "$cmd" > version.sci
  $SCILAB_EXE -nwni -f version.sci >/dev/null
  SCILAB_VERSION=`cat version.incl`
  rm -f version.sci version.incl 

  scilab_tmp_version=`expr match $SCILAB_VERSION '.*\(branch\).*'`

  if test "x$scilab_tmp_version" = "xbranch"
  then
    SCILAB_VERSION_MAJOR=-1
    SCILAB_VERSION_MINOR=-1
    SCILAB_VERSION_MICRO=-1
  else
    SCILAB_VERSION_MAJOR=`expr match "$SCILAB_VERSION" " *\([0-9]*\)"`
    SCILAB_VERSION_MINOR=`expr match "$SCILAB_VERSION" " *[0-9]*.\([0-9]*\)"`
    SCILAB_VERSION_MICRO=`expr match "$SCILAB_VERSION" " *[0-9]*.[0-9]*.\([0-9]*\)"`

    if test $SCILAB_VERSION_MAJOR -lt $REQUIRED_SCILAB_MAJOR
    then
      AC_MSG_ERROR([scilab major version does not match])
    else
      if test $SCILAB_VERSION_MINOR -lt $REQUIRED_SCILAB_MINOR
      then
        AC_MSG_ERROR([scilab minor version does not match])
      else
        if test $SCILAB_VERSION_MICRO -lt $REQUIRED_SCILAB_MICRO
        then
          AC_MSG_ERROR([scilab micro version does not match])
        fi
      fi
    fi
  fi

  if test "xwith_scilab_toolbox_dir" != "xyes"
  then
    SCILAB_TOOLBOX_DIR="$with_scilab_toolbox_dir"
  else
    SCILAB_TOOLBOX_DIR="$SCILAB_DIR/contrib/$PACKAGE_NAME-$PACKAGE_VERSION"
  fi

  AC_SUBST(HAVE_SCILAB)
  AC_SUBST(SCILAB_EXE)
  AC_SUBST(SCILAB_DIR)
  AC_SUBST(SCILAB_TOOLBOX_DIR)
  AC_SUBST(SCILAB_VERSION_MAJOR)
  AC_SUBST(SCILAB_VERSION_MINOR)
  AC_SUBST(SCILAB_VERSION_MICRO)
  AC_SUBST(SCILAB_VERSION)
])
