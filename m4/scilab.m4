AC_DEFUN([AC_CHECK_SCILAB],
[
  AH_TEMPLATE(HAVE_SCILAB,
   [Defined to 1 if Scilab is present on the system])

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

  AC_ARG_WITH(toolbox_install_prefix,
		AC_HELP_STRING([--with-toolbox-install-prefix=DIR],[Set the path to the toolbox installation directory]),
		[with_toolbox_install_prefix=$withval],
		[with_toolbox_install_prefix='yes']
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
    SCIEXE="$with_scilab_prefix/bin/scilab"
    AC_DEFINE(HAVE_SCILAB)
  else
    AC_CHECK_PROG([has_scilab],[scilab],[yes],[no], , )
    if test x$has_scilab = xno; then
      AC_MSG_ERROR([[Scilab binary program was found in your PATH], your PATH is $PATH])
    fi
    SCIEXE="scilab"
    AC_DEFINE(HAVE_SCILAB)
  fi

  dnl test if SCI is defined. If not, test if scilab is in path,
  dnl then get SCI from scilab -nwni.
  if test -z "$SCI"; then
    cmd='F=mopen("getpath.incl","w");
         mfprintf(F,SCI);
         mclose(F);exit;'
    echo "$cmd" > getpath.sci
    $SCIEXE -nw -f getpath.sci >/dev/null
    SCIDIR=`cat getpath.incl`
    rm -f getpath.sci getpath.incl 
  else
    SCIDIR="$SCI"
  fi

  dnl get scilab version
  cmd='F=mopen("version.incl","w");
       ver=getversion();
       mfprintf(F,ver);
       mclose(F);exit;'
  echo "$cmd" > version.sci
  $SCIEXE -nwni -f version.sci >/dev/null
  SCIVERSION=`cat version.incl`
  rm -f version.sci version.incl 

  scilab_tmp_version=`expr match $SCIVERSION '.*\(branch\).*'`

  if test "x$scilab_tmp_version" = "xbranch"
  then
    SCILAB_VERSION_MAJOR=-1
    SCILAB_VERSION_MINOR=-1
    SCILAB_VERSION_MICRO=-1
  else
    SCILAB_VERSION_MAJOR=`expr match "$SCIVERSION" " *\([0-9]*\)"`
    SCILAB_VERSION_MINOR=`expr match "$SCIVERSION" " *[0-9]*.\([0-9]*\)"`
    SCILAB_VERSION_MICRO=`expr match "$SCIVERSION" " *[0-9]*.[0-9]*.\([0-9]*\)"`

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

  if test "xwith_toolbox_install_prefix" != "xyes"
  then
    TOOLBOX_INSTALL_DIR="$with_toolbox_install_prefix"
  else
    TOOLBOX_INSTALL_DIR="$SCIDIR/contrib/$PACKAGE_NAME-$PACKAGE_VERSION"
  fi

  AC_SUBST(HAVE_SCILAB)
  AC_SUBST(SCIEXE)
  AC_SUBST(SCIDIR)
  AC_SUBST(TOOLBOX_INSTALL_DIR)
  AC_SUBST(SCILAB_VERSION_MAJOR)
  AC_SUBST(SCILAB_VERSION_MINOR)
  AC_SUBST(SCILAB_VERSION_MICRO)
])
